"""
McOpt (Multicommodity Optimal Transport) -- https://github.com/aleable/McOpt

Contributors:
    Alessandro Lonardi
    Enrico Facca
    Caterina De Bacco
"""

import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import identity
from scipy.sparse.linalg import spsolve


def tdensinit(self):
    """Initialization of the conductivities: mu_e ~ U(0,1)

    Returns:
        self.tdens: np.array, initialized conductivities
    """

    prng = np.random.RandomState(seed=self.seed)
    self.tdens = np.array([prng.uniform(0, 1) for i in range(self.g.number_of_edges())])

    return self.tdens


def dyn(self):
    """Execute dynamics

    Returns:
        self.tdens: np.array, conductivities at convergence
        pot: np.array, potentials at convergence
        cost_stack: np.array, cost at each iteration
    """

    print("* running [dynamics]")

    ####################################################################
    # INITIALIZATION
    ####################################################################

    time_iteration = 0
    # only needed if spsolve has problems (inside update)
    prng = np.random.RandomState(seed=self.seed)

    # initialization quantities of dynamics
    nnode = self.g.number_of_nodes()
    # incidence matrix
    inc_mat = csr_matrix(nx.incidence_matrix(self.g, nodelist=list(range(nnode)), oriented=True))
    inc_transpose = csr_matrix(inc_mat.transpose())
    inv_len_mat = diags(1/self.length, 0)
    td_mat = diags(self.tdens, 0)
    # network weighted Laplacian
    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose

    # spsolve
    stiff_relax = stiff + self.relax_linsys * identity(nnode)  # avoid zero kernel
    pot = spsolve(stiff_relax, self.forcing)                   # pressure vector

    # executing dynamics
    convergence_achieved = False
    cost = 0
    cost_stack = []

    td_mat = np.diag(self.tdens)

    ####################################################################

    while not convergence_achieved and time_iteration <= self.tot_time:

        ####################################################################
        # RUNNING THE DYNAMICS
        ####################################################################

        time_iteration += 1

        # update tdens-pot system
        tdens_old = self.tdens
        pot_old = pot

        # equations update
        self.tdens, pot, info = update(self, pot, inc_mat, inc_transpose, inv_len_mat)

        # singular stiffness matrix
        if info != 0:
            self.tdens = tdens_old + prng.rand(*tdens_old.shape) * np.mean(tdens_old) / 1000.
            pot = pot_old + prng.rand(*pot.shape) * np.mean(pot_old) / 1000.

        # convergence criteria
        var_tdens = max(np.abs(self.tdens - tdens_old)) / self.time_step
        convergence_achieved, cost, abs_diff_cost = cost_convergence(self,
                                                                     pot,
                                                                     time_iteration,
                                                                     cost,
                                                                     var_tdens,
                                                                     inc_mat,
                                                                     inv_len_mat,
                                                                     convergence_achieved)

        cost_stack.append(cost)
        if self.verbose and time_iteration % 10 == 0:
            print('\tit=%3d, err=%5.8f, J_diff=%5.8e' % (time_iteration, var_tdens, abs_diff_cost))

        elif time_iteration >= self.tot_time:
            convergence_achieved = True
            self.tdens = tdens_old
            print("\tERROR: dyn NOT converged [iteration > maxit]")

    if convergence_achieved:
        print("\tconvergence achieved [dynamics]")
        cost_stack = np.array(cost_stack)
        return self.tdens, pot, cost_stack
    else:
        print("\tERROR: dyn NOT converged")

        ####################################################################


def update(self, pot, inc_mat, inc_transpose, inv_len_mat):
    """One step update

    Parameters:
        pot: np.array, potential matrix on nodes
        inc_mat: sparse.matrix, oriented incidence matrix
        inc_transpose: sparse.matrix, oriented incidence matrix transposed
        inv_len_mat: sparse.matrix, diagonal matrix 1/l_e

    Returns:
        tdens: np.array, updated conductivities
        pot: np.array, updated potential matrix on nodes
        info: bool, sanity check flag spsolve
        """

    # update mu
    grad = inv_len_mat * inc_transpose * pot
    if self.coupling == "l2":
        rhs_ode = (self.tdens ** self.pflux) * ((grad ** 2).sum(axis=1)) - self.tdens
    if self.coupling == "l1":
        rhs_ode = (self.tdens ** self.pflux) * (np.abs(grad).sum(axis=1)) ** 2 - self.tdens

    self.tdens = self.tdens + self.time_step*rhs_ode
    td_mat = diags(self.tdens, 0)
    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose

    # spsolve
    stiff_relax = stiff + self.relax_linsys * identity(self.g.number_of_nodes())
    pot = spsolve(stiff_relax, self.forcing)

    # sanity check
    if np.any(np.isnan(pot)):
        info = -1
        pass
    else:
        info = 0

    return self.tdens, pot, info


def cost_convergence(self, pot, time_iteration, cost, var_tdens, inc_mat, inv_len_mat, convergence_achieved):
    """Evaluating convergence

    Parameters:
        pot: np.array, potential matrix on nodes
        time_iteration: int, iteration number
        cost: float, cost
        var_tdens: time step difference conductivities
        inc_mat: sparse.matrix, oriented incidence matrix
        inv_len_mat: sparse.matrix, diagonal matrix 1/l_e
        convergence_achieved: bool, convergence flag

    Returns:
        convergence_achieved: bool, updated convergence flag
        cost_update: float, updated cost
        abs_diff_cost: float, difference cost

    """

    td_mat = np.diag(self.tdens)
    flux_mat = np.matmul(td_mat * inv_len_mat * np.transpose(inc_mat), pot)

    if self.coupling == "l2":
        flux_norm = np.linalg.norm(flux_mat, axis=1)**2
    if self.coupling == "l1":
        flux_norm = np.linalg.norm(flux_mat, axis=1, ord=1)**2

    cost_update = np.sum(self.length * (flux_norm ** ((2 - self.pflux) / (3 - self.pflux))))
    abs_diff_cost = abs(cost_update - cost)

    # only the cost is calculated to evaluate convergence if beta >= 1.0
    if self.pflux >= 1.0:
        if abs_diff_cost < self.tau_cost_dyn and time_iteration > 30:
            convergence_achieved = True

    else:
        if abs_diff_cost < self.tau_cost_dyn and var_tdens < self.tau_cond_dyn:
            convergence_achieved = True

    return convergence_achieved, cost_update, abs_diff_cost
