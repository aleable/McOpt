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


def opt(self):
    """Execute fixed-point iteration

    Returns:
        flux: np.array, fluxes at convergence
        cost_stack: np.array, cost at each iteration
    """

    print("* running [fixed-point]")

    ####################################################################
    # INITIALIZATION
    ####################################################################

    time_iteration = 0

    # initialization quantities of dynamics
    nnode = self.g.number_of_nodes()
    inc_mat = csr_matrix(nx.incidence_matrix(self.g, nodelist=list(range(nnode)), oriented=True))   # incidence matrix
    inc_transpose = csr_matrix(inc_mat.transpose())
    inv_len_mat = diags(1/self.length, 0)
    td_mat = diags(self.tdens, 0)
    # network weighted Laplacian
    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose

    # spsolve
    stiff_relax = stiff + self.relax_linsys * identity(nnode)  # avoid zero kernel
    pot = spsolve(stiff_relax, self.forcing)                   # pressure vector
    flux_mat = td_mat * inv_len_mat * inc_transpose * pot      # flux initialization

    # executing dynamics
    convergence_achieved = False
    cost = 0
    cost_stack = []

    ####################################################################

    while not convergence_achieved and time_iteration <= self.tot_time:

        ####################################################################
        # RUNNING FIXED-POINT ITERATION
        ####################################################################

        time_iteration += 1

        flux_mat_new = update_iterative(self, inc_mat, inc_transpose, inv_len_mat, flux_mat)

        # convergence criteria
        norm_old = np.linalg.norm(flux_mat, axis=1)
        norm_new = np.linalg.norm(flux_mat_new, axis=1)
        abs_diff = max(np.abs(norm_new - norm_old))
        convergence_achieved, cost, abs_diff_cost = cost_convergence(self,
                                                                     time_iteration,
                                                                     cost,
                                                                     flux_mat,
                                                                     abs_diff,
                                                                     convergence_achieved)

        cost_stack.append(cost)
        if self.verbose and time_iteration % 10 == 0:
            print('\tit=%3d, err=%5.8f, J_diff=%5.8e' % (time_iteration, abs_diff, abs_diff_cost))

        # convergence achieved
        if time_iteration >= self.tot_time:
            convergence_achieved = True
            flux_mat = flux_mat_new
            print("\tERROR: fixed-point NOT converged [iteration > maxit]")

        elif time_iteration < self.tot_time:
            flux_mat = flux_mat_new

    if convergence_achieved:
        print("\tconvergence achieved [fixed-point]")
        cost_stack = np.array(cost_stack)
        return flux_mat, cost_stack
    else:
        print("\tERROR: fixed-point NOT converged")

        ####################################################################


def update_iterative(self, inc_mat, inc_transpose, inv_len_mat, flux_mat):
    """One step update

    Parameters:
        inc_mat: sparse.matrix, oriented incidence matrix
        inc_transpose: sparse.matrix, oriented incidence matrix transposed
        inv_len_mat: sparse.matrix, diagonal matrix 1/l_e
        flux_mat: np.array, fluxes

    Returns:
        flux: np.array, updated fluxes
        """

    if self.coupling == "l2":
        flux_norm = np.linalg.norm(flux_mat, axis=1)**2
    if self.coupling == "l1":
        flux_norm = np.linalg.norm(flux_mat, axis=1, ord=1)**2

    # computing scaling ad updating conductivities
    temp = (np.sum(self.length*flux_norm**((2-self.pflux)/(3-self.pflux))))**(1/(2-self.pflux))
    self.tdens = (1/temp)*flux_norm**(1/(3-self.pflux))
    td_mat = diags(self.tdens, 0)

    # computing fluxes
    temp_pinv = np.linalg.pinv(csr_matrix.todense(inc_mat*td_mat*inv_len_mat*inc_transpose))
    lagrange_mult = temp_pinv*self.forcing
    flux = td_mat*inv_len_mat*inc_transpose*lagrange_mult

    return flux


def cost_convergence(self, time_iteration, cost, flux_mat, abs_diff, convergence_achieved):
    """Evaluating convergence

    Parameters:
        time_iteration: int, iteration number
        cost: float, old cost
        flux_mat: np.array, fluxes matrix
        abs_diff: max flux fluxes difference
        convergence_achieved: bool, convergence flag

    Returns:
        convergence_achieved: bool, updated convergence flag
        cost_update: float, updated cost
        abs_diff_cost: float, difference cost

    """

    if self.coupling == "l2":
        flux_norm = np.linalg.norm(flux_mat, axis=1)**2
    if self.coupling == "l1":
        flux_norm = np.linalg.norm(flux_mat, axis=1, ord=1)**2

    cost_update = np.sum(self.length * (flux_norm ** ((2 - self.pflux) / (3 - self.pflux))))
    abs_diff_cost = abs(cost_update - cost)

    if self.pflux >= 1.0:  # only cost for beta >= 1.0
        if abs_diff_cost < self.tau_cost_opt and time_iteration > 30:
            convergence_achieved = True
    else:
        if abs_diff_cost < self.tau_cost_opt and abs_diff < self.tau_cond_opt:
            convergence_achieved = True

    return convergence_achieved, cost_update, abs_diff_cost
