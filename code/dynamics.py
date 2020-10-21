"""
McOpt (Multi-commodity Optimal Transport)

Alessandro Lonardi
Enrico Facca
Caterina De Bacco

root_file: main.py
branch_file: -initialization.py
             -dynamics.py
             -optimization.py
             -export.py
"""

#######################################
# PACKAGES
#######################################

import pickle,time, warnings
import numpy as np
import networkx as nx
import random
import copy
import scipy as sp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import identity
from scipy.sparse.linalg import spsolve

#######################################

warnings.filterwarnings("ignore", message="Matrix is exactly singular")

def tdensinit(nedge, seed=10):
    """initialization of the conductivities: mu_e ~ U(0,1)"""
    prng = np.random.RandomState(seed=seed)
    tdens_0 = [prng.uniform(0, 1) for i in range(nedge)]
    tdens_0 = np.array(tdens_0)

    return tdens_0


def dyn(g, tdens_0, pflux, length, forcing, tol_var_tdens, comm_list,seed=10,verbose=False):
    """dynamics method"""

    print("\ndynamics...")

    relax_linsys = 1.0e-5       # relaxation for stiffness matrix
    tot_time = 1000            # upper bound on number of time steps
    time_step = 0.5             # delta_t (reduce as beta gets close to 2)
    time_iteration = 0
    threshold_cost = 1.0e-6   # threshold for stopping criteria using cost
    prng = np.random.RandomState(seed=seed) # only needed if spsolve has problems (inside update)

    nnode = g.number_of_nodes()

    # initialization quantities of dynamics
    inc_mat = csr_matrix(nx.incidence_matrix(g, nodelist=list(range(nnode)), oriented=True))  # B
    inc_transpose = csr_matrix(inc_mat.transpose())     # B^T
    inv_len_mat = diags(1/length, 0)    # diag[1/l_e]
    tdens = tdens_0
    td_mat = diags(tdens, 0)    # matrix M
    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose  # B diag[mu] diag[1/l_e] B^T

    # spsolve
    stiff_relax = stiff + relax_linsys * identity(nnode)    # avoid zero kernel
    pot = spsolve(stiff_relax, forcing.transpose())     # pressure

    # dynamics
    convergence_achieved = False
    cost = 0

    fmax = forcing.max()

    while not convergence_achieved and time_iteration < tot_time:

        time_iteration += 1

        # update tdens-pot system
        tdens_old = tdens
        pot_old = pot

        # equations update
        tdens, pot, grad, info = update(tdens, pot, inc_mat, inc_transpose, inv_len_mat, forcing, time_step, pflux, relax_linsys, nnode)

        # singular stiffness matrix
        if info != 0:
            tdens = tdens_old + prng.rand(*tdens.shape) * np.mean(tdens_old)/1000.
            pot = pot_old + prng.rand(*pot.shape) * np.mean(pot_old)/1000.

        # 1) convergence with conductivities
        # var_tdens = max(np.abs(tdens - tdens_old))/time_step
        # print(time_iteration, var_tdens)

        # 2) an alternative convergence criteria: using total cost and maximum variation of conductivities
        var_tdens = max(np.abs(tdens - tdens_old))/time_step
        convergence_achieved, cost, abs_diff_cost = cost_convergence(threshold_cost, cost, tdens, pot, inc_mat, inv_len_mat, length, pflux, convergence_achieved, var_tdens)
        if verbose: 
            # print(time_iteration, var_tdens/forcing.max(), abs_diff_cost)

            # print('\r','It=',it,'err=', abs_diff,'J-J_old=',abs_diff_cost,sep=' ', end='', flush=True)
            print('\r','it=%3d, err/max_f=%5.2f, J_diff=%8.2e' % (time_iteration,var_tdens/fmax,abs_diff_cost),sep=' ', end=' ', flush=True)
            time.sleep(0.05)

        if var_tdens < tol_var_tdens:
            convergence_achieved = True

        elif time_iteration >= tot_time:
            convergence_achieved = True
            tdens = tdens_old
            print("ERROR: convergence dynamics not achieved, iteration time > maxit")

    if convergence_achieved:
        return tdens, pot
    else:
        print("ERROR: convergence dynamics not achieved")


def update(tdens, pot, inc_mat, inc_transpose, inv_len_mat, forcing, time_step, pflux, relax_linsys, nnode):
    """dynamics update"""

    grad = inv_len_mat*inc_transpose*pot    # discrete gradient
    rhs_ode = (tdens**pflux)*((grad**2).sum(axis=1)) - tdens

    # update conductivity
    tdens = tdens + time_step*rhs_ode
    td_mat = sp.sparse.diags(tdens, 0)
    # update stiffness matrix
    stiff = inc_mat*td_mat*inv_len_mat*inc_transpose

    # spsolve
    stiff_relax = stiff + relax_linsys*identity(nnode)  # avoid zero kernel

    # update potential
    
    pot = spsolve(stiff_relax, forcing.transpose())   # pressure
    if np.any(np.isnan(pot)): # or np.any(pot != pot)
        info = -1
        pass
    else:
        info = 0

    return tdens, pot, grad, info

def cost_convergence(threshold_cost, cost, tdens, pot, inc_mat, inv_len_mat, length, pflux, convergence_achieved, var_tdens):
    """computing convergence using total cost: setting a high value for maximum conducivity variability"""

    td_mat = np.diag(tdens)
    flux_mat = np.matmul(td_mat*inv_len_mat*np.transpose(inc_mat), pot)
    flux_norm = np.linalg.norm(flux_mat, axis=1)
    cost_update = np.sum(length*(flux_norm**(2*(2-pflux)/(3-pflux))))
    abs_diff_cost = abs(cost_update - cost)

    if pflux > 1.45:
        if abs_diff_cost < threshold_cost :
            convergence_achieved = True
        
    else: 
        if abs_diff_cost < threshold_cost and var_tdens < 1:
            convergence_achieved = True

    return convergence_achieved, cost_update, abs_diff_cost

def abs_trimming_dyn(g, opttdens, optpot, length, output_path, tau):
    """obtaining the trimmed graph using an absolute threshold"""

    print("trimming dynamics graph...\n")

    nnode = len(g.nodes())
    inc_mat = csr_matrix(nx.incidence_matrix(g, nodelist=list(range(nnode)), oriented=True))    # delta
    inv_len_mat = diags(1/length, 0)    # diag[1/l_e]

    mu_mat = np.diag(opttdens)
    flux_mat = np.matmul(mu_mat*inv_len_mat*np.transpose(inc_mat), optpot)
    flux_norm = np.linalg.norm(flux_mat, axis=1)

    edges_list = list(g.edges)

    g_trimmed = copy.deepcopy(g)
    g_final = copy.deepcopy(g)

    index_to_remove = []

    # remove flux below trimming threshold
    for i in range(len(flux_norm)):
        g_trimmed.remove_edge(edges_list[i][0], edges_list[i][1])
        if flux_norm[i] < tau:
            index_to_remove.append(i)
            # iteratively trim the edges
            node_1 = edges_list[i][0]
            node_2 = edges_list[i][1]
            g_final.remove_edge(node_1, node_2)

    opttdens_final = np.delete(opttdens, index_to_remove)
    length_final = np.delete(length, index_to_remove)

    temp = {}
    i = 0
    for e in g_final.edges:
        temp[(e[0], e[1])] = i
        i += 1

    # dumping for graph visualization
    pickle.dump(index_to_remove, open(output_path + "index_to_remove_dyn.pkl", "wb"))

    return g_final, opttdens_final, length_final
