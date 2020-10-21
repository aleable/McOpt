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

import pickle, time
import numpy as np
import networkx as nx
import copy
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import identity
from scipy.sparse.linalg import spsolve

#######################################


def opt(g, tdens_0, pflux, length, forcing,verbose=False):
    """optimization method"""

    print("\noptimization...")

    relax_linsys = 1.0e-5   # relaxation for stiffness matrix
    conv_thresh = 1.0e-5    # convergence threshold
    maxit = 1000
    threshold_cost = 1.0e-6 # threshold for stopping criteria using cost

    nnode = g.number_of_nodes()

    tdens = tdens_0
    td_mat = diags(tdens, 0)  # matrix M initialization

    inc_mat = csr_matrix(nx.incidence_matrix(g, nodelist=list(range(nnode)), oriented=True))    # B
    inc_transpose = csr_matrix(inc_mat.transpose())     # delta^T
    inv_len_mat = diags(1 / length, 0)  # diag[1/l_e]

    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose  # B diag[mu] diag[1/l_e] B^T
    stiff_relax = stiff + relax_linsys * identity(nnode)    # avoid zero kernel
    pot = spsolve(stiff_relax, forcing.transpose())         # pressure

    flux_mat = td_mat * inv_len_mat * inc_transpose * pot  # flux initialization

    convergence_achieved = False  # flag for convergence

    it = 0 # iteration time, to check for convergence
    cost = 0

    while not convergence_achieved:
        flux_mat_new = update_iterative(nnode, relax_linsys, inc_mat, inc_transpose, length, inv_len_mat, forcing, flux_mat, pflux)

        # 1) converence using fluxes
        # norm_old = np.linalg.norm(flux_mat, axis=1)
        # norm_new = np.linalg.norm(flux_mat_new, axis=1)
        # abs_diff = max(np.abs(norm_new - norm_old))
        # if verbose: print(it, abs_diff)

        # 2) an alternative convergence criteria: using total cost and maximum variation of fluxes
        norm_old = np.linalg.norm(flux_mat, axis=1)
        norm_new = np.linalg.norm(flux_mat_new, axis=1)
        abs_diff = max(np.abs(norm_new - norm_old))
        convergence_achieved, cost, abs_diff_cost = cost_convergence(cost, norm_old, length, pflux, threshold_cost, convergence_achieved, abs_diff)
        if verbose: 
            # print('\r','It=',it,'err=', abs_diff,'J-J_old=',abs_diff_cost,sep=' ', end='', flush=True)
            print('\r','it=%3d, err=%8.2f, J_diff=%8.2e ' % (it,abs_diff,abs_diff_cost),sep=' ', end='', flush=True)
            time.sleep(0.05)

        # convergence achieved
        if it >= maxit:
            convergence_achieved = True
            flux_mat = flux_mat_new
            print("ERROR: convergence optimization not achieved, iteration time > maxit")
            return flux_mat

        elif abs_diff < conv_thresh:
            flux_mat = flux_mat_new
            convergence_achieved = True

        elif it < maxit:
            flux_mat = flux_mat_new

        it += 1

    if convergence_achieved:
        return flux_mat
    else:
        print("ERROR: convergence optimization not achieved")


def update_iterative(nnode, relax_linsys, inc_mat, inc_transpose, length, inv_len_mat, forcing, flux_mat, pflux):
    """optimization update"""

    flux_norm = np.linalg.norm(flux_mat, axis=1)

    temp = (np.sum(length*flux_norm**(2*(2-pflux)/(3-pflux))))**(1/(2-pflux))   # computing mu as function of F
    tdens = (1/temp)*flux_norm**(2/(3-pflux))
    td_mat = diags(tdens, 0)

    # 2) pinv
    temp_pinv = np.linalg.pinv(csr_matrix.todense(inc_mat*td_mat*inv_len_mat*inc_transpose))
    lagrange_mult = -temp_pinv*np.transpose(forcing)

    return -td_mat*inv_len_mat*inc_transpose*lagrange_mult

def cost_convergence(cost, flux_norm, length, pflux, threshold_cost, convergence_achieved, abs_diff):
    """computing convergence using total cost: setting a high value for maximum flux variability"""

    cost_update = np.sum(length*(flux_norm**(2*(2-pflux)/(3-pflux))))
    abs_diff_cost = abs(cost_update - cost)

    if abs_diff_cost < threshold_cost and abs_diff < 1:
        convergence_achieved = True

    return convergence_achieved, cost_update, abs_diff_cost

def abs_trimming_opt(g, flux_mat_opt, length, output_path, tau):
    """obtaining the trimmed optimized graph using an absolute threshold"""

    print("\n\ntrimming optimization graph...")

    flux_norm = np.linalg.norm(flux_mat_opt, axis=1)

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

    length_final = np.delete(length, index_to_remove)
    flux_mat_opt_final = np.delete(flux_mat_opt, index_to_remove, 0)

    temp = {}
    i = 0
    for e in g_final.edges:
        temp[(e[0], e[1])] = i
        i += 1

    # for graph visualization
    pickle.dump(index_to_remove, open(output_path + "index_to_remove_opt.pkl", "wb"))

    return g_final, length_final, flux_mat_opt_final
