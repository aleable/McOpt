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

import numpy as np
import networkx as nx
import pickle

#######################################


def data_analysis(metric_mode, output_path, g, nnode, comm_list, length, opttdens, optpot, pflux, g_opt, flux_mat_opt, length_opt):
    """switcher metric to export"""

    nnode = int(len(g.nodes))

    def continuefunc():
        print("metrics: no export")
        return 0

    switcher = {
        "0": lambda: continuefunc(),
        "1": lambda: graph_visual(output_path, g, nnode, length, comm_list, opttdens, optpot, pflux, g_opt, flux_mat_opt, length_opt),
        "2": lambda: cost(output_path, g, nnode, length, opttdens, optpot, pflux, flux_mat_opt, length_opt),
        "3": lambda: loops(output_path, g, g_opt, pflux),
        "4": lambda: idle_edges(output_path, pflux)
    }

    return switcher.get(metric_mode, lambda: print("ERROR: invalid flag (export)"))()


def graph_visual(output_path, g, nnode, length, comm_list, opttdens, optpot, pflux, g_opt, flux_mat_opt, length_opt):
    """necessary quantities for graph visualization"""

    print("metrics: graph visualization")

    inc_mat = nx.incidence_matrix(g, nodelist=list(range(nnode)), oriented=True)
    inv_len_mat = np.diag(1/length)
    mu_mat = np.diag(opttdens)

    flux_mat = np.matmul(mu_mat*inv_len_mat*np.transpose(inc_mat), optpot)
    flux_norm = np.linalg.norm(flux_mat, axis=1)

    flux_norm_opt = np.linalg.norm(flux_mat_opt, axis=1)

    pickle.dump(flux_norm, open(output_path + "flux_norm_dyn.pkl", "wb"))
    pickle.dump(flux_norm_opt, open(output_path + "flux_norm_opt.pkl", "wb"))
    pickle.dump(comm_list, open(output_path + "comm_list.pkl", "wb"))
    pickle.dump(length, open(output_path + "length_edges_dyn.pkl", "wb"))
    pickle.dump(length_opt, open(output_path + "length_edges_opt.pkl", "wb"))
    nx.write_gpickle(g, output_path + "graph_dynamics_dyn.gpickle")
    nx.write_gpickle(g_opt, output_path + "graph_dynamics_opt.gpickle")
    
    # print('Output files: \n')
    # print(output_path + "flux_norm_dyn.pkl")
    # print(output_path +"flux_norm_opt.pkl")
    # print(output_path +"comm_list.pkl")
    # print(output_path +"length_edges_dyn.pkl")
    # print(output_path +"length_edges_opt.pkl")
    # print(output_path +"graph_dynamics_opt.gpickle")
    # print(output_path +"flux_norm_dyn.pkl")

    print('Output folder: ', output_path)
    print('Output files:')
    print("comm_list.pkl")
    print("flux_norm_dyn.pkl")
    print("flux_norm_opt.pkl")
    print("length_edges_dyn.pkl")
    print("length_edges_opt.pkl")
    print("graph_dynamics_dyn.gpickle")
    print("graph_dynamics_opt.gpickle")
 


def cost(output_path, g, nnode, length, opttdens, optpot, pflux, flux_mat_opt, length_opt):
    """computing J_gamma"""

    inc_mat = nx.incidence_matrix(g, nodelist=list(range(nnode)), oriented=True)
    inv_len_mat = np.diag(1/length)
    mu_mat = np.diag(opttdens)

    flux_mat = np.matmul(mu_mat*inv_len_mat*np.transpose(inc_mat), optpot)
    flux_norm = np.linalg.norm(flux_mat, axis=1)
    flux_norm_opt = np.linalg.norm(np.array(flux_mat_opt), axis=1)

    cost = np.sum(length*(flux_norm**(2*(2-pflux)/(3-pflux))))
    cost_opt = np.sum(length_opt*(flux_norm_opt**(2*(2-pflux)/(3-pflux))))

    with open(output_path + "cost_dyn.dat", "a") as loops_file:
        loops_file.write(str(pflux) + " " + str(cost) + "\n")
    with open(output_path + "cost_opt.dat", "a") as loops_opt_file:
        loops_opt_file.write(str(pflux) + " " + str(cost_opt) + "\n")


def loops(output_path, g, g_opt, pflux):
    """computing number of basis loops"""

    print("metrics: exporting number of basis loops")

    cycle = len(nx.cycle_basis(g))
    cycle_opt = len(nx.cycle_basis(g_opt))

    with open(output_path + "cycles_dyn.dat", "a") as cost_file:
        cost_file.write(str(pflux) + " " + str(cycle) + "\n")
    with open(output_path + "cycles_opt.dat", "a") as cost_opt_file:
        cost_opt_file.write(str(pflux) + " " + str(cycle_opt) + "\n")


def idle_edges(output_path, pflux):
    """computing number of idle edges"""

    print("metrics: exporting idle edges")

    with open(output_path + "index_to_remove_dyn.pkl", "rb") as idle_file:
        idle_index = pickle.load(idle_file)
    with open(output_path + "index_to_remove_opt.pkl", "rb") as idle_opt_file:
        idle_opt_index = pickle.load(idle_opt_file)

    with open(output_path + "idle_dyn.dat", "a") as idle_file:
        idle_file.write(str(pflux) + " " + str(len(idle_index)) + "\n")
    with open(output_path + "idle_opt.dat", "a") as idle_opt_file:
        idle_opt_file.write(str(pflux) + " " + str(len(idle_opt_index)) + "\n")
