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

import os
import pickle
import numpy as np
import networkx as nx
from argparse import ArgumentParser
from initialization import topology_generation
from initialization import file2graph
from initialization import coord_generation
from initialization import eucledian_bias
from initialization import rhs_generation
from initialization import file2forcing
from dynamics import tdensinit
from dynamics import dyn
from dynamics import abs_trimming_dyn
from optimization import opt
from optimization import abs_trimming_opt
from export import data_analysis

#######################################
# FLAGS
#######################################

p = ArgumentParser()

# flags to create or read files with different names
p.add_argument("-fg", "--file_graph", type=str, default="graph.dat")    # graph parameters (change to "graph_generated.dat" if -topol = 0)
p.add_argument("-fa", "--file_adj", type=str, default="adj.dat")        # adjacency list file name (change to "adj_generated.dat" if -topol = 0)
p.add_argument("-fs", "--file_rhs", type=str, default="rhs.dat")        # rhs file name (change to "rhs_generated.dat" if -topol = 0)
p.add_argument("-fc", "--file_coord", type=str, default="coord.dat")    # coordinate file name (change to "coord_generated.dat" if -topol = 0)

p.add_argument("-topol", "--topol_gen", type=str, default="1")          # 0 = generate topology (Waxman graph), 1 = import topology
p.add_argument("-N", "--n_node", type=int, default=30)                  # number of nodes in Waxman graph, ignored if -topol = 1
p.add_argument("-S", "--n_sosi", type=int, default=10)                  # number of commodities in Waxman graph (ignored if -topol = 1)
p.add_argument("-mass", "--mass_gen", type=str, default="1")            # 0 = generate {g^i} mass file, 1 = import mass file (ignored if -a = im)
p.add_argument("-coord", "--coord_gen", type=str, default="1")          # 0 = generate coordinate file, 1 = import file (default = 1 if Waxman graph is generated)

p.add_argument("-a", "--assignation_mass", type=str, default="ia")
# fr = Fake Receiver
# ia = Influence Assignment
# im = importing forcing from previous run

p.add_argument("-M", "--tot_mass", type=int, default=1000)              # total mass to distribute with Fake Receivers or Influence Assignment (ignored if -mass = 1 or -a = im)
p.add_argument("-l", "--bias_eucl", type=str, default="eucl")           # bias =  1 + 1.0e-3 * U(0,1), eucl = eucledian distance

p.add_argument("-b", "--exponent", type=float, default=0.5)             # beta
p.add_argument("-t", "--tau_trimming", type=float, default=1.0e-3)      # trimming threshold

p.add_argument("-d", "--metric", type=str, default="0")                 # metric to export for analysis
# 0 = no export
# 1 = graph visualization
# 2 = J_gamma
# 3 = number of basis loops
# 4 = number of idle edges
p.add_argument("-v", "--verbose", type=int, default=0)                  # whether to show comments or not

# python main.py -fg "graph_generated.dat" -fa "adj_generated.dat" -fs "rhs_generated.dat" -fc "coord_generated.dat" -topol 0 -N 200 -S 200 -a fr -mass 0 -M 10000 -b 1.5 -d 1
args = p.parse_args()

#######################################
# MAIN
#######################################

#######################################
# initialization
#######################################

path = os.getcwd()
input_path = os.getcwd() + "/../data/input/"
output_path = os.getcwd() + "/../data/output/"

# create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# generating or importing topology
nnode = args.n_node
ncomm = args.n_sosi
topol_mode = args.topol_gen

topology_generation(topol_mode, input_path, nnode, ncomm)   # generating

graph_file_name = input_path + args.file_graph
adj_file_name = input_path + args.file_adj
nnode, ncomm, nedge, graph = file2graph(graph_file_name, adj_file_name)     # importing

# assign node coordinates and choose eucledian/uniform edge lengths
coord_mode = args.coord_gen
length_mode = args.bias_eucl
coord_file_name = input_path + args.file_coord

if topol_mode == "0":   # coordinates are pre-assigned by the Waxman model
    coord_mode = "1"

coord_generation(coord_mode, input_path, coord_file_name, graph, nnode)

length = np.zeros(len(graph.edges()))
length = eucledian_bias(length_mode, graph, length)

# generate/import rhs
rhs_mode = args.mass_gen
tot_mass = args.tot_mass

rhs_generation(rhs_mode, input_path, graph, ncomm, tot_mass)

# generate rhs: Fake Receviver or Influence Assignment
assignation_mode = args.assignation_mass
rhs_file_name = input_path + args.file_rhs

forcing, comm_list, transit_list = file2forcing(assignation_mode, rhs_file_name, input_path, graph, nnode, ncomm)
forcing = np.array(forcing)

# pickling in case same forcing wants to be used again
pickle.dump(comm_list, open(input_path + "comm_list.pkl", "wb"))
pickle.dump(forcing, open(input_path + "forcing.pkl", "wb"))

#######################################
# dynamics and optimization
#######################################

pflux = args.exponent   # beta
tol_var_tdens = 10e-3   # tau_dyn (relax por change convergence criteria for beta > 1.45)

# running dynamics
tdens_0 = tdensinit(nedge)  # conductivities initialization
opttdens, optpot = dyn(graph, tdens_0, pflux, length, forcing, tol_var_tdens, comm_list , verbose = args.verbose)

# running optimization
flux_mat_opt = opt(graph, tdens_0, pflux, length, forcing, verbose = args.verbose)

#######################################
# trimming
#######################################

tau = args.tau_trimming

graph_opt, length_opt, flux_mat_opt = abs_trimming_opt(graph, flux_mat_opt, length, output_path, tau)
graph, opttdens, length = abs_trimming_dyn(graph, opttdens, optpot, length, output_path, tau)

print("connected components dyn:", nx.number_connected_components(graph))
print("connected components opt:", nx.number_connected_components(graph_opt))

#######################################
# exporting metrics for data analysis
#######################################

metric_mode = args.metric

data_analysis(metric_mode, output_path, graph, nnode, comm_list, length, opttdens, optpot, pflux, graph_opt, flux_mat_opt, length_opt)
