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

import pickle
import numpy as np
import networkx as nx
import random
from scipy.spatial import distance

#######################################


def topology_generation(topol_mode, input_path, nnode, ncomm):
    """create or import topology file"""

    def continuefunc():
        print("topology file: imported")
        return 0

    switcher = {
        "0": lambda: generate_topol_file(input_path, nnode, ncomm),
        "1": lambda: continuefunc()
    }

    return switcher.get(topol_mode, lambda: print("ERROR: invalid flag (topology)"))()


def generate_topol_file(input_path, nnode, ncomm):
    """generate topology file"""

    print("topology file: generated")

    graph_file_path = input_path + "/graph_generated.dat"
    adjacency_file_path = input_path + "/adj_generated.dat"
    graph_coord_path = input_path + "/coord_generated.dat"

    # generate Waxman graph
    g = nx.waxman_graph(nnode, alpha=0.25, beta=0.25, L=1.0, domain=(0, 0, 1, 1))
    coord = {i: (g.nodes[i]["pos"][0], g.nodes[i]["pos"][1]) for i in range(nnode)}
    nedge = len(g.edges)

    # print graph topology in file
    with open(graph_file_path, "w") as f_topol:
        f_topol.write(str(nnode) + "\n")
        f_topol.write(str(ncomm) + "\n")
        f_topol.write(str(nedge) + "\n")

    # print adjacency list in file
    with open(adjacency_file_path, "w") as f_adj:
        i = 0
        for e in g.edges():
            f_adj.write(str(e[0]) + " " + str(e[1]) + "\n")
            i += 1

    # print node coordinates (automatically generated if Waxman model)
    with open(graph_coord_path, "w") as f_coord:
        for i in range(nnode):
            f_coord.write(str(coord[i][0]) + " " + str(coord[i][1]) + "\n")
            i += 1

    return 0


def file2graph(graph_file_name, adj_file_name):
    """generate graph from a topology file"""

    graph_file = open(graph_file_name, "r")
    lines = graph_file.readlines()

    nnode = int(lines[0][:-1])
    ncomm = int(lines[1][:-1])
    nedge = int(lines[2])

    adj_file = open(adj_file_name, "r")
    lines = adj_file.readlines()

    # generate graph adjacency
    topol = np.zeros([nedge,2], dtype=int) # (tail node, head node)
    iedge = 0
    for line in lines:
        topol[iedge][:] = [int(w) for w in line.split()[0:2]]
        iedge += 1

    g = nx.Graph()
    g.add_edges_from(topol)

    return nnode, ncomm, nedge, g


def eucledian_bias(length_mode, g, length):
    """length assignation: bias or eucledian"""

    # length type
    switcher = {
        "bias": lambda: bias(g, length),
        "eucl": lambda: eucledian(g, length)
    }

    return switcher.get(length_mode, lambda: print("ERROR: invalid flag (length)"))()


def eucledian(g, length):
    """eucledian lengths assigned to edges"""

    print("length: eucledian")

    i = 0
    for edge in g.edges():
        length[i] = distance.euclidean(g.nodes[edge[0]]["pos"], g.nodes[edge[1]]["pos"])
        i += 1

    return length


def bias(g, length):
    """fake constant lengths assigned to edges"""

    print("length: bias")

    for i in range(len(g.edges())):
        length[i] = 1 + 0.001 * random.uniform(0, 1)

    return length


def coord_generation(coord_mode, input_path, coord_file_name, g, nnode):
    """create or import coordinate file"""

    switcher = {
        "0": lambda: coordgenerating(input_path, g, nnode),
        "1": lambda: coordimporting(coord_file_name, g, nnode)
    }

    return switcher.get(coord_mode, lambda: print("ERROR: invalid flag (coordinates)"))()


def coordgenerating(input_path, g, nnode):
    """coordinates generated in square [0,1]x[0,1]"""

    print("coordinates: generated")

    # generating coordinates
    coord = np.zeros([nnode, 2])
    for inode in range(nnode):
        coord[inode][:] = [random.uniform(0, 1),random.uniform(0, 1)]
        inode += 1

    # assigning coord as attributes and printing in file
    for i in range(len(coord)):
        g.nodes[i]["pos"] = coord[i]

    coord_file = open(input_path + "/coord_generated.dat", "w")
    for i in range(int(nnode)):
        coord_file.write(str(coord[i][0]) + " " + str(coord[i][1]) + "\n")
    coord_file.close()

    return 0


def coordimporting(coord_file_name, g, nnode):
    """coordinates imported from file"""

    print("coordinates: imported")

    coord_file = open(coord_file_name, "r")
    input_lines = coord_file.readlines()

    coord = np.zeros([nnode, 2])
    inode = 0
    for line in input_lines:
        coord[inode][:] = [float(w) for w in line.split()[0:2]]
        inode += 1

    for i in range(nnode):
        g.nodes[i]["pos"] = (coord[i][0],coord[i][1])

    return 0


def rhs_generation(rhs_mode, input_path, g, ncomm, tot_mass):
    """mass generation mode: 0 = import, 1 = generate"""

    def continuefunc():
        print("rhs: imported")
        return 0

    switcher = {
        "0": lambda: rhsgeneration(input_path, g, ncomm, tot_mass),
        "1": lambda: continuefunc()
    }

    return switcher.get(rhs_mode, lambda : print("ERROR: invalid flag (rhs)"))()


def rhsgeneration(input_path, g, ncomm, tot_mass):
    """generate an artificial forcing file"""

    print("rhs: generated")

    #generate sources/sinks list
    comm_list = random.sample(g.nodes, ncomm)

    f_mass = np.zeros(int(ncomm))
    g_mass = np.zeros(int(ncomm))

    # generating g and assigning h to keep mass balance (h are equal to 0, can be modified in future implementation)
    for i in range(0, len(comm_list)):
        assigned_mass = random.randint(1, tot_mass)
        indexes = list(range(0, int(ncomm)))
        indexes.pop(i)
        g_mass[i] = assigned_mass

    # write rhs file
    rhs_file = open(input_path + "/rhs_generated.dat", "w")
    for i in range(int(ncomm)):
        rhs_file.write(str(comm_list[i]) + " " + str(g_mass[i]) + " " + str(f_mass[i]) + "\n")

    rhs_file.close()

    return 0


def file2forcing(assignation_mode, file_name, input_path, g, nnode, ncomm):
    """switcher mass generation mode: 0 = import, 1 = generate"""

    switcher = {
        "fr": lambda: file2forcing_fakereceiver(file_name, g, nnode, ncomm),
        "ia": lambda: file2forcing_impassign(file_name, g, nnode, ncomm),
        "im": lambda: forcing_importing(input_path, g)
    }

    return switcher.get(assignation_mode, lambda : print("ERROR: invalid flag (forcing assignment)"))()


def file2forcing_fakereceiver(file_name, g, nnode, ncomm):
    """generate forcing using Fake Receiver method: assigning each g^i to ONE random z_v^i"""

    print("forcing: Fake Receiver")

    # upload file and generate lists: icomm, g^i, h_i
    with open(file_name) as f:
        lines = f.readlines()
        temp_list = []
        for line in lines:
            temp_array = line.strip().split(" ")
            temp_list.append(np.array([int(float(element)) for element in temp_array]))
        temp_list = np.array(temp_list)

        comm_list = temp_list[:,0]  # list of sources and sinks
        temp_g = temp_list[:,1]     # list {g^i}i

    # creation of junction list
    transit_list = [node for node in g.nodes if node not in comm_list]

    # generation of the rhs matrix S
    forcing = np.zeros((ncomm,nnode)).tolist()

    # synthetic initialization of S with Fake Receivers method
    j = 0
    for i in comm_list:
        # in-flowing mass
        forcing[j][i] = float(temp_g[j])
        # out-flowing mass: fake receiver
        usable_nodes = list(set(list(range(int(nnode)))) - set(transit_list) - set([i]))
        chosen_index = random.choice(usable_nodes)
        forcing[j][chosen_index] = - float(temp_g[j])
        j += 1

    return forcing, comm_list, transit_list


def file2forcing_impassign(file_name, g, nnode, ncomm):
    """generate forcing using Influence Assignment method"""

    print("forcing: Influence Assignment")

    # upload file and generate lists: icomm, g^i, h_i
    with open(file_name) as f:
        lines = f.readlines()
        temp_list = []
        for line in lines:
            temp_array = line.strip().split(" ")
            temp_list.append(np.array([int(float(element)) for element in temp_array]))
        temp_list = np.array(temp_list)

        comm_list = temp_list[:,0]  # list of sources and sinks
        temp_g = temp_list[:,1]     # list {g^i}i

    # creation of junction list
    transit_list = [node for node in g.nodes if node not in comm_list]

    # generation of the rhs matrix S
    norm_g = []
    for i in range(ncomm):
        sum_g = np.sum(np.array(temp_g))
        sum_g -= temp_g[i]
        norm_g.append(sum_g)

    norm_mat = []
    for i in range(ncomm):
        norm_mat.append(np.array(temp_g)/norm_g[i])

    # synthetic initialization of S with Influence Assignment method
    mat_g = np.transpose(np.tile(temp_g, (ncomm, 1)))

    j = 0
    usable_nodes = []
    forcing = np.zeros((ncomm,nnode))
    forcing[:,comm_list] = - np.floor(np.multiply(mat_g, norm_mat))  # temp_g
    for i in comm_list:
        forcing[j][i] = float(temp_g[j])
        usable_nodes.append(list(set(list(range(int(nnode)))) - set(transit_list) - set([i])))
        j += 1

    residual = np.sum(np.array(forcing), axis=1)

    for j in range(len(comm_list)):
        chosen_index = random.choice(usable_nodes[j])
        forcing[j][chosen_index] -= residual[j]

    return forcing, comm_list, transit_list


def forcing_importing(input_path, g):
    """import forcing from last run"""

    print("forcing: imported")

    with open(input_path + "comm_list.pkl", "rb") as f_comm_list:
        comm_list = pickle.load(f_comm_list)
    with open(input_path + "forcing.pkl", "rb") as f_forcing:
        forcing = pickle.load(f_forcing)

    transit_list = [node for node in g.nodes if node not in comm_list]

    return forcing, comm_list, transit_list
