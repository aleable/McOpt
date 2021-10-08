"""
McOpt (Multicommodity Optimal Transport) -- https://github.com/aleable/McOpt

Contributors:
    Alessandro Lonardi
    Enrico Facca
    Caterina De Bacco
"""

import numpy as np
import networkx as nx
import os
import pickle5 as pkl
from scipy.spatial import distance


def waxman_topology(self):
    """Generation of the Waxman graph topology

    Parameters:
        nnode: int, number of nodes

    Returns:
        self.g: nx.Graph(), Waxman graph topology
        self.length: np.array, lengths of edges
    """

    # graph topology construction
    self.g = nx.waxman_graph(n=30, alpha=1, beta=0.25, L=1.0, domain=(0, 0, 1, 1), seed=0)

    self.length = np.zeros(self.g.number_of_edges())
    for i, edge in enumerate(self.g.edges()):
        self.length[i] = distance.euclidean(self.g.nodes[edge[0]]["pos"], self.g.nodes[edge[1]]["pos"])

    # right hand side construction
    dummy_path = " "
    self.forcing = forcing_generation(self, dummy_path)

    return self.g, self.length, self.forcing


def paris_topology(self, input_path):
    """Generation of the Paris metro network topology

    Parameters:
        input_path: string, input folder path

    Returns:
        self.g: nx.Graph(), Waxman graph topology
        self.length: np.array, lengths of edges
    """

    adj_file = open(input_path + "adj.dat", "r")
    lines = adj_file.readlines()

    # graph adjacency list
    topol = np.zeros([len(lines), 2], dtype=int)
    for iedge, line in enumerate(lines):
        topol[iedge][:] = [int(w) for w in line.split()[0:2]]

    self.g.add_edges_from(topol)

    # coordinate of nodes
    coord_file = open(input_path + "coord.dat", "r")
    lines = coord_file.readlines()
    for inode, line in enumerate(lines):
        self.g.nodes[inode]["pos"] = tuple([float(w) for w in line.split()[0:2]])

    # length of edges
    self.length = np.zeros(self.g.number_of_edges())
    for i, edge in enumerate(self.g.edges()):
        self.length[i] = distance.euclidean(self.g.nodes[edge[0]]["pos"], self.g.nodes[edge[1]]["pos"])

    # right hand side construction
    forcing_path = input_path + "rhs.dat"
    self.forcing = forcing_generation(self, forcing_path)

    return self.g, self.length, self.forcing


def forcing_generation(self, forcing_path):
    """Construct forcings using the "Influence Assignment" method

    Parameters:
        forcing_path, str, path of the rhs file (dummy string if Waxman model)

    Returns:
        self.forcing: np.array(nnode, ncomm), forcing matrix
    """

    # paris metro influence assignment
    if self.method == "paris":

        with open(forcing_path) as forcing_f:
            lines = forcing_f.readlines()
            ncomm = len(lines)
            temp_forcing = np.zeros((ncomm, 3))
            for i, line in enumerate(lines):
                temp_forcing[i, :] = np.array([int(float(element)) for element in line.strip().split(" ")])

        # sources and sinks
        comm_list = temp_forcing[:, 0].astype(int)
        pkl.dump(comm_list, open(os.getcwd() + "/../data/input/comm_list.pkl", "wb"))
        # inflowing mass
        g_forcing = temp_forcing[:, 1]

        self.forcing = rhs_construction(self, ncomm, comm_list, g_forcing)

        return self.forcing

    # Waxman graphs influence assignment influence assignment
    if self.method == "synth":

        ncomm = 30
        comm_list = list(self.g.nodes())
        pkl.dump(comm_list, open(os.getcwd() + "/../data/input/comm_list.pkl", "wb"))
        prng = np.random.RandomState(seed=0)
        g_forcing = np.array([prng.uniform(0, 1)*100 for i in range(self.g.number_of_nodes())])

        self.forcing = rhs_construction(self, ncomm, comm_list, g_forcing)

        return self.forcing


def rhs_construction(self, ncomm, comm_list, g_forcing):
    """Forcing matrix construction

    Parameters:
        ncomm: int, number of commodities
        comm_list: np.array, commodities list
        g_forcing: np.array, entry mass vector

    Returns:
        self.forcing: np.array(nnode, ncomm), forcing matrix
    """

    g_forcing = g_forcing - self.rho*(g_forcing - np.mean(g_forcing))

    self.forcing = np.zeros((self.g.number_of_nodes(), ncomm))

    for n in range(ncomm):
        sum_g = np.sum(g_forcing)
        sum_g -= g_forcing[n]
        coeff_r = g_forcing / sum_g
        self.forcing[comm_list, n] = -coeff_r * g_forcing[n]
        self.forcing[comm_list[n], n] = g_forcing[n]

    return self.forcing