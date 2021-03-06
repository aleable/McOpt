# McOpt (Multi-commodity Optimal Transport)

McOpt /mæk ɒpt/ is a Python implementation of the algorithms described in:

- [1] Alessandro Lonardi, Enrico Facca, Mario Putti and Caterina De Bacco. <i>Optimal transport for multi-commodity routing on networks</i>.

This is a new algorithm capable of finding optimal multi-commodity flows on networks solving a dynamical system of equations, or solving an optimization problem by means of a fast iterative update of the flows. Our implementations proved to be more efficient than state-of-the-art methods employed to solve optimal transport problems.

If you use this code please cite [1].

The paper can be found <a href="https://arxiv.org/abs/2010.14377">here</a> (preprint).

See ```code/data_visualization.ipynb``` for an example result on the Paris metro.

## What's included

- ```code```: contains the all the Python scripts necessary to run McOpt and a Jupyter notebook ```data_visualization.ipynb``` with which it is possible to visualize the converged networks.
- ```data/input```: contains a network (Paris metro) on which the algorithm can be tested. The network topology has been pre-processed and extracted using [2]. Passengers usage data have been taken from [3].
- ```data/output```: repository used to store some of the metrics we are able to study with McOpt.

[2] R. Kujala, C. Weckström, R. K. Darst, M. N. Mladenović, and J. Saramäki, Scientific data <b>5</b>, 180089 (2018).<br/>
[3]  <a href="https://data.ratp.fr/explore/dataset/trafic-annuel-entrant-par-station-du-reseau-ferre-2019/information/"> “Trafic annuel entrant par station du réseau ferré 2019”</a>, (2019), accessed: 2020-08-28.

## Requirements

The project has been developed using Python 3.7 with the packages contained in ```requirements.txt```. It is possible to install the dependencies using pip:
```bash
pip install -r requirements.txt
```

## How to use

To download this repository, copy and paste the following:

```bash
git clone https://github.com/aleable/McOpt
```

You are ready to test the program on the given dataset, to do so type in the algorithm repository:

```bash
cd code
python main.py
```

This will result in running the scripts with all parameters initialized as default. For more information see in ```code``` repo.

### Flags

A list of all the flags can be produced with the command
```bash
python main.py --help
```

- ```-fg```: imported graph file name (```default="graph.dat"```)
- ```-fa```: imported adjacency list file name (```default="adj.dat"```)
- ```-fs```: imported right hand side file name (```default="rhs.dat"```)
- ```-fc```: imported coordinate file name (```default="coord.dat"```)
- ```-topol```: ```0``` generate topology (Waxman graph), ```1``` import topology from file (```default=1```); <br/>
  if ```-topol 0``` change the following flags as follows:
  - ```-fg "graph_generated.dat"```
  - ```-fa "adj_generated.dat"```
  - ```-fs "rhs_generated.dat"```
  - ```-fc "coord_generated.dat"```
- ```-N```: number of nodes in Waxman graph, ignored if ```-topol 1``` (```default=30```)
- ```-S```: number of commodities in Waxman graph, ignored if ```-topol 1``` (```default=10```)
- ```-mass```: ```0```/```1``` generate/import entering mass file, ignored if ```-a im``` (```default=1```)
- ```-coord```: ```0```/```1``` generate/import node coordinates, forced to ```1``` if ```-topol 0``` (```default=1```)
- ```-a```: method used to build matrix S (```default=ia```):
  - ```fr``` Fake Receiver method
  - ```ia``` Influence Assignment
  - ```im``` import forcing from previous run
- ```-M```: total entering mass randomly assigned to the commodities, ignored if ```-mass 1``` or ```-a im``` (```default=1000```)
- ```-l```: ```bias```/```eucl``` edge lengths computed as  1 + 0.001 * U(0,1)/using euclidean distance between nodes (```default=eucl```)
- ```-b```: beta (```default=0.5```)
- ```-t```: trimming thresholds (```default=0.001```)
- ```-d```: exporting quantities for data analysis (```default=0```):
  - ```0``` no export
  - ```1``` useful variable for graph visualization
  - ```2``` cost of the converged graphs
  - ```3``` number of basis loops
  - ```4``` number of idle edges
- ```-v```: ```0```/```1``` do not print/print running script comments (```default=0```)

## I/O format

### Input

Four input files are necessary to run the scripts. These should be formatted as follows.
- ```graph.dat```: contains the number of nodes, commodities and edges of the input graph. This has to be formatted with these tree quantities in subsequent rows:<br/>
```nnode```<br/>
```nncomm```<br/>
```nedge```
- ```adj.dat```: contains the adjacency list of the graph. Each row has to be of the form ```e[0] e[1]```, with ```e[0]```/```e[1]``` nodes connected by an edge ```e```.
- ```coord.dat```: contains the nodes coordinates. Each row has to be of the form ```x[inode] y[inode]```, with ```x[inode]```/```y[inode]``` x and y cartesian coordinates of node ```inode```.
- ```rhs.dat```: contains information of mass entering each commodity. Each row has to be formatted as: ```inode g[inode] 0```, where the first argument is the commodity index and the second is the mass entering ```inode```. The third dummy column, conventionally fixed to 0, may contain mass exiting from each commodity. These have been always synthetically generated in our implementation (using ```fr```/```ia```), however future developments of the code where this choice can be modified are left to practitioners.
- ```paris_nodes.csv```: file containing a list of the names of each Paris metro station. Each row is formatted as ```inode name[inode]```, with ```inode``` index of the station as used in ```adj.dat``` and ```name[inode]``` its name.
### Output

Changing value for the flag ```-d``` we can generate the following files:
- ```1```:
  - ```flux_norm_dyn.pkl```, ```flux_norm_opt.pkl``` pickled files containing the fluxes of the trimmed networks.
  - ```length_edges_dyn.pkl```, ```length_edges_opt.pkl``` pickled files containing the edge lengths of the trimmed networks.
  - ```comm_list.pkl``` pickled files containing a list with the commodities indexes.
  - ```graph_dynamics_dyn.gpickle```, ```graph_dynamics_opt.gpickle``` pickled files containing the topologies of the trimmed networks.
- ```2```: ```cost_dyn.dat```, ```cost_opt.dat``` files containing the cost of the trimmed networks. Formatted as: ```beta cost```.
- ```3```: ```cycles_dyn.dat```, ```cycles_opt.dat``` files containing the number of basis cycles of the trimmed networks. Formatted as: ```beta #cycles```.
- ```4```: ```idle_dyn.dat```, ```idle_opt.dat``` files containing the number of idle edge (links that gets trimmed). Formatted as:  ```beta idle edges```.

## Contacts

For any issues or questions, feel free to contact us sending an email to <a href="alessandro.lonardi@tuebingen.mpg.de
">alessandro.lonardi@tuebingen.mpg.de</a>.

## License

Copyright (c) 2020 <a href="https://aleable.github.io/">Alessandro Lonardi</a>, <a href="https://www.cdebacco.com/">Caterina De Bacco</a> and Enrico Facca.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
