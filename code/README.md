![alt text](https://github.com/aleable/McOpt/blob/main/misc/logo.svg)

## Implementation details

### Table of Contents  

- [What's included](#whats-included)  
- [Requirements](#requirements)  
- [How to use](#how-to-use)  
    - [Parameters](#parameters)  
- [I/O format](#io-format)  
    - [Input](#input)
    - [Output](#output)
- [Additional implementation details](#additional-implementation-details)
    - [Convergence criteria](#convergence-criteria)
- [Usage examples](#usage-examples)  
- [Contacts](#contacts)
- [License](#license)


## What's included

- ```dashboard.ipynb```: Jupyter notebook containing an easy-to-use interface with McOpt
- ```dashboard_misc.ipynb```: complementary functions needed by ```dashboard.ipynb```
- ```main.py```: main function containing the McOpt class
- ```initialization.py```: initialization of the multicommodity problem, i.e. construction of the graph topology and on the forcings
- ```dynamics.py```: finite difference scheme (Dynamics)
- ```optimization.py```: fixed-point scheme (Optimization)

## Requirements

All the dependencies needed to run the algorithms can be installed using ```setup.py```.
This script can be executed with the command:

```bash
python setup.py
```

Now, you are ready to use the code! To do so, you can simply use the notebook ```dashboard.ipynb```, from which you can access our solvers. <br/>

## How to use

### Parameters

The parameters you can pass to McOpt are described in detail in ```dashboard.ipynb```. They are:

- *Problem construction*
    - ```method``` = ```"synth"``` or ```"paris"```: choose which graph topology you want to build, use ```"synth"``` to generate a [Waxman graph](https://networkx.org/documentation/stable/reference/generated/networkx.generators.geometric.waxman_graph.html), ```"paris"``` to build the Paris metro network
    - ```coupling``` = ```"l1"``` or ```"l2"```: choose which function you want to use to couple the fluxes, ```"l1"``` the 1-norm squared, ```"l2"``` to the 2-norm squared
    - ```museed``` (```type=int```): seed for the conductivities random initialization
    - ```pflux``` (```type=float```): 0 < β < 2, regulatory parameter for traffic congestion
    - ```rho``` (```type=float```): 0 ≤ ρ ≤ 1, regulatory parameter for the input loads

- *Dynamics parameters*
    - ```time_step``` (```type=float```): time step for the finite difference discretization of the dynamics
    - ```tot_time``` (```type=int```): upper bound on the number of time steps (safety variable)
    - ```relax_linsys``` (```type=float```): relaxation for the weighted network Laplacian
    - ```tau_cond_dyn``` (```type=float```): threshold for convergence of Dynamics, using the conductivities
    - ```tau_cost_dyn``` (```type=float```):  threshold for convergence of Dynamics, using the cost

- *Fixed-point (Optimization) scheme parameters*
    - ```tau_cond_opt``` (```type=float```): threshold for convergence of Optimization, using the conductivities
    - ```tau_cost_opt``` (```type=float```): threshold for convergence of Optimization, using the cost

- *Misc*
    - ```VERBOSE``` (```type=bool```): print info while running the schemes

## I/O format

### Input

If you want to test the code on the Paris Metro network you need these input files:

- ```adj.dat```: contains the adjacency list of the graph. Each row is of the form ```u v```, with ```u``` and ```v``` nodes connected by an edge ```e```
- ```coord.dat```: contains the nodes coordinates. Each row has to be of the form ```x[v] y[v]```, with ```x[v]``` and ```y[v]``` cartesian coordinates of ```v```
- ```rhs.dat```: contains information on the mass entering each commodity. Each row has to be formatted as: ```v g[v] 0```, where the first argument is the commodity (node) index and the second is the mass entering in it. The third dummy column, conventionally fixed to 0, may contain mass exiting from each commodity. These quantities have been always synthetically generated in our implementation, however this is only a conventional choice

### Output

The outputs returned by our schemes are:

- ```flux_dyn```: optimal fluxes returned by Dynamics
- ```flux_opt```: optimal fluxes returned by Optimization, out fixed-point iteration scheme
- ```graph```: graph topology
- ```length```: length of the edges
- ```forcing```: input forcing matrix
- ```cost_dyn```: evolution of the cost along a solution trajectory of Dynamics
- ```cost_opt```: evolution of the cost along a solution trajectory of Optimization

These can be serialized using the functions ```export_flux()``` and ```export_cost()``` defined in ```main.py```.


## Additional implementation details

### Convergence criteria

The convergence criteria we choose for our algorithms is the following. We stop Dynamics when the cost difference between two consecutive time steps is below a certain threshold *and* the maximum absolute difference between the conductivities is under another threshold (resp. the norm of the fluxes for Optimization). Precisely:
```python
# abs_diff_cost is the difference, in absolute value, of the cost in two consecutive time steps
# var_tdens is the difference of two consecutive conductivities/fluxes
if abs_diff_cost < self.tau_cost_dyn and var_tdens < self.tau_cond_dyn:
    convergence_achieved = True
```
As β approaches 2 we relax this criteria. Due to the strong concavity of the cost it is in fact progressively harder to find a unique minimum, and while the energy gets minimized we may jump between different feasible solutions of the fluxes. Thus, we impose:
```python
# we do not set any condition on the conductivities/fluxesm the lower bound on the time steps is set for safety
if self.pflux >= 1.0:
    if abs_diff_cost < self.tau_cost_dyn and time_iteration > 30:
        convergence_achieved = True
```

## Usage examples

For a basic usage example of the code you can simply take a look look at ```dashboard.ipynb```. <br/>
The execution of McOpt is performed in three steps:
- *Initialization*: first, you need to pass the necessary parameters to the McOpt class. Similarly to what done in ```dashboard.ipynb```, you can simply run:

```python
from main import *  # import McOpt class
# [...] you may want to load other packages

# Problem construction
method = "paris"    # here we choose to construct the Paris metro topology
coupling = "l2"     # we couple the commoditis with the 2-norm of their fluxes
museed = 0
pflux = 1.0
rho = 0.0

# Dynamics parameters
time_step = 0.5
tot_time = 10000
relax_linsys = 1.0e-5
tau_cond_dyn = 1.0e-1
tau_cost_dyn = 1.0e-1

# Fixed-point (Optimization) scheme parameters
tau_cond_opt = 1.0e-1
tau_cost_opt = 1.0e-1

VERBOSE = False

# create an object
mcopt = McOpt(method, coupling, museed, pflux, VERBOSE, rho,
              relax_linsys, tau_cond_dyn, tau_cost_dyn, time_step, tot_time,
              tau_cond_opt, tau_cost_opt)

# construct the topology
mcopt.ot_setup()
```
- *Execution*: now, you can execute our schemes. You just need to run:
```python
mcopt.dyn_exec()    # execute the finite difference discretization of the dynamics
mcopt.opt_exec()    # run the fixed-point iterations
```

- *Serialization*: once the algorithms reached convergence, you can export the results with:

```python
# export the fluxes at convergence with the two schemes, the graph topology,
# the length of the edges, and the forcing matrix
flux_dyn, flux_opt, graph, length, forcing = mcopt.export_flux()

# export the cost evolution
cost_dyn, cost_opt = mcopt.export_cost()
```

## Contacts

For any issues or questions, feel free to contact us sending an email to <a href="alessandro.lonardi@tuebingen.mpg.de">alessandro.lonardi@tuebingen.mpg.de</a>.

## License

Copyright (c) 2021 <a href="https://aleable.github.io/">Alessandro Lonardi</a>, <a href="https://www.cdebacco.com/">Caterina De Bacco</a> and <a href="https://enricofacca.github.io/">Enrico Facca</a>

<sub>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:</sub>

<sub>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.</sub>

<sub>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</sub>
