![alt text](https://github.com/aleable/McOpt/blob/main/misc/logo.svg)

___

McOpt /mæk ɒpt/ (Multicommodity Optimal Transport) is a Python implementation of the algorithms used in:

- [1] Alessandro Lonardi, Enrico Facca, Mario Putti and Caterina De Bacco. <i>Designing optimal networks for multicommodity transport problem</i>. <a href="https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.043010">Phys. Rev. Research <b>3</b>, 043010</a> [<a href="https://arxiv.org/abs/2010.14377">arXiv</a>].
- [2] Alessandro Lonardi, Mario Putti and Caterina De Bacco. <i>Multicommodity routing optimization for engineering networks</i>. [<a href="https://arxiv.org/abs/2110.06171">arXiv</a>].

This is a scheme capable of finding optimal multicommodity flows on networks solving a dynamical system of equations, or with a fixed-point scheme for the update of the flows. Our implementations is competitive with many state-of-the-art methods employed to solve optimal transport problems.

**If you use this code please cite [1].**

## What's included

- ```code```: contains the all the scripts necessary to run McOpt, and a Jupyter notebook (```dashboard.ipynb```) with which it is possible to easily interact with them and visualize the results
- ```data/input```: contains the data needed to build the Paris metro network, where the algorithm can be tested. The network topology has been pre-processed and extracted using [3]. Passengers usage data have been taken from [4]
- ```misc```: files used for the README.md
- ```setup.py```: setup file to build the Python environment

[3] Rainer Kujala, Christoffer Weckström, Richard K. Darst, Miloš N. Mladenović and Jari Saramäki, <a href="https://www.nature.com/articles/sdata201889">Scientific data <b>5</b>, 180089 (2018)</a>.<br/>
[4]  <a href="https://data.ratp.fr/explore/dataset/trafic-annuel-entrant-par-station-du-reseau-ferre-2019/information/"> “Trafic annuel entrant par station du réseau ferré 2019”</a>, (2019), accessed: 2020-08-28.

## How to use

To download this repository, copy and paste the following:

```bash
git clone https://github.com/aleable/McOpt
```


**You are ready to test the code! But if you want to know how click [HERE](https://github.com/aleable/McOpt/tree/main/code)**.

## Contacts

For any issues or questions, feel free to contact us sending an email to <a href="alessandro.lonardi@tuebingen.mpg.de">alessandro.lonardi@tuebingen.mpg.de</a>.

## License

Copyright (c) 2021 <a href="https://aleable.github.io/">Alessandro Lonardi</a>, <a href="https://www.cdebacco.com/">Caterina De Bacco</a> and <a href="https://enricofacca.github.io/">Enrico Facca</a>

<sub>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:</sub>

<sub>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.</sub>

<sub>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</sub>
