Noda is a Python package to simulate diffusion in multicomponent alloys.
It provides tools to explore thermodynamic and mobility properties of alloys
and to solve the diffusion problem using the finite difference method in 1D.

Our aim is to provide researchers and any interested people a means to study
diffusion in alloys, and to study and modify the underlying models. For help
using Noda or other inquiries, please contact the maintainer by email.

Noda is distributed under the GPLv3 license (see LICENSE.txt).

Dependencies:

* Python 3.11+
* numpy, scipy, pandas, matplotlib, odfpy, openpyxl, numdifftools

The documentation is available at
[onera.github.io/Noda](https://onera.github.io/Noda); it contains
installation instructions, tutorials and the package reference.

Quick install
-------------

```
git clone https://github.com/onera/Noda.git
cd Noda
pip install .
```

Citing
------

If you use Noda, please cite:

T. Gheno, V. Szczepan, C. Salsi, C. Desgranges, D. Monceau,
*Simulation of diffusion with non-equilibrium vacancies, Kirkendall shift and
porosity in single-phase alloys*, Computational Materials Science 215 (2022)
111785.  
Publisher version: [10.1016/j.commatsci.2022.111785](https://doi.org/10.1016/j.commatsci.2022.111785)   
Accepted version: [hal-03792936](https://hal.science/hal-03792936)