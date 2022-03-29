# poLYG
Polygon area minimization and maximization as in [CG:SHOP challenge 2019](https://cgshop.ibr.cs.tu-bs.de/competition/cg-shop-2019/#problem-description). This algorithm got the 2nd place in the challenge. The paper descring the results is available at http://dx.doi.org/10.1145/3503999 and also at https://arxiv.org/abs/2106.14728

Team:
- Loic Crombez
- [Yan Gerard](https://yangerard.wordpress.com/)
- [Guilherme Dias da Fonseca](https://pageperso.lis-lab.fr/guilherme.fonseca/)

```
pypy3 poLYG.py [args] inputfile

Inputfile extensions are ignored. Args are of the form variable=value (no space, no dash).
In fact, they will be executed with python exec(arg), so feel free to add formulas that depend on n or other parameters.

Some variables and their default values are:

maximize=True   % Maximum or minimum area
pen=90          % Parameter 1/alpha of the weight function
hood=2          % Neighborhood kappa of an edge
opt=True        % Apply local search optimization
hops=1          % Value of ell for the local search
multirun=False  % Run many times
sigma=0         % Gaussian noise of the weight function
seed=1          % Random seed
timeout=150     % Maximum number of seconds to start a new run
```
