# Results (CI)

## Known answers

| | L/h | Euler-Bernoulli (mm) | Timoshenko, Cowper k (mm) | Timoshenko, k = 5/6 (mm) | shear share |
|---|---|---|---|---|---|
| beam | 10 | 0.952381 | 0.959667 | 0.959810 | 0.76% |
| stocky | 4 | 0.060952 | 0.063867 | 0.063924 | 4.56% |

Bracket, pull load: net-section stress 41.667 MPa, Heywood Kt 2.4219.

## Cantilevers: tip deflection (centre of the tip face)

| case | solver | mesh | dofs | tip (mm) | vs Timoshenko | vs Euler-Bernoulli | solve (s) |
|---|---|---|---|---|---|---|---|
| beam | CalculiX C3D10 | coarse | 3,030 | 0.950655 | -0.939% | -0.181% | 0.12 |
| beam | CalculiX C3D10 | medium | 19,737 | 0.952146 | -0.784% | -0.025% | 0.60 |
| beam | CalculiX C3D10 | fine | 118,878 | 0.952856 | -0.710% | +0.050% | 9.95 |
| beam | FEniCSx P2 | coarse | 3,030 | 0.950655 | -0.939% | -0.181% | 0.04 |
| beam | FEniCSx P2 | medium | 19,737 | 0.952146 | -0.784% | -0.025% | 0.43 |
| beam | FEniCSx P2 | fine | 118,878 | 0.952856 | -0.710% | +0.050% | 6.35 |
| beam | CalculiX C3D4 | coarse | 573 | 0.467809 | -51.253% | -50.880% | 0.01 |
| beam | CalculiX C3D4 | medium | 3,222 | 0.798296 | -16.815% | -16.179% | 0.08 |
| beam | CalculiX C3D4 | fine | 17,148 | 0.906703 | -5.519% | -4.796% | 0.66 |
| beam | FEniCSx P1 | coarse | 573 | 0.467809 | -51.253% | -50.880% | 0.00 |
| beam | FEniCSx P1 | medium | 3,222 | 0.798296 | -16.815% | -16.179% | 0.02 |
| beam | FEniCSx P1 | fine | 17,148 | 0.906703 | -5.519% | -4.796% | 0.28 |
| stocky | CalculiX C3D10 | coarse | 1,518 | 0.062368 | -2.347% | +2.322% | 0.03 |
| stocky | CalculiX C3D10 | medium | 8,259 | 0.062655 | -1.897% | +2.793% | 0.22 |
| stocky | CalculiX C3D10 | fine | 48,216 | 0.062777 | -1.707% | +2.993% | 3.57 |
| stocky | FEniCSx P2 | coarse | 1,518 | 0.062368 | -2.347% | +2.322% | 0.02 |
| stocky | FEniCSx P2 | medium | 8,259 | 0.062655 | -1.897% | +2.793% | 0.16 |
| stocky | FEniCSx P2 | fine | 48,216 | 0.062777 | -1.707% | +2.993% | 2.16 |
| stocky | CalculiX C3D4 | coarse | 285 | 0.037996 | -40.507% | -37.663% | 0.01 |
| stocky | CalculiX C3D4 | medium | 1,371 | 0.052557 | -17.708% | -13.773% | 0.04 |
| stocky | CalculiX C3D4 | fine | 7,056 | 0.059571 | -6.727% | -2.267% | 0.26 |
| stocky | FEniCSx P1 | coarse | 285 | 0.037996 | -40.507% | -37.663% | 0.00 |
| stocky | FEniCSx P1 | medium | 1,371 | 0.052557 | -17.708% | -13.773% | 0.01 |
| stocky | FEniCSx P1 | fine | 7,056 | 0.059571 | -6.727% | -2.267% | 0.10 |

## Bracket

| solver | mesh | dofs | hole size (mm) | shelf: tip uz (mm) | pull: tip ux (mm) | Kt | Kt vs Heywood | Kt at the faces | shelf: fillet s1 (MPa) | shelf: hole s1 (MPa) | solve (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CalculiX C3D10 | coarse | 18,951 | 2.0 | -0.664070 | 0.0244338 | 2.6297 | +8.58% | 2.2335 | 93.95 | 61.50 | 0.96 |
| FEniCSx P2 | coarse | 18,951 | 2.0 | -0.664058 | 0.0244337 | 2.6553 | +9.64% | 2.1582 | 95.26 | 62.39 | 0.82 |
| CalculiX C3D10 | medium | 80,898 | 1.0 | -0.665197 | 0.0244492 | 2.5521 | +5.38% | 2.2399 | 95.30 | 62.32 | 10.82 |
| FEniCSx P2 | medium | 80,898 | 1.0 | -0.665196 | 0.0244492 | 2.5641 | +5.87% | 2.2180 | 95.77 | 62.80 | 6.40 |
| CalculiX C3D10 | fine | 278,649 | 0.6 | -0.665479 | 0.0244522 | 2.5523 | +5.38% | 2.2277 | 94.63 | 61.87 | 95.10 |
| FEniCSx P2 | fine | 278,649 | 0.6 | -0.665479 | 0.0244522 | 2.5510 | +5.33% | 2.2359 | 94.75 | 61.86 | 47.67 |

## Experiment: the beams with Poisson's ratio 0 (CalculiX C3D10)

| case | mesh | tip (mm) | Timoshenko, nu = 0 (mm) | vs Timoshenko | vs Euler-Bernoulli |
|---|---|---|---|---|---|
| beam | coarse | 0.957710 | 0.958095 | -0.040% | +0.560% |
| beam | medium | 0.958047 | 0.958095 | -0.005% | +0.595% |
| beam | fine | 0.958054 | 0.958095 | -0.004% | +0.596% |
| stocky | coarse | 0.063098 | 0.063238 | -0.221% | +3.520% |
| stocky | medium | 0.063195 | 0.063238 | -0.069% | +3.679% |
| stocky | fine | 0.063198 | 0.063238 | -0.064% | +3.684% |

## Experiment: hole Kt against plate thickness (pull load, fine meshes)

| case | t/d | solver | Kt, peak | Kt, mean of the 8 through-thickness bins | Kt at the faces | Heywood |
|---|---|---|---|---|---|---|
| thin | 0.2 | CalculiX C3D10 | 2.4609 (+1.6%) | 2.4462 (+1.0%) | 2.4047 | 2.4219 |
| thin | 0.2 | FEniCSx P2 | 2.4608 (+1.6%) | 2.4472 (+1.0%) | 2.4086 | 2.4219 |
| bracket | 0.8 | CalculiX C3D10 | 2.5523 (+5.4%) | 2.4939 (+3.0%) | 2.2277 | 2.4219 |
| bracket | 0.8 | FEniCSx P2 | 2.5510 (+5.3%) | 2.4975 (+3.1%) | 2.2359 | 2.4219 |

## Same mesh, two solvers: tip displacement (largest component)

| case | mesh | load | CalculiX (mm) | FEniCSx (mm) | relative difference |
|---|---|---|---|---|---|
| beam | coarse | tip | -0.950654600 | -0.950654628 | 3.0e-08 |
| beam | medium | tip | -0.952145700 | -0.952145720 | 2.1e-08 |
| beam | fine | tip | -0.952855900 | -0.952855925 | 2.7e-08 |
| stocky | coarse | tip | -0.062367640 | -0.062367641 | 2.2e-08 |
| stocky | medium | tip | -0.062655080 | -0.062655085 | 7.5e-08 |
| stocky | fine | tip | -0.062776720 | -0.062776718 | 3.8e-08 |
| bracket | coarse | shelf | -0.664070100 | -0.664057876 | 1.8e-05 |
| bracket | coarse | pull | 0.024433780 | 0.024433712 | 2.8e-06 |
| bracket | medium | shelf | -0.665196700 | -0.665196035 | 1.0e-06 |
| bracket | medium | pull | 0.024449200 | 0.024449197 | 1.1e-07 |
| bracket | fine | shelf | -0.665479300 | -0.665479209 | 1.4e-07 |
| bracket | fine | pull | 0.024452170 | 0.024452167 | 1.4e-07 |
| thin | fine | shelf | -43.900690000 | -43.900258734 | 9.8e-06 |
| thin | fine | pull | 0.196970400 | 0.196956332 | 7.1e-05 |

Plots: cantilever_convergence.png, bracket_hole_kt.png
