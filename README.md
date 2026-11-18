# HELM-RT: Regulating transformers in the holomorphic embedding load flow

[![DOI](https://zenodo.org/badge/1463191471.svg)](https://doi.org/10.5281/zenodo.24354088)

MATLAB code and data for the paper

> G. R. S. Lima, E. J. Oliveira, and L. W. Oliveira, "Regulating Transformers as State Variables of the Holomorphic Embedding Load Flow: Joint Voltage and Active Power Flow Control with Operating Limits and Discrete Taps," submitted to *IEEE Open Access Journal of Power and Energy*, 2026.

The code solves the power flow with regulating transformers by the Holomorphic Embedding Load Flow Method (HELM) and by the Newton–Raphson (NR) method, with the same model and the same outer loop:

- each controlled transformer keeps its nominal complex ratio in the admittance matrix, and the deviations of the tap magnitude and of the phase angle are power series of the embedding;
- one transformer can regulate the voltage of a bus (mode V), the active-power flow through itself (mode F), or both at the same time (mode VF); parallel units at the same bus share one tap (master–follower);
- the flow equation uses the same branch current that the embedded network assigns to the transformer, so the controlled flow and the nodal balance are coherent at every order;
- the recursion matrix is sparse, constant for all orders, and factorized once per pass;
- one outer loop handles the reactive limits of the generators and the limits of tap and phase, with release of the fixed controls, and then moves the taps and phases to discrete positions with a voltage deadband;
- all power series (voltages, taps and phases) are evaluated in a single step per order by an incremental form of Wynn's epsilon algorithm.

## Requirements

- MATLAB R2018a or later. The code also runs in GNU Octave 8.
- No toolbox and no MATPOWER installation are required.

## Repository structure

| Folder / file | Content |
| --- | --- |
| `Main_*.m` | Main scripts that reproduce the results of the paper |
| `funcoes/` | Functions (one per file) |
| `dados/` | Test systems |
| `referencia/` | Data files of the paper, used by `Main_Todos.m` to check the results |
| `resultados/` | Output data files (created on the first run) |

Code comments are written in Portuguese.

## Reproducing the results

Run the main scripts from the root folder of the repository; each script adds `funcoes` and `dados` to the path and writes its data files to `resultados`.

| Script | Study |
| --- | --- |
| `Main_Modos.m` | IEEE 30, 57 and 118 buses, 24 h, two wind scenarios, modes V, F and VF, without limits |
| `Main_Limites.m` | The same 24 h studies with reactive, tap and phase limits; continuous and discrete taps |
| `Main_ACTIVSg.m` | ACTIVSg2000 with 208 load tap changers (298 transformers): loss of each of the ten most loaded 500-kV corridors at heavy (hour 5368) and light (hour 8504) load, also solved with the taps fixed to measure the action of the LTCs |
| `Main_Convergencia.m` | Convergence of Wynn's epsilon algorithm versus the plain partial sums (IEEE systems and ACTIVSg2000) |
| `Main_Benchmark.m` | Computing times of one solution (median of repeated runs) |
| `Main_Todos.m` | Runs all scripts above except `Main_Benchmark.m` and checks every output file against `referencia/` |

`Main_Todos.m` takes a few minutes and reports, for each data file, the largest difference from the values of the paper. In every study the script also checks that HELM and NR reach the same set of saturated generators, taps and phases, with the same number of passes of the outer loop.

Computing times depend on the MATLAB version and on the computer; they are printed on the screen and are not part of the reference files.

## Main functions

```matlab
[V, a, f, info] = HELMT_Resolver(D, ctrl, lim, metodo, tol, nmax, discreto)
```

- `D`: system data built by `Prepara_dados` from a case in MATPOWER format;
- `ctrl`: controlled transformers (fields `ramo`, `modo` = `'V'`, `'P'` or `'VP'`, `barra`, `Vesp`, `Psp`);
- `lim`: limits of tap and phase, steps, voltage deadband and maximum number of step moves;
- `metodo`: `'HELM'` or `'NR'`; `discreto`: `true` for discrete taps and phases;
- `V`: complex bus voltages (p.u.); `a`, `f`: tap magnitude and phase (rad) of each control;
- `info`: convergence flag, passes of the outer loop, number of HELM orders or NR iterations, saturated generators, taps and phases, and the event log.

`HELMT_Nucleo` and `NRT_Nucleo` solve one pass with fixed bus types and are called by `HELMT_Resolver`.

## Test systems

| File | Description | Source |
| --- | --- | --- |
| `case30.m`, `case57.m`, `case118.m` | IEEE 30, 57 and 118-bus systems | MATPOWER |
| `case_ACTIVSg2000_h5368.m`, `case_ACTIVSg2000_h8504.m` | Synthetic 2000-bus Texas system, hours 5368 (annual peak) and 8504 (light load) of the 2016 time series | [HELM-AIC](https://github.com/HELM-UFJF/HELM-AIC) ([doi:10.5281/zenodo.24318005](https://doi.org/10.5281/zenodo.24318005)) |

The controlled transformers, specified voltages, wind farms and the 24-h load curve of the IEEE systems are set in `funcoes/Config_Sistema.m`. The load tap changers of the ACTIVSg2000 are built by `funcoes/Caso_ACTIVSg.m`.

## License

The code is released under the MIT License (see `LICENSE`). The test system data in `dados/` keep the licenses of their original sources: the IEEE case files are distributed with MATPOWER, and the ACTIVSg2000 data are licensed under CC BY 4.0.

## Citation

If you use this code, please cite the paper above and the archived code ([doi:10.5281/zenodo.24354088](https://doi.org/10.5281/zenodo.24354088)). The reference will be updated after publication. The metadata for citing this repository are in `CITATION.cff`.

## Contact

Glaucus R. S. Lima, Graduate Program in Electrical Engineering, Federal University of Juiz de Fora (UFJF), Brazil.
