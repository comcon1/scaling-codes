# Simple phenomenological models

Extremly fast and accurate 1D C++ solver of Reaction-Diffusion equation is used here for parametric scan 
of three phenomenological models of embryonic scaling: Expansion-Repression, Contraction-Induction (original) [1], and Contraction-Induction (improved) [2].

1. D. Ben-Zvi, N. Barkai, Scaling of morphogen gradients by an expansion-repression integral feedback control, Proc. Natl. Acad. Sci. 107 (2010) 6924–6929. doi:[10.1073/pnas.0912734107](https://doi.org/10.1073/pnas.0912734107).
2. Orlov, et al., **in press**

## Install prerequisites

You need to install *boost* headers, *boost_program_options* libraries, and *LAPACK* C libraries.

First, clone the repository:

`git clone https://github.com/comcon1/scaling-codes`

Then init all the submodules:

`git submodule update --init`

Install required pacakges. For ubuntu/debian use the following packages:

`apt install g++ libboost-dev liblapack3-dev`

The last package is named `liblapack-dev` in some distributions. For scan codes install also *python3*:

`apt install python3`

## Build codes and run the scanner

Enter model directory, e.g. ER model: `phen-models/er`.

If you have *boost* installed into non-standard place, correct **CPPFLAGS** and/or **LIBS** variables as well. 
In other case, just run `make` to build and `make test` to verify if executable is OK.

Before starting the scan, go to `phen-models/scanner` and modify **NCORES** variable in the *Makefile* according to your processor core numbers. Next run scan with `make erscan`. 
If you want to modify parametric set, change *paramset-er.py*.

# Complex model of mesoderm patterning

Lorem ipsum dolor sit amet..

## Install prerequisites

Install the same packages as required for simple models (see above). C++ code of complex model does not require any additional libraries. At the same time, code of the scanner requires `scipy` and `numpy`.
The scanner uses **dual-annealling** stochastic optimizer that requires *scipy* v. 1.2.0+. If your *scipy* version is lower than upgrade it, e.g., via pip:

```
apt install python3-pip
pip3 install pip --upgrade
pip3 install numpy scipy
```

## Build codes and run the scanner

Make C++ solver in the same way as described in phenomonological model section.

Next, modify **NCORES** variable in the *Makefile* accordingly, and perform the scan in `scanner` directory with the command:

`make scan`

We strongly recommend to perform the scan in the detached mode. For example, when using *screen*:

`screen make scan`

To terminate the process accurately, please terminate *optimizer.py* process directly:

`pkill optimizer`

