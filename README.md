# scaling-codes

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

For scan codes install also *python3*:

`apt install python3`

## Build codes

Enter model directory, e.g. ER model: `phen-models/er`.

Modify **NCORES** variable in the *Makefile* according to your processor core numbers.
If you have *boost* installed into non-standard place, correct **CPPFLAGS** and/or **LIBS** variables as well. 
In other case, just run `make` to build and `make test` to verify if executable is OK.

Next, go to `phen-models/scanner` and run scan with `make erscan`. 
If you want to modify parametric set, change *paramset-er.py*.

## Good luck!
