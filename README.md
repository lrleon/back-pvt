# `back-pvt`

`back-pvt` is a C++ library and backend tools for reservoir/production
fluids characterization through the Black Oil Method.

## Main Features

. Runtime entirely written in C++14.

. Fully reentrant code, which makes it multithreaded.

. Characterization of the following reservoir fluids: black oils, dry
gases and wet gases. 

. Twenty-six PVT properties of crude oil, natural gas and formation
brines.

. Comprehensive set of correlations, methods and standard equations. 

. Backend in command line interface for generating high-resolution,
high-precision petroleum fluid datasets of twenty-six PVT properties
of crude oil, natural gas and formation brines through
high-performance processing. 

. Backend command line interface for generating high-resolution,
high-precision petroleum fluid datasets by applying model calibration
to Constant Composition Expansion (CCE) and Differential Liberation
(DL) fluid experimental data at one or various temperatures. 

## Building

### Dependencies

`back-pvt` requires having installed the following dependencies:

. GSL (GNU Scientific Library).

. GNU MPFR Library.

. Imake (`xutils-dev` on ubuntu distributions). On others distro you
could follow [these instrutions](http://www.snake.net/software/imake-stuff/imake-faq.html#where-to-get).

. `Aleph-w` <https://github.com/lrleon/Aleph-w>.

. `uconv` units conversion library <https://github.com/lrleon/uconv>.

. Ruby and the following gems: `bibtex-ruby`, `calculus`, `citeproc`,
`csl`, `csl-styles` and `citeproc-ruby`.

### Building

On the directory where you have the sources, perform:

	xmkmf
	make Makefiles
	make depend
	
The command `make depend` generates compilation "dependencies", which
allow you to detect when the sources have to be recompiled. Each time
you add a new header file to the `include` directory, execute `make
depend`. 

To generate the backend executables:

	make all
	
After building, the backends will be in the directories `tests` and
`backend`. The first directory contains the backends without
optimization. In the second one the backend are compiled with the most
possible agressive optimizations. 

## List of backends

### `eval-corr`

It contains all the defined correlations and allows to evaluate them.

Examples:

	
