# `back-pvt`

`back-pvt` is a C++ library and backend tools for reservoir/production
fluids characterization through the Black Oil Method.

## Main Features

- Runtime entirely written in C++17.

- Fully reentrant code, which makes it multithreaded.

- Characterization of the following reservoir fluids: black oils, dry
gases and wet gases. 

- Twenty-six PVT properties of crude oil, natural gas and formation
brines.

- Comprehensive set of correlations, methods and standard equations. 

- Backend in command line interface for generating high-resolution,
high-precision petroleum fluid datasets of twenty-six PVT properties
of crude oil, natural gas and formation brines through
high-performance processing. 

- Backend command line interface for generating high-resolution,
high-precision petroleum fluid datasets by applying model calibration
to Constant Composition Expansion (CCE) and Differential Liberation
(DL) fluid experimental data at one or various temperatures. 

- Tested on gcc, clang and Intel (icc), compilable on any platform
that supports any of these.

- Comprehensive set of correlations, methods and standard equations.

- Compilable in Windows, OSX, Linux, and Solaris operating systems.

## Building

### Dependencies

`back-pvt` requires having installed the following dependencies:

- GSL (GNU Scientific Library).

- GNU MPFR Library.

- Imake (`xutils-dev` on Ubuntu distributions). On others distros you
could follow [these instrutions](http://www.snake.net/software/imake-stuff/imake-faq.html#where-to-get).

- `Aleph-w` <https://github.com/lrleon/Aleph-w>.

- `uconv` units conversion library <https://github.com/lrleon/uconv>.

- Niels Lohmann (nlohmann) json library <https://github.com/nlohmann/json>.

- Ruby and the following gems: `bibtex-ruby`, `calculus`, `citeproc`,
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
possible aggressive optimizations. 

## List of backends

### `test-corr`

It contains all the defined correlations and allows to evaluate them.

Examples:

In order to see the full list of correlations:

     ./test-corr -l
     
Given a correlation, for knowing the expected parameters and its technical note:

    ./test-corr -C BobMillanArcia -P

    MILLÁN-ARCIA CORRELATION, CALCULATION OF SATURATED OIL
    FORMATION VOLUME FACTOR

    Millán-Arcia
    OilCorrelation : SaturatedOilVolumeFactor : BobMillanArcia
      Return type = RB/STB
      Min result  = 1.014
      Max result  = 1.248

      Parameters (4):
         name    unit       min           max number 
                                                 
          api     api     9 api      20.2 api      1 
          rsb scf/STB 0 scf/STB 1e+09 scf/STB      2 
            p    psia    0 psia  29985.3 psia      3 
           pb    psia  222 psia   3432.7 psia      4 

      DATA BANK:

        Venezuelan heavy crudes were correlated.

      References:

        Millán-Arcia, E. A. (1984). Correlaciones para estimar el comportamiento PVT de crudos pesados venezolanos. Corpoven.

        Pérez, V. L., Heny, C. A., & Lago, M. E. (2001). Evaluación y generación de correlaciones para propiedades PVT y viscosidades de crudos extrapesados de la Faja Petrolífera del Orinoco. Mérida, Venezuela: Facultad de Ingeniería Química, Universidad de Los Andes.
	
In order to evaluate the correlation given some parameters:

    ./test-corr -C BobMillanArcia 10 200 600 800
    
    BobMillanArcia(10 api, 200 scf/STB, 600 psia, 800 psia) = 1.24475 RB/STB
    
If you want to use another unit, then you can use the `-u` option. For
example, if you want to specify the pressure in psig, then perform

    ./test-corr -C BobMillanArcia 10 200 600 800 -u "3 psig"
    
    BobMillanArcia(10 api, 200 scf/STB, 600 psig, 800 psia) = 1.24616 RB/STB
    
### Fluid tuner

### Z factor tuner

### Grid generator
