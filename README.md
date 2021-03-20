# PPDyn (Plasma Particle Dynamics)
[![CI](https://github.com/sayanadhikari/PPDyn/actions/workflows/main.yml/badge.svg)](https://github.com/sayanadhikari/PPDyn/actions/workflows/main.yml)
[![build](https://github.com/sayanadhikari/PPDyn/actions/workflows/make.yml/badge.svg)](https://github.com/sayanadhikari/PPDyn/actions/workflows/make.yml)

A python code to simulate plasma particles with Molecular Dynamics Algorithm. [Numba JIT compiler](https://numba.pydata.org/) for Python has been implemented for faster performance. 

## Problem
<!--Rayleigh Problem = gas between 2 plates ([Alexander & Garcia, 1997](https://doi.org/10.1063/1.168619)) >>

## Contributors
- [Sayan Adhikari](https://github.com/sayanadhikari), UiO, Norway. [@sayanadhikari](https://twitter.com/sayanadhikari)
- [Rupak Mukherjee](https://github.com/RupakMukherjee), PPPL, USA.

Installation
------------
#### Prerequisites
1. [GNU Make](https://www.gnu.org/software/make/)
2. [python3 or higher](https://www.python.org/download/releases/3.0/)
3. [git](https://git-scm.com/)

#### Procedure
First make a clone of the master branch using the following command
```shell
git clone https://github.com/sayanadhikari/PPDyn.git
```
Then enter inside the *PPDyn* directory 
```shell
cd PPDyn
```
Now complile and built the *DSMCPy* code
```shell
make all
``` 
Usage
-----
Upon successful compilation, run the code using following command
```shell
ppdyn --i input.ini
```
Parameter Setup
----------------------
Edit the _input.ini_ and run the code again. The basic structure of _input.ini_ is provided below,
```ini
;
; @file		input.ini
; @brief	DSMCPy inputfile.
;
scope = default

scope = default


[simbox]
Lx  = 10.0    ; System length in X
Ly  = 10.0    ; System length in Y
Lz  = 10.0    ; System length in Z

[particles]
N     = 100     ; Number of particles
Vxmax = 1.0     ; Maximum velocity in X
Vymax = 1.0     ; Maximum velocity in Y
Vzmax = 1.0     ; Maximum velocity in Z
Temp  = 0.010   ;

[boundary]
btype = periodic ; Type of boundary

[time]
tmax  = 50.0    ; Final time
dt    = 0.010   ; time step size

[diagnostics]
dumpPeriod  = 10    ; Data dump period
useHDF5     = False ; Use HDF5 library to write data (Not implemented yet)

[options]
plotFigure    = False  ;
useNumba      = True  ;set to false to disable Numba
```
