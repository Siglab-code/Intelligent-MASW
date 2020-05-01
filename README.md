# Intelligent MASW

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3776875.svg)](https://doi.org/10.5281/zenodo.3776875) 

## Install Required Library

The multichannel analysis of surface waves (MASW) forward solver is implemented in Fortran. The MASW inversion is done through Python Scipy library. 

```
$ pip install tensorflow
```

A smaller CPU-only package is also available:

```
$ pip install tensorflow-cpu
```
 
 

MASW inversion, as simple as this: 

Test = IMASW.masw(filename, initial_guess, lower_bound, upper_bound)   

Test.inverse() 


![Demo](/TRF.gif)

