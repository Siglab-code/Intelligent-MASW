# Intelligent MASW

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3776875.svg)](https://doi.org/10.5281/zenodo.3776875) 

## Install 
The multichannel analysis of surface waves (MASW) forward solver is implemented in Fortran. The MASW inversion is done through Python Scipy library. The attached wrapper for Fortran code can only run in Linux system. For the Windows user, the Fortran code nees to be complied throough Numpy f2py function. The following dependencies are required: 

gfortran compiler:
```
$ sudo apt install gfortran-9
```

blas and lapack library: 

```
$ sudo apt-get install libblas-dev checkinstall 
$ sudo apt-get install libblas-doc checkinstall 
$ sudo apt-get install liblapacke-dev checkinstall 
$ sudo apt-get install liblapack-doc checkinstall
```

wrapper function: 
```
$ f2py -m name -c fortran.f90
```



code demonstration:  

```
$ Test = IMASW.masw(filename, initial_guess, lower_bound, upper_bound)    
$ Test.inverse() 
```

![Demo](/TRR.gif)

