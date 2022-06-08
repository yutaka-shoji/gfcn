# gfunc
G-function calculation library written in Fortran.
- Infinite Line Source Model (Ingersoll et al., 1954; Carslaw and Jaeger, 1992)
- Infinite Cylinder Source Model (Ingersoll et al., 1954; Carslaw and Jaeger, 1992)
- Composite-Medium Infinite Line Source Model (Li and Lai, 2012)

## Quick start
### On Google Colab
Just run quickstart.ipynb on Google Colab.  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yutaka-shoji/gfunc/blob/main/quickstart.ipynb)

## Local Use (Mac, Linux)
### Dependencies
Confirm the following dependencies are installed on your system.
- [gfortran](https://gcc.gnu.org/fortran/): GNU Fortran compiler
- [BLAS](http://www.netlib.org/blas/) & [LAPACK](http://www.netlib.org/lapack/): Mathematical computing library (recommend [OpenBLAS](https://github.com/xianyi/OpenBLAS))

### Clone
Git users may clone this repository with the following command:
```sh
git clone https://github.com/yutaka-shoji/gfunc.git
```

Or download this repository.

### Compile
Enter the directory and compile with following:
```sh
make all
```

### Validate
Python interface bundled.
```sh
python quickstart.py
```

## Acknowledgment
A part of the following Fortran library is included.
- [QUADPACK](http://www.netlib.org/quadpack/)
- [SLATEC](http://www.netlib.org/slatec/)
