# corr4py

Correlation functions speed up by Fortran

利用Fortran加速的互相关函数

## Dependencies: 
* FFTW3
* OpenMP
* numpy
* scipy(only for benchmark)

## Usage:

* `make build` for corrFFT.so file from fortran code.
* `python corr.py` for benchmark


## Matrix test:

### Random F64：1000 x 100 x 10000 :
value in (0,1)

|THREADS_NUM|fortran|corr4py|scipy(single thread)|
|-|-|-|-|
|2|32s| 42s| 88s|
|8|9s| 18s| -|
|16|6s| 14s| -|

Diff in sum: 0.0007627892564816459

### F32: 1000 x 100 x 10000 :
value in (0,1)


|THREADS_NUM|fortran|corr4py|scipy(single thread)|
|-|-|-|-|
|2|19s| 26s| 72s|
|8| 5s| 12s| - |
|16| 3s| 10s| - |

Diff in sum: 405401.25


corr4py take additional 10s(F64)/7s(F32) to transfer data from python to fortran.



