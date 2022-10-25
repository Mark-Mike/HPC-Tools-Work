#!/bin/bash

# gcc dgesv_gauss-elim.c -o dgesv_gauss-elim_gcc -lopenblas
LDLIBS="-lopenblas -fopenmp" make dgesv_gauss-elim_openmp_gcc CFLAGS="-march=native"
./dgesv_gauss-elim_openmp_gcc 100
