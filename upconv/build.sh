gcc -o upconv -O2 -march=haswell -fopenmp -ftree-vectorize -lm -lfftw3 upconv.c fileio.c makepath.c -lm -lfftw3
