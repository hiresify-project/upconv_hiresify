#!/bin/sh

gcc -Os -lfftw3 -msse2 -fopenmp -ftree-vectorize -o upconv ../upconv.c ../fileio.c 
#gcc -O3 -march=core2 -pipe -fomit-frame-pointer -fopenmp -ftree-vectorize -o upconv_core2.exe ../upconv.c -llibfftw3-3
#gcc -Os -march=core2  -pipe -fomit-frame-pointer -fopenmp -ftree-vectorize -o upconv.exe ../upconv.c -llibfftw3-3
#gcc -O3 -march=athlon-xp -pipe -fomit-frame-pointer -fopenmp -ftree-vectorize -o upconv_athlon_xp.exe ../upconv.c -llibfftw3-3
#cp upconv.exe ../../upconvfe
