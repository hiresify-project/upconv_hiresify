#!/bin/sh

gcc -Os -msse2 -fopenmp -ftree-vectorize -o raw2wav raw2wav.c ./../upconv/fileio.c ./../PLG_AUDIO_IO/PLG_AUDIO_IO.c ../upconv/makepath.c
#cp raw2wav.exe ../upconvfe/
