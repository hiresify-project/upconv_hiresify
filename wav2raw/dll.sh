#!/bin/sh

gcc -Os -msse2 -fopenmp -ftree-vectorize -o wav2raw wav2raw.c ./../upconv/fileio.c ./../PLG_AUDIO_IO/PLG_AUDIO_IO.c
cp wav2raw.exe ../upconvfe/
