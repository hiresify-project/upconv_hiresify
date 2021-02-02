gcc -O0 -msse2 -march=native -fopenmp -ftree-vectorize -o wav2raw wav2raw.c ./../upconv/fileio.c ./../PLG_AUDIO_IO/PLG_AUDIO_IO.c ../upconv/makepath.c
