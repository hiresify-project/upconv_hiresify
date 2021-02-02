gcc -O0 -march=haswell -msse2 -fopenmp -ftree-vectorize -o raw2wav raw2wav.c ./../upconv/fileio.c ./../PLG_AUDIO_IO/PLG_AUDIO_IO.c -lm -static -lstdc++ -lgcc -lwinpthread -Wl,--stack,10485760
