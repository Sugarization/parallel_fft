CC = g++
all:
	$(CC) main.cc fft_base.cc fft_parallel.cc -DNDEBUG -fopenmp -std=c++17 -g -mfma  -march=native -O3 -o fft
debug:
	$(CC) main.cc fft_base.cc fft_parallel.cc -fopenmp -std=c++17 -mfma -O1 -g -o fft