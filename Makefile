CC = g++
all:
	$(CC) main.cc fft_base.cc fft_parallel.cc -DNDEBUG -fopenmp -std=c++17 -O3 -o fft
debug:
	$(CC) main.cc fft_base.cc fft_parallel.cc -fopenmp -std=c++17 -O2 -g -o fft