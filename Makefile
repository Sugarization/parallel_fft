CC = g++
OBJ = fft_base.o fft_para.o timer.o
.PHONY: all debug

all: timer fft
	$(CC) main.cc -fopenmp $(OBJ) -o fft

timer:
	$(CC) timer.cc -c -o timer.o

fft: fft_base fft_parallel

fft_parallel: 
	$(CC) -c fft_parallel.cc -DNDEBUG -fopenmp -std=c++17 -g -mfma -march=native -O3 -o fft_base.o

fft_base:
	$(CC) -c fft_base.cc -std=c++17 -g -O3 -mfma -o fft_para.o

check: timer fft 
	$(CC) main.cc -DFFT_CHECK $(OBJ) -o fft

clean:
	rm *.o fft