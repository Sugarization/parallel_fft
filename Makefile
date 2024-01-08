CC = g++
EXEC = fft
FFT_OBJ = fft_parallel.o fft_base.o transpose.o
TEST_OBJ = timer.o tests.o
.PHONY: all clean

all: $(FFT_OBJ) $(TEST_OBJ) main.cc
	$(CC) $^ -O3 -fopenmp -std=c++17 -mfma  -o $(EXEC)

%.o: %.cc 
	$(CC) -c -O3 -std=c++17 $^ -o $@

fft_base.o:
	$(CC) -c fft_base.cc -std=c++17 -O3 -mfma -o fft_base.o

fft_parallel.o:
	$(CC) -c fft_parallel.cc -std=c++17 -O3 -mfma -fopenmp -march=native -o fft_parallel.o

transpose.o:
	$(CC) -c transpose.cc -std=c++17 -O3 -mfma -fopenmp -o transpose.o

clean:
	rm -f *.o $(EXEC)