#include "tests.h"

int main(int argc, char **argv)
{
    int N = 1 << 10;
    if (argc >= 2) {
        sscanf(argv[1], "%d", &N);
    }
    
    sanityCheck(FFT_Type::cooley, FFT_Type::embed, N, 1);
    return 0;
}
