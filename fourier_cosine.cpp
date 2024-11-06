//*****************************************************************************************************************
//
//   fourier_cosine.cpp:
//
//    This is an example program that calculates the FFT of a cosine function
//  using the FFTW package. The goal of this to compare the numerical solution
//  with the analytical well known solution. We selected a cosine function with
//  a frequency of 3, which means that we expect a peak at i=3 and at i=N-3.
//   We observe those expected peaks and also we calculate the inverse FFT to 
//  compare that result with the analytical cosine function, getting a perfect match
//  as expected.
//  
//	flags 
//      g++ -o fourier_cosine.out fourier_cosine.cpp -lfftw3 -lm 
//      g++ -o fourier_cosine.out fourier_cosine.cpp -lfftw3 -lm && ./fourier_cosine.out 
//
//          By Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (11/06/2024)
//
//*****************************************************************************************************************
//*****************************************************************************************************************

#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <fstream>
#define N 16

using namespace std;

int main(void) {
    fftw_complex in[N], out[N], in2[N]; /* double [2] */
    fftw_plan p, q;
    int i;

    /* prepare a cosine wave */
    for (i = 0; i < N; i++) {
        in[i][0] = cos(3 * 2*M_PI*i/N);
        in[i][1] = 0;
    }

    /* forward Fourier transform, save the result in 'out' */
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Write the results (frequency, real part, imaginary part) to text files
    ofstream outfile_re("fftw_output_re.txt");
    ofstream outfile_im("fftw_output_im.txt");

    for (i = 0; i < N; i++){
        outfile_re << i << ", " << out[i][0] << '\n';  // Real part
        outfile_im << i << ", " << out[i][1] << '\n';  // Imaginary part
    }
    fftw_destroy_plan(p);

    /* backward Fourier transform, save the result in 'in2' */
    printf("\nInverse transform:\n");
    q = fftw_plan_dft_1d(N, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(q);
    /* normalize */
    for (i = 0; i < N; i++) {
        in2[i][0] *= 1./N;
        in2[i][1] *= 1./N;
    }
    for (i = 0; i < N; i++)
        printf("recover: %3d %+9.5f %+9.5f I vs. %+9.5f %+9.5f I\n",
        i, in[i][0], in[i][1], in2[i][0], in2[i][1]);
    fftw_destroy_plan(q);

    fftw_cleanup();
    return 0;
}