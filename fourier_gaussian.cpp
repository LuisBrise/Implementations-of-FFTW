//*****************************************************************************************************************
//
//   fourier_gaussian.cpp:
//
//    This is an example program that calculates the FFT of a gaussian function
//  using the FFTW package. The goal of this to compare the numerical solution
//  with the analytical well known solution. We implemented subrutine that saves
//  the result in some .txt files.
//  
//	flags 
//      g++ -o fourier_gaussian.out fourier_gaussian.cpp -lfftw3 -lm 
//      g++ -o fourier_gaussian.out fourier_gaussian.cpp -lfftw3 -lm && ./fourier_gaussian.out 
//
//          By Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (11/06/2024)
//
//*****************************************************************************************************************
//*****************************************************************************************************************

#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <fstream>

using namespace std;

// Define a Gaussian function in time domain
double gaussian(double t, double t0, double sigma) {
    return exp(-0.5 * pow((t - t0) / sigma, 2.0));
}

int main() {
    const int n = pow(2,6);  // Number of points (must be a power of 2)
    const double t_min = -10.0;  // Minimum time
    const double t_max = 10.0;   // Maximum time
    const double delta_t = (t_max - t_min) / n;  // Time step
    const double sigma = 1.0;  // Width of the Gaussian
    const double t0 = 0.0;     // Center of the Gaussian

    // Arrays to hold the input and output data
    fftw_complex* fft_input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n); 
    fftw_complex* fft_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    // Initialize the time-domain data with the Gaussian function
    for (int i = 0; i < n; ++i) {
        double t = t_min + i * delta_t;
        fft_input[i][0] = gaussian(t, t0, sigma); //real part of the input
        fft_input[i][1] = 0;                      //imaginary part
    }

    // Plan for FFT
    // One may use "fftw_plan_dft_r2c_1d" for real-to-complex transforms instead
    // You'd get N / 2 + 1 complex outputs for N real inputs 
    // (the redundant symmetric outputs are not generated)
    fftw_plan plan = fftw_plan_dft_1d(n, fft_input, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);

    // Perform the FFT
    fftw_execute(plan);

    // Create frequency array
    const double df = 1.0 / (n * delta_t);  // Frequency resolution
    double freqs[n];  // Only n/2 + 1 frequencies are returned by FFTW (real-to-complex)
    for (int i = 0; i < n + 1; ++i) {
        freqs[i] = -1.0 / (2.0 * delta_t) + i * df;
    }

    // Write the results (frequency, real part, imaginary part) to text files
    ofstream outfile_re("fftw_output_re.txt");
    ofstream outfile_im("fftw_output_im.txt");

    if (outfile_re.is_open() && outfile_im.is_open()) {
        for (int i = 0; i < n; ++i) {
            outfile_re << i+1 << ", " << fft_result[i][0] << '\n';  // Real part
            outfile_im << freqs[i] << ", " << fft_result[i][1] << '\n';  // Imaginary part
        }
        outfile_re.close();
        outfile_im.close();
        cout << "FFTW output written to fftw_output_re.txt and fftw_output_im.txt" << endl;
    } else {
        cerr << "Unable to open files for writing." << endl;
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(fft_result);

    return 0;
}

