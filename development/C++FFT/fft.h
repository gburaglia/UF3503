//
//
//  Fast Fourier Transform
//  Using Cooley-Tukey decimation-in-time radix-2 algorithm
//  Main information sources:
//  https://www.codeproject.com/Articles/9388/How-to-implement-the-FFT-algorithm
//  http://www.katjaas.nl/bitreversal/bitreversed.html
//  http://www.ti.com/lit/an/spna071a/spna071a.pdf
//  Main code source:
//  https://www.nayuki.io/page/free-small-fft-in-multiple-languages
//  Created by Gabriela Buraglia on 3/13/17.
//  Copyright © 2017 Gabriela Buraglia. All rights reserved.



//Bit Reversal
static size_t reverseBits(size_t x, unsigned int n);

//Computes the Discrete Fourier Transform (DFT) of the given vector
//Stores the result back into the vector
//The vector's length must be a power of 2
//Uses the Cooley-Tukey decimation-in-time radix-2 algorithm
void transformRadix2(std::vector<double> &real, std::vector<double> &imag);

//Testing Fft
static void testFft(int n);


//Filling vector with random numbers (-1,1)
static void randomReals(std::vector<double> &vec);

