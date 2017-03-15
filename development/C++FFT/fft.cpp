
#include <cmath>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <iostream>
#include "fft.h"

using std::size_t;
using std::vector;


static size_t reverseBits(size_t x, unsigned int n) {
    size_t result = 0;
    unsigned int i;
    
    //x & 1 is bitwise AND operation
    //if the last bit in x is 1, result is 1 otherwise 0
    
    //x>>= 1 shift 1 bit to the right and assign new value to x
    //x>>=1 is x/(2^1)
    
    //Note: logical bit shifting rounds down to lowest integer
    
    //result<<1 shift 1 bit to the right
    //result<<1 is result*(2^1)
    
    //note that the indexes start from 0 which means the
    //real part of the complex is on the even-indexes
    //and the complex part is on the odd-indexes
   
    for (i = 0; i < n; i++, x >>= 1)
        result = (result << 1) | (x & 1);
    return result;
}

void transformRadix2(vector<double> &real, vector<double> &imag) {
 
    if (real.size() != imag.size())
        throw "Mismatched lengths";
    
    //n is size of real vector, which should be same as imaginary vector
    size_t n = real.size();
    
    unsigned int levels;
    //Compute levels = floor(log2(n))
    //floor(log2(n)) is largest int <= (log2(n))
    //i.e. 2^x, where x would be levels
    
    {
        size_t temp = n;
        levels = 0;
        while (temp > 1) {
            levels++;
            //temp>>=1
            //shift one bit to the right and assigns new value to temp
            //temp>>=1 is temp/(2^1)
            temp >>= 1;
    }
        
        //1u<<levels is 1*(2^levels)
        //Compares result to n
        //If not =, it means length was not a power of 2
        if ((1u << levels) != n)
            throw "Length is not a power of 2";
    }
    
    //...so far just verifying that vector length is a power of 2
    
    
    // Trignometric tables
    //Calculates "twiddle factors"
    //e^(-2*i*pi/N)
    vector<double> cosTable(n / 2);
    vector<double> sinTable(n / 2);
    for (size_t i = 0; i < n / 2; i++) {
        cosTable[i] = cos(2 * M_PI * i / n);
        sinTable[i] = sin(2 * M_PI * i / n);
    }
    
    // Bit-reversed addressing permutation
    for (size_t i = 0; i < n; i++) {
        //i goes from 0 to n
        //levels is power of exponent
        size_t j = reverseBits(i, levels);
        if (j > i) {
            //swap the real part
            double temp = real[i];
            real[i] = real[j];
            real[j] = temp;
            //swap the complex part
            temp = imag[i];
            imag[i] = imag[j];
            imag[j] = temp;
        }
    }
    
    // Cooley-Tukey decimation-in-time radix-2 FFT
    //the input is in bit-reversed order (hence "decimation-in-time")
    //tp is twiddle product
    //tpre-for real, tpim-for imaginary
    for (size_t size = 2; size <= n; size *= 2) {
        size_t halfsize = size / 2;
        size_t tablestep = n / size;
        for (size_t i = 0; i < n; i += size) {
            for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
                double tpre =  real[j+halfsize] * cosTable[k] + imag[j+halfsize] * sinTable[k];
                double tpim = -real[j+halfsize] * sinTable[k] + imag[j+halfsize] * cosTable[k];
            
                //odd index calcuations
                real[j + halfsize] = real[j] - tpre;
                imag[j + halfsize] = imag[j] - tpim;
                
                //even index calculations
                real[j] += tpre;
                imag[j] += tpim;
            }
        }
        // Prevent overflow in 'size *= 2'
        if (size == n)
            break;
    }
    
   
}

int main()
{
        testFft(8);
}

static void testFft(int n)
{
    vector<double> inputreal(n);
    vector<double> inputimag(n);
    randomReals(inputreal);
    //Commented out here, asssuming there is no imaginary input portion
    //randomReals(inputimag);
    
    vector<double> actualoutreal(inputreal);
    vector<double> actualoutimag(inputimag);
    transformRadix2(actualoutreal, actualoutimag);
    
    //Size
    std::cout << "size =" << std::setw(2) << std::setfill(' ') << n<< std::endl<< std::endl;
    
    //Column titles
    std::cout<< std::left << std::setw(10)<<std::setfill(' ')<< "IN-REAL" << std::left << std::setw(10) <<std::setfill(' ')<< "IN-IMG" << std::left << std::setw(10) <<std::setfill(' ')<<"OUT-REAL"<< std::left << std::setw(10)<<std::setfill(' ')<<"OUT-IMG"<<std::endl;

    //Printing out original and transformed data
    for (int i = 0; i < n; i++)
    {
        std::cout<< std::left << std::setw(10)<<std::setfill(' ')<< inputreal[i] << std::left << std::setw(10) <<std::setfill(' ')<< inputimag[i]<< std::left << std::setw(10) <<std::setfill(' ')<<actualoutreal[i]<< std::left << std::setw(10)<<std::setfill(' ')<<actualoutimag[i]<<std::endl;
    }
    
    
}

static void randomReals(vector<double> &vec)
{
    //using iterators, which are limited pointers
    for (vector<double>::iterator it = vec.begin(); it != vec.end(); ++it)
        //*it dereferences iterator
        //assignment inserts value in place iterator points to
        
        //returns from (-1, 1)
        *it = (rand() / (RAND_MAX + 1.0)) * 2 - 1;
}


