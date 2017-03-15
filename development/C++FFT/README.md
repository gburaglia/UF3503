# FFT in C++

## What is FFT?
* Optimized computational algorithm to implement the Discreet Fourier Transform
* Factorizes DFT matrix into a product of
* Many was to implement


### This specific implementation:
##### Using Cooley-Tukey decimation-in-time radix-2 FFT
* Simplest and most common form of the Cooley–Tukey algorithm
* Divides DFT of size N into two interleaved DFTs of size N/2
  * By Danielson-Lancoz lemma
* 1st compute even-indexed inputs
* 2nd compute odd-indexed inputs
* 3rd combines 2 results
* Performed recursively reduces runtime to O(N log N)
  * These subsets then divided into their subsets
  * Continues until 2 members per subset
  * This is computational result of the Danielson-Lancoz lemma O(N2) -> O(N log N)

Limitation: Array size must equal 2^x

Solution: If it's not fill with next size and fill remaining spaces with 0s

##### Code breakdown
1. Make sure data size is 2^x
2. Reverse bits
3. Compute "twiddle factors" with Trigonometric tables and store for future use
4. Carry out Cooley-Tukey algorithm explained above

### Further definitions:

** Twiddle factor ** :
* In FFT algorithms, refers to trigonometric constant coefficients that are multiplied by the data in the course of the algorithm ![Figure 1](Images/Factors.png?raw=true)
* Step 3 in

** Bit Reversal ** :
* Makes the mathematical calculations of the second part easier.
![Figure 1](Images/BitReversal2.png?raw=true)

 ![Figure 1](Images/BitReversal.png?raw=true)




** Extra Picture Explanations for Cooley-Tukey ** :
* Splitting
 ![Figure 1](Images/Splitting.png?raw=true)
* Final Step
![Figure 1](Images/PseudoCode.png?raw=true)
