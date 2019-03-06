# FFT_Cooley-Tukey_LargePrimeField aka NTT or Number Theoretic Transform

Fast Fourier Transform using the *Cooley-Tukey* algorithm for a large prime field. <br>
So far it is coded in Python, may try it out with Rust if time permits. <br>


**Input Conditions:** <br>
    a prime p of the order 256 or 512 bits. <br>
    a primitive 2^n(th) root of unity w and number n, <br>
    the set of evaluation points will be w,w^2,...,w^(2^n) <br>
    a 2^n vector of coefficients for a polynomial p <br>
    
**Compute:** p(w),...,p(w^(2^n)) <br>

**Test Case Details:** <br>
1st line test file is a 507-bit prime. <br>
2nd line test file is a 2^10th root of unity. <br>
3rd line test file is the value of n. <br>
After that there are 49 test cases. <br>
1st line of test case is the input vector. <br>
2nd line of test case is the target output. <br>
3rd line of test case is the 2^10th root of unity. <br>

**Implementations** <br>
1. **naiveNTT:** is the basic Cooley-Tukey NTT algorithm <br>
2. **modexpNTT**: is the modular exponentiation version of **naiveNTT** <br>
3. **NSquaredNTT**: is the most naive version without any optimization, runs in O(n^2) time <br>
