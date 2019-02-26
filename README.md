# FFT_Cooley-Tukey_LargePrimeField

I will try to implement Fast Fourier Transform using the Cooley-Tukey algorithm for a large prime field. <br>
I will give it a go with Python first, if all goes well I may try it out with Rust so that this can be a
learning opportunity for that as well.


*Input Conditions:* <br>
    a prime p of the order 256 or 512 bits. <br>
    a primitive 2^n(th) root of unity w and number n, <br>
    the set of evaluation points will be w,w^2,...,w^(2^n) <br>
    a 2^n vector of coefficients for a polynomial p <br>
    
*Compute:* p(w),...,p(w^(2^n))
