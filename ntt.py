import numpy as np
import random
import gmpy2
from gmpy2 import mpz
import math
from sympy import ntt as pac

hex_512_string = "F"*128
hex_256_string = "F"*64
MAX_512_integer = int(hex_512_string, 16)
MAX_256_integer = int(hex_256_string, 16)

# Combination of Carol Primes and Centered Traingular Primes
PRIMES = [409, 571, 631, 829, 1489, 1999, 2341, 2971, 3529, 4621, 4789,\
		 7039, 7669, 8779, 9721, 10459, 10711, 13681, 14851, 16069, 16381,\
		 17659, 20011, 20359, 23251, 25939, 27541, 29191, 29611, 31321, 34429,\
		 36739, 40099, 40591, 42589, 16127, 1046527, 16769023, 1073676287, \
		 68718952447, 274876858367, 4398042316799, 1125899839733759, 18014398241046527,]


def bigBitArrayGenerator(len_arr, p, start=0):
	"""
	The goal of this function is to fill a certain length of array with random
	integer from a minimum to maximum, where maximum is a prime number given
	as the second input.

	Parameters
	----------
	len_arr 	: number of elements wanted in array
	p			: a prime number
	start 		: minimum start for random integer; optional
	
	Returns
	-------
	array 		: array with random integers up untill prime number
	
	"""
	tmp = list()
	start = start
	k = 0
	lenTmp = 0
	while (lenTmp is not len_arr):
		k = random.randint(start, p)
		if len(tmp) == len_arr:
			break
		if k in tmp:
			lenTmp -= 1
		else:
			tmp.append(k)
		lenTmp += 1
		print (lenTmp)
	tmp.sort()
	return tmp


def conv_check_prime(n):
	"""
	A convenience method to check if the number is a prime.
	Returns True or False, based on if prime or not.

	Parameters
	----------
	n 		: integer
	
	Returns
	-------
	b 		: boolean
	
	"""
	return gmpy2.is_prime(n)

def conv_rou(k, prime):
	"""
	A convenience method to find out root of unity. If found more than one
	value returns a random value. With the current state of the function
	it will consider 1 as a valid root of unity contestant.

	Parameters
	----------
	k 		: A power of 2, the value of 2**n
	prime 	: prime number
	
	Returns
	-------
	k 		: root of unity
	
	"""
	rou = list()
	a = range(0, prime)
	for i in a:
		b = ( pow( mpz(i),mpz(k) ) % mpz(prime) == 1)
		if b:
			rou.append(i)
		else:
			pass
	if len(rou) == 1:
		return rou[0]
	else:
		print (rou)
		return random.choice(rou)

def calcWKS(k, w_k):
	wk_list = list()
	for i in range(k):
		wk_list.append(mpz(w_k)**mpz(i))
	return wk_list
def calcWKSI(k, w_k):
	wk_list = list()
	for i in range(k):
		wk_list.append(1//(mpz(w_k)**mpz(i)))
	return wk_list

def naiveNTT(p, k, w_k, A):
	"""
	A naive approach of the Cooley-Tukey algorithm, which is as provided in the
	lecture notes here:
	pdfs.semanticscholar.org/7a2a/4e7f8c21342a989d253704cedfb936bee7d7.pdf
	While Python 3 and above is capable of handling decently large integers,
	even then at one point overflow happens for big integer arithmetics. So I
	have used mpz type from gmpy2 library. This handles pretty much all kinds
	of big integer a 64 bit computer can compute in any decent time without any
	GPU capabilities. Please notice that, following the algorithm in the
	lecture I take input k instead of n. So, the input in place of n
	will be in the form of 2**n. Since, roots of unity is meant for complex
	numbers I have assumed that, what was actually meant here is a form of
	S = range(0 .. prime -1) and if (element mod prime)^k mod prime == 1
	then that element would be the root of unity.
	Also, this algorithm is more precisely known as Number Theoretic
	Transformation or NTT hence the reason for the naming instead of FFT.

	Parameters
	----------
	p 		: prime number (int)
	k 		: a power of 2 (int) / (2**n == k)
	w_k 	: kth root of unity (int)
	A 		: a int vecotr with max len of k ([int])
	
	Returns
	-------
	array 	: k length array
	
	"""
	if len(A) < k:
		raise ValueError("Length of vector(A) is not equal to \'k\'. Or if you\
			are entering as number \'n\'' then length of A is not equal to\
			2**n")
	if not isinstance(A, list):
		raise TypeError("Please enter a list of integers")
	if not np.all([isinstance(c, int) for c in A]):
		raise ValueError("Please make sure that the list (input A) only \
			contains integers")
	if not isinstance(p, int):
		raise ValueError("Please enter an integer for p")
	if not gmpy2.is_prime(p):
		raise ValueError("Please enter a prime number")
	if p == MAX_512_integer or p > MAX_512_integer:
		raise TypeError("Prime number is only supported up to 512 bits")
	if not isinstance(k, int):
		raise ValueError("Please enter an integer for k")



	retA = np.zeros(k)
	if k == 1:
		return A

	X_even = naiveNTT(p, k//2, w_k[0], A[::2])
	X_odd = naiveNTT(p, k//2, w_k[0], A[1::2])

	for i in range(k//2):

		retA[i] = (mpz( X_even[i] ) + ( ( w_k[i] ) * mpz( X_odd[i] ) )) % p
		
		retA[i+k//2] = (mpz( X_even[i] ) + ( ( w_k[i+k//2] ) * mpz( X_odd[i] ) )) % p
		
	return retA



if __name__ == '__main__':
	print ("Vector length 2")
	a = (calcWKS(2, 6))
	print (naiveNTT(7, 2**1, a, [2, 4]))
	print (naiveNTT(7, 2**1, a, [2, 6]))
	print ("Vector length 4")
	a = (calcWKS(4, 2))
	print (naiveNTT(11, 2**2, a, [2, 4, 6, 8]))
	a = (calcWKS(4, 2))
	print (naiveNTT(5, 2**2, a, [1, 4, 0, 0]))
	#test with sympy ntt
	print (pac([2,4], 7))
