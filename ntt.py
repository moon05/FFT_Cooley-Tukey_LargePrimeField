import numpy as np
import random
import gmpy2
from gmpy2 import mpz
import math
from sympy import ntt as pac
import timeit
from itertools import chain
from collections import Counter

hex_512_string = "F"*128
hex_256_string = "F"*64
MAX_512_integer = int(hex_512_string, 16)
MAX_256_integer = int(hex_256_string, 16)

BIG_512_bit_PRIME = 1448710878603695229409412310732144622848441888753223448575120702981577542281415018652993561619563625768081760802164016332664868818287202677226268571177357443196189721711411095564550863916701456326047703907558905556769006452749568850893726964097212712534033673887482435571016493398875714661616695973805482073975013312603832180072496240990965064767802796335226649504512522482656310868574294571090894743593101737100996294180550347627991120402954521019478319594965583602635139406639554109668397505671888878121704968823589493631448107907114217361283249827317880226419193187177259109633935404301712464708823348328976497742932566121346081365525408595640468945290937758700039629256982005716182859284737841217145224064176914658186450842162070805962848366278043048270579195785338465402255311543317952316743397150239782344358067398630187697581974758378383164399938534007082923964046532198511555664094030236510446638373837514965653315445276132434900279652476100332956459885822293290144417509773544510283867826006224349987445172866722671938741506492171575030480586949013381712605577261029640924163986974030623538195577104997915377544237531482486437582261243814911374465420432886979360185662951101417810621647359742220813855302270079227324459422907061626612236531675332230997981115552142564940071988563230306482080972374909297566710962638028801
TEST_CASE_PRIME = 0
TEST_CASE_UNITY = 0
TEST_CASE_K = 1024
TEST_CASES = dict()
TEST_ANSWERS = dict()
# Combination of Carol Primes and Centered Traingular Primes
PRIMES = [409, 571, 631, 829, 1489, 1999, 2341, 2971, 3529, 4621, 4789,\
		 7039, 7669, 8779, 9721, 10459, 10711, 13681, 14851, 16069, 16381,\
		 17659, 20011, 20359, 23251, 25939, 27541, 29191, 29611, 31321, 34429,\
		 36739, 40099, 40591, 42589, 16127, 1046527, 16769023, 1073676287, \
		 68718952447, 274876858367, 4398042316799, 1125899839733759, 18014398241046527,]



###############################################################################
####################### Test File Processing Functions ########################
###############################################################################

def hasNewLineInList(l):
	for i in l:
		if "\n" in i:
			return True
	return False

def convertStringsToIntsInList(l):
	j = l.split()
	for i in range(len(j)):
		j[i] = mpz(int(j[i]))
	return j

def processCASE(strSuperList):
	global TEST_CASES
	global TEST_ANSWERS
	s = "CASE_"
	try:
		assert len(strSuperList%2) == 0
	except Exception:
		print("Something is wrong with the test file")
	for i in range(0,len(strSuperList),3):
		if i is 0:
			s+= str(i)
		else:
			s+=str((i+1)//3)
		if s not in TEST_CASES.keys():
			TEST_CASES[s] = convertStringsToIntsInList(strSuperList[i])
		if s not in TEST_ANSWERS.keys():
			TEST_ANSWERS[s] = convertStringsToIntsInList(strSuperList[i+1])
		s = "CASE_"
	return




###############################################################################
############################ NTT Helper Functions #############################
###############################################################################

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


def check_prime(n):
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

def generate_rou(k, prime):
	"""
	A convenience method to find out roots of unity. If found more than one
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
	rou = [[] for c in range(k)]
	a = range(0, prime)
	for i in a:
		print ("prime "+ str(i))
		for j in range(k):
			print ("power" + str(k))
			b = ( (mpz(i)**mpz(j)) % mpz(prime) == 1)
			if b:
				rou[j].append(i)
			else:
				pass
	counts = Counter(chain.from_iterable(rou))
	uniques = [k for k, c in counts.items() if c==1]
	print (uniques)
	return uniques

def calcWKS(k, w_k):
	"""
	Helper method to pre computer all powers of roots of unity

	Parameters
	----------
	k 		: A power of 2, the value of 2**n
	prime 	: n_th root of unity
	
	Returns
	-------
	k 		: vector of range (w_k**0 .. w_k**(k-1))
	
	"""
	wk_list = list()
	for i in range(k):
		wk_list.append(mpz(w_k)**mpz(i))
	return wk_list



###############################################################################
############################ Main NTT Function/s ##############################
###############################################################################


def naiveNTT(p, k, w_k, A):
	"""
	A naive approach of the Cooley-Tukey algorithm, which is as provided in the
	lecture notes here:
	pdfs.semanticscholar.org/7a2a/4e7f8c21342a989d253704cedfb936bee7d7.pdf
	While Python 3 and above is capable of handling decently large integers,
	even then at a certain point overflow happens for big integer arithmetics.
	So I have used mpz from gmpy2 library. This handles all kinds of big
	integer a 64 bit computer can compute in any decent time without any
	GPU capabilities. Please notice that, following the algorithm in the
	lecture I take input k instead of n. So the value of k will be in the form
	of 2**n.

	Also, this algorithm is more precisely known as Number Theoretic
	Transformation or NTT hence the naming, instead of FFT.

	Parameters
	----------
	p 		: prime number (int)
	k 		: a power of 2 (int) / (2**n == k)
	w_k 	: kth root of unity (int)
	A 		: a int vector with k elements; each element is a type: mpz(int)
	
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
	if not gmpy2.is_prime(p):
		raise ValueError("Please enter a prime number")
	if p == MAX_512_integer or p > MAX_512_integer:
		raise TypeError("Prime number is only supported up to 512 bits")
	if not isinstance(k, int):
		raise ValueError("Please enter an integer for k")



	retA = np.zeros(k)
	if k == 1:
		return A

	X_even = naiveNTT(p, k//2, w_k[0]**2, A[::2])
	X_odd = naiveNTT(p, k//2, w_k[0]**2, A[1::2])

	for i in range(k//2):

		retA[i] = (mpz( X_even[i] ) + ( ( w_k[i] ) * mpz( X_odd[i] ) )) % p
		
		retA[i+k//2] = (mpz( X_even[i] ) + ( ( w_k[i+k//2] ) * mpz( X_odd[i] ) )) % p
		
	return retA

def modexpNTT(p, k, w_k, A):
	retA = np.zeros(k)
	if k == 1:
		return A

	X_even = modexpNTT(p, k//2, w_k[0]**2, A[::2])
	X_odd = modexpNTT(p, k//2, w_k[0]**2, A[1::2])

	for i in range(k//2):

		retA[i] = pow((mpz( X_even[i] ) + ( ( w_k[i] ) * mpz( X_odd[i] ) )), 1, p)
		
		retA[i+k//2] = pow((mpz( X_even[i] ) + ( ( w_k[i+k//2] ) * mpz( X_odd[i] ) )), 1, p)
		
	return retA

def NSquaredNTT(p, k, w_k, A):
	W = [[] for c in range(k)]
	for i in range(k):
		tmp = []
		for j in range(k):
			tmp.append((2**(i*j))%p)
		W[i] = tmp
	retA = np.zeros(k)
	for j in range(len(A)):
		for i in range(k):
			retA[j] += (W[j][i]*A[i]) % p

	return retA

###############################################################################
############################# Timing Functions ################################
###############################################################################

def func_naive():
	a = (calcWKS(4, 2))
	naiveNTT(5, 2**2, a, [1, 4, 0, 0])
def func_modexp():
	a = (calcWKS(4, 2))
	modexpNTT(5, 2**2, a, [1, 4, 0, 0])
def func_nsquared():
	a = (calcWKS(4, 2))
	NSquaredNTT(5, 2**2, 2, [1, 4, 0, 0])
def func_calcWKS():
	calcWKS(TEST_CASE_K, TEST_CASE_UNITY)



###############################################################################
############################## Read Test File #################################
###############################################################################


def readTestFile(filename):
		with open("./"+filename) as fp:
			allLINES = fp.read().splitlines()
			totLines = 0
			for line in fp:
				totLines += 1
			global TEST_CASE_PRIME
			global TEST_CASE_UNITY
			TEST_CASE_PRIME = mpz(int(allLINES.pop(0)))
			TEST_CASE_UNITY = mpz(int(allLINES.pop(0)))
			allLINES.pop(0)
			print ("Length of file before processing {}".format(len(allLINES)))
			processCASE(allLINES)
			print ("Processing done")
			return


if __name__ == '__main__':
	print ("Vector length 2")
	# a = (calcWKS(2, 6))
	# print (naiveNTT(7, 2**1, a, [2, 4]))
	# print (naiveNTT(7, 2**1, a, [2, 6]))
	# print ("Vector length 4")
	# a = (calcWKS(4, 2))
	a = (calcWKS(4, 2))
	# print (naiveNTT(5, 2**2, a, [1, 4, 0, 0]))
	#test with sympy ntt
	# print (pac([2,4], 7))
	# print (NSquaredNTT(5, 2**2, 2, [1, 4, 0, 0]))
	# print (conv_rou(20, 159))
	# print (timeit.timeit(func_naive, number=10000))
	# print (timeit.timeit(func_pow, number=10000))
	# print (timeit.timeit(func_nsquared, number=10000))
	print (readTestFile("test.txt"))
	# print (TEST_CASE_UNITY)
	# print (b)
	# print (TEST_CASES.keys())
	# print (TEST_ANSWERS.keys())
	# print (TEST_CASES["CASE_48"])
	# print (len(CASE_1))
	print ("STARTED")
	b = calcWKS(TEST_CASE_K, TEST_CASE_UNITY)
	print (timeit.timeit(func_calcWKS, number=10000))	
	naiveNTT(TEST_CASE_PRIME, TEST_CASE_K, b, TEST_CASES["CASE_0"])
	print ("DONE")
