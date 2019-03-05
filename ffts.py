import numpy as np
import timeit


# np.shape			--> gives how length of array if 1D array, comes as tuple
# np.arange			--> creates an array from (0, n-1)
# np.reshaped		--> (x, y): y is how many values are there in the smallest
#						array; x = k/y, given k is length or the array to be
#						reshaped, so x is how many small arrays in big array
def DFT_slow(x):
	x = np.asarray(x, dtype=float)
	N = x.shape[0]
	n = np.arange(N)
	k = n.reshape((N, 1))
	M = np.exp(-2j * np.pi * k * n / N)
	return np.dot(M, x)

def FFT(prime_number = 7, n = 1, root_of_unity = 6, vector_coef = [2, 4]):
    """A recursive implementation of the 1D Cooley-Tukey FFT"""
    x = np.asarray(vector_coef, dtype=float)
    N = x.shape[0]
    
    if N % 2 > 0:
        raise ValueError("size of x must be a power of 2")
    elif N <= 32:  # this cutoff should be optimized
        return DFT_slow(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j * np.pi * np.arange(N) / N)
        return np.concatenate([X_even + factor[:N / 2] * X_odd,
                               X_even + factor[N / 2:] * X_odd])


def FFT_Vectorized(x):
	x = np.asarray(x, dtype=float)
	N = x.shape[0]

	if np.log2(N) % 1 > 0:
		raise ValueError("size of x must be a power of 2")
	
	N_min = min(N, 32)

	n = np.arange(N_min)
	k = n[:, None]
	M = np.exp(-2j * np.pi * n * k / N_min)
	X = np.dot(M, x.reshape((N_min, -1)))

	# X is now a 2D array
	while X.shape[0] < N:
		X_even = X[:, :X.shape[1]/2]
		X_odd = X[:, X.shape[1]/2:]
		factor = np.exp(-1j * np.pi * np.arange(X.shape[0])
						 / X.shape[0])[:, None]
		X = np.vstack([X_even + factor * X_odd,
						X_even - factor * X_odd])

	return X.ravel() 


x = ([10, 12, 14, 16])
# y = [13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095, 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084094, 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095**2, 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095**3, 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095**4, 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095**5, 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095**6]
print np.allclose(FFT_Vectorized(x), np.fft.fft(x))
# print np.allclose(np.fft.fft(y))

