import ctypes
import numpy as np
import os
import time

LIB = ctypes.cdll.LoadLibrary('./build/libcombustor.so')

def main():
	ne = [10, 100, 1000, 10000]
	nq = 3; ns = 5;

	for ie in range(len(ne)):
		ne_=ne[ie]
		avg = 0.0
		for j in range(1):
			t0 = time.time()

			Uq = np.zeros([ne_, nq, ns])

			Uq[:, :, 0] = 1.
			Uq[:, :, 1] = 0.
			Uq[:, :, 2] = -5.140729783123304e+04
			Uq[:, :, 3] = 0.21
			Uq[:, :, 4] = 0.79

			P = np.zeros([ne_, nq, 1])

			# lib = ctypes.cdll.LoadLibrary('./build/libcombustor.so')
			LIB.get_pressure_multithread(
					ctypes.c_void_p(Uq.ctypes.data), 
					ctypes.c_void_p(P.ctypes.data),
					ctypes.c_int(ne_), 
					ctypes.c_int(nq), 
					ctypes.c_int(ns),
						)
			print(P)
			t1 = time.time()
			avg += (t1 - t0)
		avg /= 1
		# print('Wall clock time = ' + str(avg) + 'seconds for ne = ' + str(ne_))

		print(f"\nWall clock time = %g seconds" % avg + " for ne = " + str(ne_))
		print("--------------------------------------------------------" + \
		"-----------------------")
if __name__ == "__main__":
	main()
