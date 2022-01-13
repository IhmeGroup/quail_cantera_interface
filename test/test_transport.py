import ctypes
import numpy as np
import os
import time
import cantera as ct

LIB = ctypes.cdll.LoadLibrary('../build/libquail_cantera_interface.so')

GAS_CONSTANT = 8.3144621000e3 # [J / (K kmole)]

def get_W_from_Y(Wi, Y):
	Wi1 = 1./Wi
	return 1./np.dot(Wi1, Y)

def get_T_from_rhop(rho, P, W):
	return P / (rho * GAS_CONSTANT/W)

def get_P_from_rhoT(rho, T, W):
	return rho * GAS_CONSTANT/W * T

def get_gamma(cv, W):
	return (cv + GAS_CONSTANT / W) / cv;

def set_state_from_primitives(rho, T, u, Y, U):
	
	gas = ct.Solution('air_test.xml')

	U[:, :, 0] = rho
	U[:, :, 1] = rho*u
	U[:, :, 3] = rho * Y[0]

	W = get_W_from_Y(gas.molecular_weights, Y)
	P = get_P_from_rhoT(rho, T, W)
	gas.TPY = T, P, "O2:{},N2:{}".format(Y[0, 0], Y[1, 0])

	gamma = get_gamma(gas.cv, W)
	
	# Double check this definition
	U[:, :, 2] = rho * gas.UV[0] + 0.5*rho*u*u

	return U

def main():
	ne = [10, 100, 1000, 10000]
	nq = 3; ns = 5; nsp = 2

	for ie in range(len(ne)):
		ne_=ne[ie]
		avg = 0.0
		for j in range(1):
			t0 = time.time()

			Uq = np.zeros([ne_, nq, ns])

			# Uq[:, :, 0] = 1.
			# Uq[:, :, 1] = 0.
			# Uq[:, :, 2] = -5.140729783123304e+04
			# Uq[:, :, 3] = 0.21
			# Uq[:, :, 4] = 0.79
			rho = 0.3627
			u = 0.0
			T = 937.15

			Y = np.array([[0.21], [0.79]])

			set_state_from_primitives(rho, T, u, Y, Uq)

			mu = np.zeros([ne_, nq, 1])
			kap = np.zeros([ne_, nq, 1])
			Dk = np.zeros([ne_, nq, nsp])
			CANTERA_FILENAME = "air_test.xml"

			LIB.get_tranport_properties(
				ctypes.c_void_p(Uq.ctypes.data), 
				ctypes.c_void_p(mu.ctypes.data),
				ctypes.c_void_p(kap.ctypes.data),
				ctypes.c_void_p(Dk.ctypes.data),
				ctypes.c_int(ne_), 
				ctypes.c_int(nq), 
				ctypes.c_int(ns),
				ctypes.c_int(nsp),
				1,
				ctypes.c_char_p(CANTERA_FILENAME.encode('utf-8'))
				)
			import code; code.interact(local=locals())
			t1 = time.time()
			avg += (t1 - t0)
		avg /= 1
		# print('Wall clock time = ' + str(avg) + 'seconds for ne = ' + str(ne_))

		print(f"\nWall clock time = %g seconds" % avg + " for ne = " + str(ne_))
		print("--------------------------------------------------------" + \
		"-----------------------")
if __name__ == "__main__":
	main()
