import ctypes
import numpy as np
import os
import time
import cantera as ct

def set_state_from_conservatives(elem_ID, quad_ID, Uq, gas=None):

	if gas is None:
		gas = ct.Solution('air.xml')

	# Get energy
	e = Uq[2] / Uq[0] - 0.5 * (Uq[1]**2 / Uq[0])
	# Get specific volume
	nu = 1./Uq[0]
	# Get YO2
	YO2 = Uq[3] / Uq[0]
	# Get YN2
	YN2 = Uq[4] / Uq[0]

	# gas_elem = physics.gas_elems[elem_ID, quad_ID]
	gas.UVY = e, nu, "O2:{},N2:{}".format(YO2, YN2)
	
	return gas

def get_pressure(Uq):
	gas = None
	ne = Uq.shape[0]
	nq = Uq.shape[1]

	P = np.zeros([ne, nq, 1])
	for ie in range(ne):
		for iq in range(nq):
			gas = set_state_from_conservatives(ie, iq, Uq[ie, iq])
			P[ie, iq] = gas.P

	return P


def main():
	ne = [10, 100, 1000]
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

			P = get_pressure(Uq)

			t1 = time.time()
			avg += (t1 - t0)
		avg /= 1
		print('Wall clock time = ' + str(avg) + 'seconds for ne = ' + str(ne_))
		# print(f"\nWall clock time = %g seconds" % avg + " for ne = " + str(ne_))
		# print("--------------------------------------------------------" + \
		# "-----------------------")
if __name__ == "__main__":
	main()
