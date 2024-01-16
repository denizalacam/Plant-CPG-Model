#! /usr/bin/python

import ctypes as ct
import numpy as np
import tools as tl
import os

CUDA_ENABLED = False
NOISE_ENABLED = False

lib = ct.cdll.LoadLibrary(os.path.dirname(__file__)+'/lib/_plant.so')
PI2 = tl.PI2


# set hard:
V_0 = 0.                # coupling threshold
threshold_slope = 100.          # coupling threshold slope

N_EQ1 = 5
N_EQ3 = 3*N_EQ1
N_EQ4 = 4*N_EQ1
NP = 8


INITIAL_ORBIT = np.array([-0.62376542, 0.00650901])
dt = 0.05
stride = 50
N_integrate = 2*10**5
IDX_THRESHOLD = 0
THRESHOLD = 0.

def single_orbit(parameterArray, DT_ORBIT=dt, N_ORBIT=N_integrate, STRIDE_ORBIT=stride, V_threshold=THRESHOLD, verbose=0):
    y_initial = INITIAL_ORBIT.copy()    #passing INITIAL_ORBIT directly to integrator or copying just the
                                        # reference changes it's state at the end
#    print 'initial orbit: ', INITIAL_ORBIT, '; y_initial: ', y_initial
                                #TODO integrator should take a const value as initial state array
    X = integrate_one_rk4(parameterArray, y_initial,DT_ORBIT/float(STRIDE_ORBIT), N_ORBIT, STRIDE_ORBIT )
    x_raw, y = X[0], X[1]
    x_m, y_m = tl.splineLS1D(), tl.splineLS1D()

    try:
        ni = np.asarray(tl.crossings(x_raw, V_threshold), dtype=int) # convert to millivolts
        x, y = x_raw[ni[-2]:ni[-1]], y[ni[-2]:ni[-1]]
        t = tl.PI2*np.arange(x.size)/float(x.size-1)
        #compute smoothened curves of V,x for the values between the last two crossings
        x_m.makeModel(x, t); y_m.makeModel(y, t)

    except:
        print '# single_orbit:  No closed orbit found!'
        raise ValueError

    T = DT_ORBIT*x.size	 # in msec.

    return x_m, y_m, T

lib.integrate_one_rk4.argtypes = [ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.POINTER(ct.c_double),
					ct.c_double, ct.c_uint, ct.c_uint]
def integrate_one_rk4(parameterArray, initial_state, dt, N_integrate, stride=42):
	initial_state = np.asarray(initial_state)
	assert initial_state.size == N_EQ1

	X_out = np.zeros((N_EQ1*N_integrate), float)

	lib.integrate_one_rk4(initial_state.ctypes.data_as(ct.POINTER(ct.c_double)),
				parameterArray.ctypes.data_as(ct.POINTER(ct.c_double)),
				X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
				ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride))
	return np.reshape(X_out, (N_EQ1, N_integrate), 'F')

#computes traces and phase differences in cpp, returns both
lib.integrate_n_rk4_phasedifferences_electricExcitatoryAndInhibitoryConnections.argtypes = [ct.c_uint, ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER((ct.c_uint)),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_int),
                                                 ct.POINTER(ct.c_int),
                                                 ct.POINTER(ct.c_double),
                                                 ct.c_double, ct.c_uint, ct.c_uint, ct.c_uint]
def integrate_n_rk4_phasedifferences_electricExcitatoryAndInhibitoryConnections(networkSize, initial_states,
                                     coupling_strengths_chlorine, coupling_strengths_potassium,
                                     coupling_strengths_sodium, coupling_strengths_electric,
                                     dt, N_integrate, stride=1, parameterSet=None):
    initial_states = np.asarray(initial_states) #
    assert initial_states.size == N_EQ1*networkSize

    coupling_strengths_chlorine = np.asarray(coupling_strengths_chlorine)
    assert coupling_strengths_chlorine.size == networkSize * networkSize
    coupling_strengths_potassium = np.asarray(coupling_strengths_potassium)
    assert coupling_strengths_potassium.size == networkSize * networkSize
    coupling_strengths_sodium = np.asarray(coupling_strengths_sodium)
    assert coupling_strengths_sodium.size == networkSize * networkSize
    coupling_strengths_electric = np.asarray(coupling_strengths_electric)
    assert coupling_strengths_electric.size == networkSize * networkSize

    #for each of cells 1 to n-1, phase differences with respect to cell 0 is computed at all times when crossings occur in Cell 0
    #actual numberOfCrossings that are used in computation of phase differences is stored via 'numberOfCrossings' (after ignoring the crossings
    # at the beginning and at the end that can not be used for phase differences
    assert parameterSet.size == NP*networkSize
    X_out = np.zeros((networkSize*N_integrate), float)
    P_out = np.zeros(((networkSize-1)*(N_integrate/2+1)), float) #Maximum possible size
    numberOfPhaseDifferences = np.zeros(1,int);

    point_dimension = networkSize*N_EQ1;
    current_state = np.zeros(point_dimension, float)
    temp_state = np.zeros(point_dimension, float)
    k1 = np.zeros(point_dimension, float)
    k2 = np.zeros(point_dimension, float)
    k3 = np.zeros(point_dimension, float)
    k4 = np.zeros(point_dimension, float)
    lastCrossings = np.zeros(point_dimension/2, int)
    numberOfCrossingsToEdit = np.zeros(point_dimension/2, int)
    bm_factor = np.zeros(networkSize, float)
    initial_states_start = 0
    lib.integrate_n_rk4_phasedifferences_electricExcitatoryAndInhibitoryConnections(ct.c_uint(networkSize), initial_states.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     parameterSet.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_chlorine.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_potassium.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_sodium.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_electric.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     X_out.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     P_out.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     numberOfPhaseDifferences.ctypes.data_as(ct.POINTER(ct.c_uint)),
                                                                     current_state.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     temp_state.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     k1.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     k2.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     k3.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     k4.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     lastCrossings.ctypes.data_as(ct.POINTER(ct.c_int)),
                                                                     numberOfCrossingsToEdit.ctypes.data_as(ct.POINTER(ct.c_int)),
                                                                     bm_factor.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     ct.c_double(dt), ct.c_uint(N_integrate), ct.c_uint(stride), ct.c_uint(initial_states_start))
    #print "hasCycles = " + str(hasCycles(P_out, numberOfPhaseDifferences, networkSize));
    P_out = P_out[0:numberOfPhaseDifferences[0]*(networkSize-1)]
    #np.set_printoptions(precision=20)
    #print np.reshape(P_out, (numberOfPhaseDifferences[0], networkSize-1), 'C')
    return np.reshape(X_out, (N_integrate, networkSize), 'C'), np.reshape(P_out, (numberOfPhaseDifferences[0], networkSize-1), 'C')


lib.hasCycles.argtypes=[ ct.POINTER(ct.c_double),ct.POINTER(ct.c_uint), ct.c_uint]
def hasCycles(phase_differences_output,  numberOfPhaseDifferences, networkSize) :
    return  lib.hasCycles(  phase_differences_output.ctypes.data_as(ct.POINTER(ct.c_double)),
                            numberOfPhaseDifferences.ctypes.data_as(ct.POINTER(ct.c_uint)),
                            ct.c_uint(networkSize)  );

#function that performs a sweep across different initial conditions for phase differences and returns the final point of the trajectory
#                                                                   for each of the initial conditions given
lib.sweeptraces_electricExcitatoryAndInhibitoryConnections.argtypes = [ct.c_uint, ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
                                                 ct.POINTER(ct.c_double),
                                                 ct.c_uint,ct.c_double, ct.c_uint, ct.c_uint, ct.c_uint]
def sweeptraces_electricExcitatoryAndInhibitoryConnections(networkSize, initial_states_array, noOfInitialStates,
                coupling_strengths_chlorine, coupling_strengths_potassium,
                coupling_strengths_sodium, coupling_strengths_electric,
                dt, N_integrate, cycles, stride=1,parameterSet=None):
    initial_states_array = np.asarray(initial_states_array) #
    assert initial_states_array.size == N_EQ1*networkSize*noOfInitialStates

    coupling_strengths_chlorine = np.asarray(coupling_strengths_chlorine)
    assert coupling_strengths_chlorine.size == networkSize * networkSize
    coupling_strengths_potassium = np.asarray(coupling_strengths_potassium)
    assert coupling_strengths_potassium.size == networkSize * networkSize
    coupling_strengths_sodium = np.asarray(coupling_strengths_sodium)
    assert coupling_strengths_sodium.size == networkSize * networkSize
    coupling_strengths_electric = np.asarray(coupling_strengths_electric)
    assert coupling_strengths_electric.size == networkSize * networkSize

    trajectory_targets_phasedifferences = np.zeros((networkSize-1)*noOfInitialStates, float)
    print 'cycles : ', cycles

    assert parameterSet.size == NP*networkSize
    lib.sweeptraces_electricExcitatoryAndInhibitoryConnections(ct.c_uint(networkSize), initial_states_array.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     parameterSet.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_chlorine.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_potassium.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_sodium.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     coupling_strengths_electric.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     trajectory_targets_phasedifferences.ctypes.data_as(ct.POINTER(ct.c_double)),
                                                                     ct.c_uint(noOfInitialStates),ct.c_double(dt), ct.c_uint(N_integrate),
                                                                     ct.c_uint(cycles), ct.c_uint(stride))
    return np.reshape(trajectory_targets_phasedifferences, (noOfInitialStates, networkSize-1), 'C')



def parameters_generic(networkSize, parameter_input_array=None):
    #   I_0=0.4, epsilon_0=0.3, x_0=0., k_0=10. , E_Cl=-1.5, E_K=-1.2, E_Na=1.5, m_0=1. ; #### sigma_0=0.
    if parameter_input_array == None :
        parameter_array_onecell = np.array([0.4, 0.3, 0., 10., -1.5, -1.2, 1.5, 1.])
    else :
        parameter_array_onecell = np.array(parameter_input_array)

    parameter_array = parameter_array_onecell
    for i in range(1,networkSize):
        parameter_array = np.concatenate((parameter_array,parameter_array_onecell))
    return parameter_array


def nullcline_x(x, I, m=1., E=0., g=0.):
    return m*(x - x**3) + I + g*(E - x)


def nullcline_y(x, V_0, k):
    return 1./(1.+np.exp(-k*(x-V_0)))


if __name__ == "__main__":
    lib.main.argtypes = []
    lib.main()
