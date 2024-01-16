#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N_EQ1	7
#define NP1	4

#define C_m	1.	  // uF/cm^2
#define V_Na	30.	  // mV
#define V_Ca	140.	  // mV
#define V_K	-75.	  // mV
#define V_L	-40.	  // mV
#define g_Na	4.	// mmho/cm^2
#define g_Ca	0.01	// mmho/cm^2
#define g_K	0.3	// mmho/cm^2
#define g_K_Ca	0.03	// mmho/cm^2
#define g_L	0.003	// mmho/cm^2
#define g_h 0.001
#define C_1	127./105.
#define C_2	8265./105.
#define lambda	1./12.5
//#define tau_x	9400	// set 1 ms ??
#define tau_x	235     // set2
//#define A	0.15	// ??
//#define B	-50.	// ??
//#define rho	0.00015	// set 1 ms^-1
#define rho	0.0003  // set2
//#define K_c	0.00425	// set 1 mV^-1
#define K_c	0.0085	// set 2
#define THRESHOLD_SLOPE	100.	// mV^-1
#define THRESHOLD	0.  // mV
#define E_Cl_inh   -80
#define E_K_inh    -75
#define E_Na_Exc -10
/* passed as parameters for each cell
#define E_syn	-80  // mV
#define shift	-57
#define alpha	.01
#define beta	0.001
*/
#define th	-20

#pragma acc routine seq
double alpha_m(const double V)		{ return 0.1*(50.-V)/(exp((50.-V)/10.)-1.); }
#pragma acc routine seq
double beta_m(const double V)		{ return 4.*exp((25.-V)/18.); }
#pragma acc routine seq
double alpha_h(const double V)		{ return 0.07*exp((25.-V)/20.); }
#pragma acc routine seq
double beta_h(const double V)		{ return 1./(1.+exp((55.-V)/10.)); }
#pragma acc routine seq
double alpha_n(const double V)		{ return 0.01*(55.-V)/(exp((55.-V)/10.)-1.); }
#pragma acc routine seq
double beta_n(const double V)		{ return 0.125*exp((45.-V)/80.); }

#pragma acc routine seq
double m_inf(const double V_tilde)	{ return alpha_m(V_tilde)/(alpha_m(V_tilde)+beta_m(V_tilde)); }
#pragma acc routine seq
double h_inf(const double V_tilde)	{ return alpha_h(V_tilde)/(alpha_h(V_tilde)+beta_h(V_tilde)); }
#pragma acc routine seq
double tau_h(const double V_tilde)	{ return 1./(alpha_h(V_tilde)+beta_h(V_tilde)); }
#pragma acc routine seq
double n_inf(const double V_tilde)	{ return alpha_n(V_tilde)/(alpha_n(V_tilde)+beta_n(V_tilde)); }
#pragma acc routine seq
double tau_n(const double V_tilde)	{ return 1./(alpha_n(V_tilde)+beta_n(V_tilde)); }
//double x_inf(const double V)		{ return 1./(1.+exp(A*(B-V))); }
#pragma acc routine seq
//double x_inf(const double V)		{ return ((1/(exp(-0.3*(V+40))+1))); }     //  set 1

double x_inf(const double V)          { return (1/(exp(0.15*(-V-50))+1)); }   // set 2
 

#pragma acc routine seq
double boltzmann(const double V, const double V_0, const double k)
{
	return 1./(1.+exp(-k*(V-V_0)));
}

#pragma acc routine seq
void derivs_one_electricExcitatoryAndInhibitoryConnections(double* y, double* dydt, const double* p)
{
	double V=y[0], h=y[1], n=y[2], x=y[3], Ca=y[4], S=y[5], z=y[6];
	double V_tilde = C_1*V+C_2;// V_Ca=p[0];

    // parameter array p = [E_syn_alpha	-80, shift	-57, alpha	.01, beta	0.001]

	dydt[0] = -g_Na*pow(m_inf(V_tilde), 3)*h*(V-V_Na) - g_Ca*x*(V-V_Na) - (g_K*pow(n, 4)+g_K_Ca*Ca/(0.5+Ca))*(V-V_K) - g_L*(V-V_L) - g_h*(pow((1./(1.+exp(-(V+63.)/7.8))),3)*z*(V - 70.));	// dV/dt
        dydt[1] = lambda*(h_inf(V_tilde)-h)/tau_h(V_tilde);	// dh/dt
        dydt[2] = lambda*(n_inf(V_tilde)-n)/tau_n(V_tilde);	// dn/dt
        dydt[3] = (x_inf(V)-x)/tau_x;				// dx/dt
        dydt[4] = rho * (K_c * x * (V_Ca - V + p[1]) - Ca);	// d[Ca2+]/dt
        dydt[5] = p[2] * (1. - S)/(1.+exp(-10.*(V-th))) - p[3]*S;  //alpha syn
        dydt[6] = 0.5*((1./(1+exp((V-50.)/0.1)))-z)/(7.1+10.4/(1+exp((V+68.)/2.2)));
        //printf("%lf \n", y[5]);
}

/*
This is different from integrate_n_rk4 in that it returns both V and X whereas integrate_n_rk4 returns just V
*/
void integrate_one_rk4(double* y, const double* params, double* output, const double dt, const unsigned N, const unsigned stride)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ1], y2[N_EQ1], k1[N_EQ1], k2[N_EQ1], k3[N_EQ1], k4[N_EQ1];

	dt2 = dt/2.; dt6 = dt/6.;
	for(j=0; j<N_EQ1; j++) output[j] = y[j];

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_one_electricExcitatoryAndInhibitoryConnections(y, k1, params);
			for(k=0; k<N_EQ1; k++) y1[k] = y[k]+k1[k]*dt2;
			derivs_one_electricExcitatoryAndInhibitoryConnections(y1, k2, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y[k]+k2[k]*dt2;
			derivs_one_electricExcitatoryAndInhibitoryConnections(y2, k3, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y[k]+k3[k]*dt;
			derivs_one_electricExcitatoryAndInhibitoryConnections(y2, k4, params);
			for(k=0; k<N_EQ1; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<N_EQ1; j++) output[N_EQ1*i+j] = y[j];
	}
}

/*
Computes the following derivatives for each neuron and also adds the input from other neurons :
     x' = m (x-x^3)-y+I
     y' = eps * (Bfun(x, x_0, k_2) - y)
         n => number of neurons
         y => array with structure x0,y0,x1,y1,x2,y2,...x{n-1},y{n-1} for the above equations
         dydt => derivatives of y i.e., x0',y0',x1',y1',x2',y2'...x{n-1}',y{n-1}'
         p => parameter array storing NP1 parameters for each neuron i.e., p[0] to p[NP1-1] for x0',y0' followed by
                                p[NP] to p[2*NP-1] for x1',y1' and so on
         coupling_strengths_type => array of size n*n storing coupling strengths in the order (0->0),(0->1),(0->2) ... (0->n-1) followed by
                              (1->0),(1->1),(1->2),(1->3)...(1->n-1) and so on

     xi'_coupled = xi' + (E0_chlorine - xi) * ( coupling_strengths_chlorine(0->i)*bm_factor(x0) + coupling_strength_chlorine(1->i)*bm_factor(x1) + ... )
                    + (E0_potassium - xi) * ( coupling_strengths_potassium(0->i)*bm_factor(x0) + coupling_strength_potassium(1->i)*bm_factor(x1) + ... )
                    + (E0_sodium - xi) * ( coupling_strengths_sodium(0->i)*bm_factor(x0) + coupling_strength_sodium(1->i)*bm_factor(x1) + ... )
                    + (x0 - xi) * coupling_strengths_chlorine(0->i) + (x1 - xi) * coupling_strengths_chlorine(1->i) + ... ;
     yi'_coupled = yi'
     bm_factor here uses different k, x0 than in the computation of uncoupled y' ;
     A neuron j will influence neuron i only if it's voltage is above the threshold obtained from the boltzmann function; See bolzmann.py
     TODO : Instead of using three different arrays, should use a linked list for the coupling strengths; For each neuron 'i', the correspoinding linked list will contain nodes each having the presynaptic neuron number, synaptic strenght, and type of channel (or the channel's equilibrium potential); This improves performance from O(n^2) to O(E) where E is no.of edges or synaptic connections
*/
#pragma acc routine seq
void derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(const unsigned n, double* y, double* dydt,
                                const double* p, const double* coupling_strengths_chlorine, const double* coupling_strengths_potassium,
                                const double* coupling_strengths_sodium, const double* coupling_strengths_electric, const double* coupling_alpha_syn,
                                double* bm_factor)
{
	unsigned i, j;
	for(i=0;i <n; i++){
	    bm_factor[i] = boltzmann(y[i*N_EQ1], THRESHOLD, THRESHOLD_SLOPE);
		derivs_one_electricExcitatoryAndInhibitoryConnections(y+i*N_EQ1, dydt+i*N_EQ1, p+i*NP1);
	}
	for(i=0; i<n; i++)
	{
		double bm_factor_sum = 0;
		for(j=0;j<n;j++){
		//E_syn = p[i*NP1+0]
		   dydt[i*N_EQ1] += (p[i*NP1+0] - y[i*N_EQ1])* coupling_alpha_syn[j*n+i] * y[j*N_EQ1+5] +
					(E_Na_Exc -y[i*N_EQ1])*coupling_strengths_sodium[j*n+i] * y[j*N_EQ1+5] +
					(y[j*N_EQ1]-y[i*N_EQ1])*coupling_strengths_electric[j*n+i];
		}
	}
}


/*
    This method finds the rk4 approximation of the curve for a system of 'n' neurons obeying the following equations :
         x' = m (x-x^3)-y+I
         y' = eps * (Bfun(x, x_0, k_2) - y)
     networkSize => number of neurons
     y => array with structure x0,y0,x1,y1,x2,y2,...x{n-1},y{n-1} for the above equations
     dydt => derivatives of y i.e., x0',y0',x1',y1',x2',y2'...x{n-1}',y{n-1}'
     p => parameter array storing NP1 parameters for each neuron i.e., p[0] to p[NP1-1] for x0',y0' followed by
                            p[NP] to p[2*NP-1] for x1',y1' and so on
     coupling_strengths => array of size n*n storing coupling strengths in the order (0->0),(0->1),(0->2) ... (0->n-1) followed by
                          (1->0),(1->1),(1->2),(1->3)...(1->n-1) and so on
     output => array of size n*N in the following order x0,x1,...x{n-1} representing initial point on the curve
                        followed by x0-1,x1-1,...x{n-1}-1 representing the next point and so on
     dt => step size
     N => number of points to compute
     stride => number of strides to make while computing each of the N points
     burstLength => minimum distance between spikes in order to consider them to be from different bursts

     also in the output for each time point 1 to N, check crossings for each neuron;
     for each neuron store all the crossings and also the number of crossings present;
     compute and return phase differences of all cells w.r.t the first cell phase differences

     Whenever there is a spike and its distance from the last spike is > burstLength, it is considered a crossing for the burst
*/
#pragma acc routine seq
void integrate_n_rk4_phasedifferences_electricExcitatoryAndInhibitoryConnections(const unsigned networkSize,
                    const double* initial_state, const double* p,
                    const double* coupling_strengths_chlorine, const double* coupling_strengths_potassium,
                    const double* coupling_strengths_sodium, const double* coupling_strengths_electric,const double* coupling_alpha_syn,
                    double* output, double* phase_differences_output, unsigned* numberOfPhaseDifferences,
                    double* current_state, double* temp_state, double*k1, double*k2, double*k3, double*k4,
                    int* lastCrossings, int* numberOfCrossingsToEdit, double *bm_factor,
                    int* lastSpikes, const int burstLength,
                    const double dt, const unsigned N, const unsigned stride, const unsigned initial_states_start) {

    unsigned i, j, m, point_dimension = networkSize*N_EQ1;
	int threshhold = 0, allCellsHadZeroPhaseBefore = 0, numberOfCrossingsUnedited = 0, phase_differences_output_index=0;

	//const int burstLength = 5000;
    //create another array lastSpikeCross which refers to where the spike crossed threshhold last time, where as lastCrossings would refer to where the burst started like in the case of Fitzhugh_Nagumo

	for(m=0; m<point_dimension; m++){
	    current_state[m]=initial_state[m+initial_states_start];
	    if(m%N_EQ1==0) {
	        output[m/N_EQ1]=current_state[m];
		//printf("%lf\n",current_state[m]);
	        lastCrossings[m/N_EQ1]=N+1;
	        lastSpikes[m/N_EQ1]=N+1;
            numberOfCrossingsToEdit[m/N_EQ1] = 0;
	    }
	}

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, current_state, k1, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				temp_state[m] = current_state[m]+k1[m]*dt/2;
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, temp_state, k2, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				temp_state[m] = current_state[m]+k2[m]*dt/2;
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, temp_state, k3, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				temp_state[m] = current_state[m]+k3[m]*dt;
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, temp_state, k4, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				current_state[m] += (k1[m]+2.*(k2[m]+k3[m])+k4[m])*dt/6;
		}
		for(m=0; m<point_dimension; m++){
            if(m%N_EQ1==0) {
                //Outputs only voltages of each cell.
                output[i*point_dimension/N_EQ1+m/N_EQ1]=current_state[m];
		//printf("%lf\n",current_state[m]);
                // computing the crossings here itself

                if(m==0){ // first cell
                    if(current_state[m]>=threshhold && output[(i-1)*point_dimension/N_EQ1+m/N_EQ1] < threshhold){
                        if(i-lastSpikes[m/N_EQ1] > burstLength){
                            //implies a crossing in the first cell
                            if(!allCellsHadZeroPhaseBefore){
                                int allCellsHadZeroPhaseBefore1 = 1;
                                for(j=1;j<point_dimension/N_EQ1;j++){
                                    if(i<lastCrossings[j]) {
                                        allCellsHadZeroPhaseBefore1=0;
                                    }
                                }
                                allCellsHadZeroPhaseBefore = allCellsHadZeroPhaseBefore1;
                            }
                            if(allCellsHadZeroPhaseBefore){
                                for(j=1;j<point_dimension/N_EQ1;j++){
                                    phase_differences_output[phase_differences_output_index*(point_dimension/N_EQ1 -1) + j-1] = i - lastCrossings[j];
                                    numberOfCrossingsToEdit[j]=numberOfCrossingsToEdit[j]+1;
                                }
                                phase_differences_output_index++;
                            }
                        }
                        lastSpikes[m/N_EQ1]=i;
                    }
                } else{ // cells other than first cell
                    if(current_state[m]>=threshhold && output[(i-1)*point_dimension/N_EQ1+m/N_EQ1] < threshhold){
                        if(i-lastSpikes[m/N_EQ1] > burstLength){
                            if(!allCellsHadZeroPhaseBefore){
                                lastCrossings[m/N_EQ1]=i;
                            }else{
                                j=numberOfCrossingsToEdit[m/N_EQ1];
                                while(j>0){
                                    phase_differences_output[(phase_differences_output_index-j)*(point_dimension/N_EQ1-1)+m/N_EQ1-1]/=i-lastCrossings[m/N_EQ1];
                                    j--;
                                }
                                lastCrossings[m/N_EQ1]=i;
                                numberOfCrossingsToEdit[m/N_EQ1]=0;
                            }
                        }
                        lastSpikes[m/N_EQ1]=i;
                    }
                }

            }
		}
	}

    for(m=0; m<point_dimension; m++){
	    if(m%N_EQ1==0) {
            if(numberOfCrossingsToEdit[m/N_EQ1]>numberOfCrossingsUnedited){
                numberOfCrossingsUnedited = numberOfCrossingsToEdit[m/N_EQ1];
            }
	    }
	}
    //number of the phase differences in the array
    *numberOfPhaseDifferences = phase_differences_output_index - numberOfCrossingsUnedited;
        
}

/*
    Same as integrate_n_rk4_phasedifferences_electricExcitatoryAndInhibitoryConnections
    No need to return output array, so uses a cyclic array for output

    Need to be updated for plant_alpha in this file, currently unusable
*/
#pragma acc routine seq
void integrate_n_rk4_phasedifferences_nooutput_electricExcitatoryAndInhibitoryConnections(const unsigned networkSize,
                    const double* initial_state, const double* p,
                    const double* coupling_strengths_chlorine, const double* coupling_strengths_potassium, 
                    const double* coupling_strengths_sodium, const double* coupling_strengths_electric,const double* coupling_alpha_syn,
                    double* output, double* phase_differences_output, unsigned* numberOfPhaseDifferences,
                    double* current_state, double* temp_state, double*k1, double*k2, double*k3, double*k4,
                    int* lastCrossings, int* numberOfCrossingsToEdit, double *bm_factor,
                    const double dt, const unsigned N, const unsigned stride, const unsigned initial_states_start) {

    unsigned i, j, m, point_dimension = networkSize*N_EQ1;
	int threshhold = 0, allCellsHadZeroPhaseBefore = 0, numberOfCrossingsUnedited = 0, phase_differences_output_index=0;

	for(m=0; m<point_dimension; m++){
	    current_state[m]=initial_state[m+initial_states_start];
	    if(m%2==0) {
	        output[m/2]=current_state[m];
	        lastCrossings[m/2]=N+1;
            numberOfCrossingsToEdit[m/2] = 0;
	    }
	}

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, current_state, k1, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				temp_state[m] = current_state[m]+k1[m]*dt/2;
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, temp_state, k2, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				temp_state[m] = current_state[m]+k2[m]*dt/2;
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, temp_state, k3, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric,coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				temp_state[m] = current_state[m]+k3[m]*dt;
			derivs_n_variablecoupling_electricExcitatoryAndInhibitoryConnections(networkSize, temp_state, k4, p, coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric,coupling_alpha_syn, bm_factor);

			for(m=0; m<point_dimension; m++)
				current_state[m] += (k1[m]+2.*(k2[m]+k3[m])+k4[m])*dt/6;
		}
		for(m=0; m<point_dimension; m++){
            if(m%2==0) {
                //Outputs only voltages of each cell.
                output[(i%2)*point_dimension/2+m/2]=current_state[m];

                // computing the crossings here itself
                if(m==0){ // first cell
                    if(current_state[m]>=threshhold && output[((i-1)%2)*point_dimension/2+m/2] < threshhold){
                        //implies a crossing in the first cell
                        if(!allCellsHadZeroPhaseBefore){
                            int allCellsHadZeroPhaseBefore1 = 1;
                            for(j=1;j<point_dimension/2;j++){
                                if(i<lastCrossings[j]) {
                                    allCellsHadZeroPhaseBefore1=0;
                                }
                            }
                            allCellsHadZeroPhaseBefore = allCellsHadZeroPhaseBefore1;
                        }
                        if(allCellsHadZeroPhaseBefore){
                            for(j=1;j<point_dimension/2;j++){
                                phase_differences_output[phase_differences_output_index*(point_dimension/2 -1) + j-1] = i - lastCrossings[j];
                                numberOfCrossingsToEdit[j]=numberOfCrossingsToEdit[j]+1;
                            }
                            phase_differences_output_index++;
                        }
                    }
                } else{ // cells other than first cell
                    if(current_state[m]>=threshhold && output[((i-1)%2)*point_dimension/2+m/2] < threshhold){
                        if(!allCellsHadZeroPhaseBefore){
                            lastCrossings[m/2]=i;
                        }else{
                            j=numberOfCrossingsToEdit[m/2];
                            while(j>0){
                                phase_differences_output[(phase_differences_output_index-j)*(point_dimension/2-1)+m/2-1]/=i-lastCrossings[m/2];
                                j--;
                            }
                            lastCrossings[m/2]=i;
                            numberOfCrossingsToEdit[m/2]=0;
                        }
                    }
                }

            }
		}
	}

    for(m=0; m<point_dimension; m++){
	    if(m%2==0) {
            if(numberOfCrossingsToEdit[m/2]>numberOfCrossingsUnedited){
                numberOfCrossingsUnedited = numberOfCrossingsToEdit[m/2];
            }
	    }
	}
    //number of the phase differences in the array
    *numberOfPhaseDifferences = phase_differences_output_index - numberOfCrossingsUnedited;

}


/*
Takes phase difference output array and checks if there is any cycle;
If a cycle is found, computes the average of all the points in the cycle and
appends to phase differences array and increments the count as well by one.
*/
int hasCycles(double* phase_differences_output, unsigned* numberOfPhaseDifferences, unsigned networkSize) {
    int hasCycles=0, index=0, m;
    double cyclePhaseDifference[networkSize-1];

    for( m=0; m<networkSize-1; m++) {
        cyclePhaseDifference[m]=0;
    }

    for (index=*numberOfPhaseDifferences-2; index>=0 ; index--) {
        for(m=0; m<networkSize-1; m++) {

            double currentPhaseDifference = phase_differences_output[index*(networkSize-1)+m];
            cyclePhaseDifference[m] += currentPhaseDifference;
            hasCycles=1;

            // Matching
            if( currentPhaseDifference
                  != phase_differences_output[(*numberOfPhaseDifferences-1)*(networkSize-1)+m]) {
                  hasCycles=0;
                  break;
            }
        }
        if(hasCycles==1){
            //Compute average
            for(m=0; m<networkSize-1; m++){
                cyclePhaseDifference[m] /= *numberOfPhaseDifferences-index-1;
                phase_differences_output[*numberOfPhaseDifferences*(networkSize-1)+m] = cyclePhaseDifference[m];
            }
            *numberOfPhaseDifferences++;

            return hasCycles;
        }
    }

    return hasCycles;

}

/*
for an array of initial conditions, do everything in integrate_n_rk4_phasedifferences and return just the last point in phase differences trajectory for each initial condition
*/
void sweeptraces_electricExcitatoryAndInhibitoryConnections(const long networkSize, const double* initial_states_array, const double* p,
                    const double* coupling_strengths_chlorine, const double* coupling_strengths_potassium,
                    const double* coupling_strengths_sodium, const double* coupling_strengths_electric,const double* coupling_alpha_syn,
                    double* trajectory_targets_phasedifferences,
                    const long noOfInitialStates, const double dt, const long N, int cycles, const long stride) {

    int l,m, point_dimension = networkSize*N_EQ1;

    //double output[N*networkSize]; not working when 'output' copied privately since size exceeds 2.3GB for all threads together; As of pgcc v15.5
    //double output[N*networkSize*noOfInitialStates];

    //experimenting with stateless computations; not using full output array since we are only using phase_difference array; also making it private since
    // size won't be greater than 2.3GB anymore
    double output[2*networkSize];

    //double phase_differences_output[(networkSize-1)*(N/2+1)*noOfInitialStates]; declaring it here gives errors with pgcc, hence creating array in python.
    double phase_differences_output[(networkSize-1)*(cycles)];

    unsigned numberOfPhaseDifferences=0;
    double current_state[point_dimension], temp_state[point_dimension], k1[point_dimension], k2[point_dimension], k3[point_dimension], k4[point_dimension];
    int lastCrossings[point_dimension/2], numberOfCrossingsToEdit[point_dimension/2];
    double bm_factor[networkSize];
    struct timeval start_time, stop_time, elapsed_time;
    gettimeofday(&start_time,NULL);

    /*
	In order to avoid memory issues, the number of threads, regulated here by num_gangs(16) vector_length(128), could also be dynamic, based on the size of N*networkSize*noOfInitialStates .
	Also, for very large resolutions, say 150x150, rather than clogging memory and thereby slowing down everything, it is instead better to run in multiple batches of say, 50x50 -
	something that makes a reasonable use of GPU memory.
	*/
	/*#pragma omp parallel for \
        private(phase_differences_output,numberOfPhaseDifferences),\
        private(current_state,temp_state,k1,k2,k3,k4,lastCrossings,numberOfCrossingsToEdit),\
        private(bm_factor)*/
        //copy(output)
	/*
	private(output) ; not working for >2.3GB
	num_gangs(16) vector_length(128) ; Can set the maximum no.of parallel threads with this; Could adopt dynamically based on sizes but not worth the wait.
	                                                    Rather, if memory requirement's too great to be fulfilled, split the input array and run in batches
    private(phase_differences_output,numberOfPhaseDifferences),\
    */
    #pragma acc parallel loop independent ,\
        pcopyin(initial_states_array[0:point_dimension*noOfInitialStates],p[0:NP1*networkSize],coupling_strengths_chlorine[0:networkSize*networkSize],coupling_strengths_potassium[0:networkSize*networkSize],coupling_strengths_sodium[0:networkSize*networkSize],coupling_strengths_electric[0:networkSize*networkSize],coupling_alpha_syn[0:networkSize*networkSize]), \
        copy(trajectory_targets_phasedifferences[0:(networkSize-1)*noOfInitialStates]), \
        private(current_state,temp_state,k1,k2,k3,k4,lastCrossings,numberOfCrossingsToEdit),\
        private(bm_factor),\
        private(output,numberOfPhaseDifferences),\
        private(phase_differences_output[0:(networkSize-1)*(cycles)])
	for (l=0;l<noOfInitialStates;l++){
        /*
        integrate_n_rk4_phasedifferences_nooutput(networkSize,initial_states_array, p, coupling_strengths,
                    output+l*N*networkSize, phase_differences_output, &numberOfPhaseDifferences,
                    current_state,temp_state,k1,k2,k3,k4,
                    lastCrossings,numberOfCrossingsToEdit,bm_factor,
                    dt, N, stride, l*point_dimension);
        */
        integrate_n_rk4_phasedifferences_nooutput_electricExcitatoryAndInhibitoryConnections(networkSize,initial_states_array, p,
                    coupling_strengths_chlorine, coupling_strengths_potassium, coupling_strengths_sodium, coupling_strengths_electric, coupling_alpha_syn,
                    output, phase_differences_output, &numberOfPhaseDifferences,
                    current_state,temp_state,k1,k2,k3,k4,
                    lastCrossings,numberOfCrossingsToEdit,bm_factor,
                    dt, N, stride, l*point_dimension);

        for(m=0;m<networkSize-1;m++){
	    if(numberOfPhaseDifferences>0)
            	trajectory_targets_phasedifferences[l*(networkSize-1)+m]=
                    phase_differences_output[(numberOfPhaseDifferences-1)*(networkSize-1)+m];
        }
    }

    /*
	printf("\n");
	for(l=0; l<noOfInitialStates;l++){
		for(m=0;m<networkSize-1;m++)printf("%f ",trajectory_targets_phasedifferences[l*(networkSize-1)+m]);
		printf("\n");
	}
	fflush(stdout);
    */

    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

}


/*
int main(){
	// test case for 6 cells
	double *initial_states_array;	
	double p[] =    {0.4,0.3,0.,10.,-1.5,-1.2,1.5,1.,
	                0.4,0.3,0.,10.,-1.5,-1.2,1.5,1.,
	                0.4,0.3,0.,10.,-1.5,-1.2,1.5,1.,
	                0.4,0.3,0.,10.,-1.5,-1.2,1.5,1.,
	                0.4,0.3,0.,10.,-1.5,-1.2,1.5,1.,
	                0.4,0.3,0.,10.,-1.5,-1.2,1.5,1};
	double coupling_strengths_chlorine[] = {0.,0.001,0.001,0.001,0.001,0.,
	                                        0.,0.001,0.001,0.001,0.001,0.,
                                                0.,0.001,0.001,0.001,0.001,0.,
                                                0.,0.001,0.001,0.001,0.001,0.,
                                                0.,0.001,0.001,0.001,0.001,0.,
                                                0.,0.001,0.001,0.001,0.001,0.};
	double *trajectory_targets_phasedifferences, *phase_differences_output;
	void *ptr ;
	int i,j,n, networkSize = 6, N = 10065, cycles = 100;

	printf("Enter size of grid\n");.
	scanf("%d",&n);
	trajectory_targets_phasedifferences = (double*)malloc((networkSize-1)*n*n*sizeof(double));
	ptr = malloc(networkSize*N_EQ1*n*n*sizeof(double));
	if(!ptr)
	{
		printf("malloc returns NULL\n");
		return 1;
	}
	initial_states_array = (double*)ptr;
	double initial_state[] = {-0.026541128437339426, 0.07017957749966305, -0.026541128437339426, 0.07017957749966305, -0.026541128437339426, 0.07017957749966305, -0.026541128437339426, 0.07017957749966305, -0.026541128437339426, 0.07017957749966305,
	                            -0.026541128437339426, 0.07017957749966305};
	for(i=0;i<n*n;i++){
		for(j=0;j<networkSize*N_EQ1;j++){
			initial_states_array[i*networkSize*N_EQ1+j]=initial_state[j];
		}	
	}	
	sweeptraces_electricExcitatoryAndInhibitoryConnections(networkSize, initial_states_array, p,
                    coupling_strengths_chlorine, coupling_strengths_chlorine, coupling_strengths_chlorine,
                    coupling_strengths_chlorine,
                    trajectory_targets_phasedifferences, n*n,
					(double)3/1000, N, cycles, 30);
}


*/
