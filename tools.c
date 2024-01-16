
/*
 * Three-cell Leech heart interneuron model is implimented 
 * as in the draft  "Bifurcations of bursting polyrhythms
 * in 3-cell central pattern generators" of Jeremy Wojcik,
 * Robert Clewley, and Andrey Shilnikov, prepared for
 * Journal of Neuroscience 2012.
 */


#include <math.h>
#include <stdio.h>


#define NGLS	9
#define NP	12
#define Esyn	-0.0625  // Sajia: same
#define EmNa	-0.0305  // Sajia: 0.0305
#define EhNa	-0.0325  // Sajia: 0.0333
#define V_K2	-0.018   // Sajia: same


// for a given array 'x' of voltages and a 'threshhold', returns in an array 'ti', the indices j of x where crossings take place i.e.,
//          where x[j] < threshhold and x[j+1] > threshhold
//  TODO : should be >=
void crossings(const double* x, const unsigned N, double* ti, const double threshold)
{
	unsigned i, j;
	i = 0;
	for(j=1; j<N; j++)
	{
		if(x[j] > threshold && x[j-1] < threshold)
		{
			ti[i++] = (double)(j-1)+(threshold-x[j-1])/(x[j]-x[j-1]);	
		}
	}
}


void derivs_one(const double* y, double* dxdt, const double* p) {
	dxdt[0] = ( -p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2]) - p[4]*y[2]*y[2]*(y[0]-p[1]) - p[3]*(y[0]-p[0]) - p[9])/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];
}


void derivs_two(const double* y, double* dxdt, const double* p, const double* kij) {
	double coupling = kij[0]*(Esyn-y[0])/(1.+exp(-1000.*(y[3]-p[11])));
	dxdt[0] = (-p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2])-p[4]*y[2]*y[2]*(y[0]-p[1])-p[3]*(y[0]-p[0])-p[9]+coupling)/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];													
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];

	coupling = kij[1]*(Esyn-y[3])/(1.+exp(-1000.*(y[0]-p[11])));
	dxdt[3] = ( -p[5]*pow(1./(1.+exp(-150.*(y[3]-EmNa))), 3.)*y[4]*(y[3]-p[2]) - p[4]*y[5]*y[5]*(y[3]-p[1]) - p[3]*(y[3]-p[0]) - p[9] + coupling )/p[6];
	dxdt[4] = (1./(1.+exp(500.*(y[3]-EhNa)))-y[4])/p[7];
	dxdt[5] = (1./(1.+exp(-83.*(y[3]-V_K2+p[10])))-y[5])/p[8];
}


void derivs_three(const double* y, double* dxdt, const double* p, const double* kij) {
	double coupling = kij[0]*(Esyn-y[0])/(1.+exp(-1000.*(y[3]-p[11]))) + kij[1]*(Esyn-y[0])/(1.+exp(-1000.*(y[6]-p[11])));
	dxdt[0] = (-p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2])-p[4]*y[2]*y[2]*(y[0]-p[1])-p[3]*(y[0]-p[0])-p[9]+coupling)/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];													
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];

	coupling = kij[2]*(Esyn-y[3])/(1.+exp(-1000.*(y[0]-p[11]))) + kij[3]*(Esyn-y[3])/(1.+exp(-1000.*(y[6]-p[11])));
	dxdt[3] = ( -p[5]*pow(1./(1.+exp(-150.*(y[3]-EmNa))), 3.)*y[4]*(y[3]-p[2]) - p[4]*y[5]*y[5]*(y[3]-p[1]) - p[3]*(y[3]-p[0]) - p[9] + coupling )/p[6];
	dxdt[4] = (1./(1.+exp(500.*(y[3]-EhNa)))-y[4])/p[7];
	dxdt[5] = (1./(1.+exp(-83.*(y[3]-V_K2+p[10])))-y[5])/p[8];

	coupling = kij[4]*(Esyn-y[6])/(1.+exp(-1000.*(y[0]-p[11]))) + kij[5]*(Esyn-y[6])/(1.+exp(-1000.*(y[3]-p[11])));
	dxdt[6] = (-p[5]*pow(1./(1.+exp(-150.*(y[6]-EmNa))), 3.)*y[7]*(y[6]-p[2])-p[4]*y[8]*y[8]*(y[6]-p[1])-p[3]*(y[6]-p[0])-p[9]+coupling)/p[6];	     
	dxdt[7] = (1./(1.+exp(500.*(y[6]-EhNa)))-y[7])/p[7];												     
	dxdt[8] = (1./(1.+exp(-83.*(y[6]-V_K2+p[10])))-y[8])/p[8];											     
}

void leech_three_rk4(double* y, const double dt, const unsigned N, const unsigned stride, const double* p, const double* kij, double* output) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[NGLS], y2[NGLS], k1[NGLS], k2[NGLS], k3[NGLS], k4[NGLS], P[NP];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<NP; j++) P[j] = p[j+1];
	for(j=0; j<3; j++) output[j] = y[3*j]; 	
	for(i=1; i<N; i++) {
		for(j=0; j<stride; j++) {
			derivs_three(y, k1, P, kij);
			for(k=0; k<NGLS; k++) y1[k] = y[k]+k1[k]*dt2; 			
			derivs_three(y1, k2, P, kij);
			for(k=0; k<NGLS; k++) y2[k] = y[k]+k2[k]*dt2; 			
			derivs_three(y2, k3, P, kij);
			for(k=0; k<NGLS; k++) y2[k] = y[k]+k3[k]*dt; 			
			derivs_three(y2, k4, P, kij);
			for(k=0; k<NGLS; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<3; j++) output[3*i+j] = y[3*j]; 					
	}
}

void leech_three_signal_rk4(double* Kt, double* y, const double dt, const unsigned N, const unsigned stride, const double* p, double* kij, double* output) {
	unsigned i, j, k, nKt=0;
	double dt2, dt6;
	double y1[NGLS], y2[NGLS], k1[NGLS], k2[NGLS], k3[NGLS], k4[NGLS], P[NP];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<NP; j++) P[j] = p[j+1];
	for(j=0; j<3; j++) output[j] = y[3*j]; 	
	for(i=1; i<N; i++) {
		for(j=0; j<stride; j++) {
			kij[2] = Kt[nKt]; kij[3] = Kt[nKt++];
			derivs_three(y, k1, P, kij);
			for(k=0; k<NGLS; k++) y1[k] = y[k]+k1[k]*dt2; 			

			kij[2] = Kt[nKt]; kij[3] = Kt[nKt++];
			derivs_three(y1, k2, P, kij);
			for(k=0; k<NGLS; k++) y2[k] = y[k]+k2[k]*dt2; 			
			derivs_three(y2, k3, P, kij);
			for(k=0; k<NGLS; k++) y2[k] = y[k]+k3[k]*dt; 			

			kij[2] = Kt[nKt]; kij[3] = Kt[nKt];
			derivs_three(y2, k4, P, kij);

			for(k=0; k<NGLS; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<3; j++) output[3*i+j] = y[3*j]; 					
	}
}

void leech_three_rk4_full_output(double* y, const double dt, const unsigned N, const unsigned stride, const double* p, const double* kij, double* output) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[9], y2[9], k1[9], k2[9], k3[9], k4[9], P[NP];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<NP; j++) P[j] = p[j+1];
	for(j=0; j<9; j++) output[j] = y[j];
	for(i=1; i<N; i++) {
		for(j=0; j<stride; j++) {
			derivs_three(y, k1, P, kij);
			for(k=0; k<9; k++) y1[k] = y[k]+k1[k]*dt2; 			
			derivs_three(y1, k2, P, kij);
			for(k=0; k<9; k++) y2[k] = y[k]+k2[k]*dt2; 			
			derivs_three(y2, k3, P, kij);
			for(k=0; k<9; k++) y2[k] = y[k]+k3[k]*dt; 			
			derivs_three(y2, k4, P, kij);
			for(k=0; k<9; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<9; j++) output[9*i+j] = y[j];
	}
}


void leech_two_rk4(double* y, const double dt, const unsigned N, const double* p, const double* kij, double* output) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[6], y2[6], k1[6], k2[6], k3[6], k4[6], P[NP];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<NP; j++) P[j] = p[j+1];
	for(j=0; j<3; j++) output[j] = y[3*j]; 	
	for(i=0; i<N; i++) {
		derivs_two(y, k1, P, kij);
		for(k=0; k<6; k++) y1[k] = y[k]+k1[k]*dt2; 			
		derivs_two(y1, k2, P, kij);
		for(k=0; k<6; k++) y2[k] = y[k]+k2[k]*dt2; 			
		derivs_two(y2, k3, P, kij);
		for(k=0; k<6; k++) y2[k] = y[k]+k3[k]*dt; 			
		derivs_two(y2, k4, P, kij);
		for(k=0; k<6; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
	}
	for(j=0; j<6; j++) output[j] = y[j]; 					
}


void leech_one_rk4(double* y, const double dt, const unsigned N, const unsigned stride, const double* p, double* output) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[3], y2[3], k1[3], k2[3], k3[3], k4[3], P[NP];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<NP; j++) P[j] = p[j+1];
	for(j=0; j<3; j++) output[j] = y[j];
	for(i=1; i<N; i++) {
		for(j=0; j<stride; j++) {
			derivs_one(y, k1, P);
			for(k=0; k<3; k++) y1[k] = y[k]+k1[k]*dt2; 			
			derivs_one(y1, k2, P);
			for(k=0; k<3; k++) y2[k] = y[k]+k2[k]*dt2; 			
			derivs_one(y2, k3, P);
			for(k=0; k<3; k++) y2[k] = y[k]+k3[k]*dt; 			
			derivs_one(y2, k4, P);
			for(k=0; k<3; k++) y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<3; j++) output[3*i+j] = y[j];
	}
}






