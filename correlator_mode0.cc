#include "correlator.h"
#include <math.h>

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

Correlator::Correlator(const unsigned int numcorrin, const unsigned int pin, const unsigned int min) {
	setsize(numcorrin, pin, min);
}

Correlator::~Correlator() {

	if (numcorrelators ==0) return;
//Mode 0
	delete [] shift;
	delete [] correlation;
	delete [] ncorrelation;
	delete [] accumulator;
	delete [] naccumulator;
	delete [] insertindex;
	
	delete [] shift2;
	delete [] correlation2;
	delete [] ncorrelation2;
	delete [] accumulator2;
	delete [] naccumulator2;
	delete [] insertindex2;
	
	delete [] shift3;
	delete [] correlation3;
	delete [] ncorrelation3;
	delete [] accumulator3;
	delete [] naccumulator3;
	delete [] insertindex3;
	
	delete [] f1_mode0;
	delete [] f2_mode0;
	delete [] f3_mode0;
	
//Mode 1
	delete [] shift_aa;
	delete [] correlation_aa;
	delete [] ncorrelation_aa;
	delete [] accumulator_aa;
	delete [] naccumulator_aa;
	delete [] insertindex_aa;
	
	delete [] shift_bb;
	delete [] correlation_bb;
	delete [] ncorrelation_bb;
	delete [] accumulator_bb;
	delete [] naccumulator_bb;
	delete [] insertindex_bb;
	
	delete [] shift_cc;
	delete [] correlation_cc;
	delete [] ncorrelation_cc;
	delete [] accumulator_cc;
	delete [] naccumulator_cc;
	delete [] insertindex_cc;
	
	delete [] shift_ab;
	delete [] correlation_ab;
	delete [] ncorrelation_ab;
	delete [] accumulator_ab;
	delete [] naccumulator_ab;
	delete [] insertindex_ab;
	
	delete [] shift_ac;
	delete [] correlation_ac;
	delete [] ncorrelation_ac;
	delete [] accumulator_ac;
	delete [] naccumulator_ac;
	delete [] insertindex_ac;
	
	delete [] shift_bc;
	delete [] correlation_bc;
	delete [] ncorrelation_bc;
	delete [] accumulator_bc;
	delete [] naccumulator_bc;
	delete [] insertindex_bc;

	delete [] f1_mode1;
	delete [] f2_mode1;
	delete [] f3_mode1;
	delete [] f4_mode1;
	delete [] f5_mode1;
	delete [] f6_mode1;
}

void Correlator::setsize (const unsigned int numcorrin, const unsigned int pin, const unsigned int min) {
	numcorrelators = numcorrin;
	p = pin;
	m = min;
	dmin = p/m;

	length = numcorrelators*p; 
	
//Mode 0

	shift = new double* [numcorrelators];
	correlation = new double* [numcorrelators];
	ncorrelation = new unsigned long int* [numcorrelators];
	accumulator = new double[numcorrelators];
	naccumulator = new unsigned int [numcorrelators];
	insertindex = new unsigned int [numcorrelators];
	
	shift2 = new double* [numcorrelators];
	correlation2 = new double* [numcorrelators];
	ncorrelation2 = new unsigned long int* [numcorrelators];
	accumulator2 = new double[numcorrelators];
	naccumulator2 = new unsigned int [numcorrelators];
	insertindex2 = new unsigned int [numcorrelators];
	
	shift3 = new double* [numcorrelators];
	correlation3 = new double* [numcorrelators];
	ncorrelation3 = new unsigned long int* [numcorrelators];
	accumulator3 = new double[numcorrelators];
	naccumulator3 = new unsigned int [numcorrelators];
	insertindex3 = new unsigned int [numcorrelators];
	
//Mode 1

	shift_aa = new double* [numcorrelators];
	correlation_aa = new double* [numcorrelators];
	ncorrelation_aa = new unsigned long int* [numcorrelators];
	accumulator_aa = new double[numcorrelators];
	naccumulator_aa = new unsigned int [numcorrelators];
	insertindex_aa = new unsigned int [numcorrelators];
	
	shift_bb = new double* [numcorrelators];
	correlation_bb = new double* [numcorrelators];
	ncorrelation_bb = new unsigned long int* [numcorrelators];
	accumulator_bb = new double[numcorrelators];
	naccumulator_bb = new unsigned int [numcorrelators];
	insertindex_bb = new unsigned int [numcorrelators];
	
	shift_cc = new double* [numcorrelators];
	correlation_cc = new double* [numcorrelators];
	ncorrelation_cc = new unsigned long int* [numcorrelators];
	accumulator_cc = new double[numcorrelators];
	naccumulator_cc = new unsigned int [numcorrelators];
	insertindex_cc = new unsigned int [numcorrelators];
	
	shift_ab = new double* [numcorrelators];
	correlation_ab = new double* [numcorrelators];
	ncorrelation_ab = new unsigned long int* [numcorrelators];
	accumulator_ab = new double[numcorrelators];
	naccumulator_ab = new unsigned int [numcorrelators];
	insertindex_ab = new unsigned int [numcorrelators];
	
	shift_ac = new double* [numcorrelators];
	correlation_ac = new double* [numcorrelators];
	ncorrelation_ac = new unsigned long int* [numcorrelators];
	accumulator_ac = new double[numcorrelators];
	naccumulator_ac = new unsigned int [numcorrelators];
	insertindex_ac = new unsigned int [numcorrelators];
	
	shift_bc = new double* [numcorrelators];
	correlation_bc = new double* [numcorrelators];
	ncorrelation_bc = new unsigned long int* [numcorrelators];
	accumulator_bc = new double[numcorrelators];
	naccumulator_bc = new unsigned int [numcorrelators];
	insertindex_bc = new unsigned int [numcorrelators];

	for (unsigned int j=0; j<numcorrelators;++j) {
		//Mode 0
		shift[j] = new double [p];
		shift2[j] = new double [p];
		shift3[j] = new double [p];

		// It can be optimized: Apart from correlator 0, correlation and ncorrelation arrays only use p/2 values
		correlation[j] = new double [p];
		ncorrelation[j] = new unsigned long int [p];
		
		correlation2[j] = new double [p];
		ncorrelation2[j] = new unsigned long int [p];
		
		correlation3[j] = new double [p];
		ncorrelation3[j] = new unsigned long int [p];
		
		//Mode 1
		shift_aa[j] = new double [p];
		shift_bb[j] = new double [p];
		shift_cc[j] = new double [p];
		shift_ab[j] = new double [p];
		shift_ac[j] = new double [p];
		shift_bc[j] = new double [p];

		// It can be optimized: Apart from correlator 0, correlation and ncorrelation arrays only use p/2 values
		correlation_aa[j] = new double [p];
		ncorrelation_aa[j] = new unsigned long int [p];
		
		correlation_bb[j] = new double [p];
		ncorrelation_bb[j] = new unsigned long int [p];
		
		correlation_cc[j] = new double [p];
		ncorrelation_cc[j] = new unsigned long int [p];
		
		correlation_ab[j] = new double [p];
		ncorrelation_ab[j] = new unsigned long int [p];
		
		correlation_ac[j] = new double [p];
		ncorrelation_ac[j] = new unsigned long int [p];
		
		correlation_bc[j] = new double [p];
		ncorrelation_bc[j] = new unsigned long int [p];
	}

	t = new double[length];
	f1_mode0 = new double[length];
	f2_mode0 = new double[length];
	f3_mode0 = new double[length];
	
	f1_mode1 = new double[length];
	f2_mode1 = new double[length];
	f3_mode1 = new double[length];
	f4_mode1 = new double[length];
	f5_mode1 = new double[length];
	f6_mode1 = new double[length];
	
	Gt = new double[length];
	wt = new double[length];
	Gw = new double[length];
	Gw_storage = new double[length];
	Gw_loss = new double[length];
}

void Correlator::initialize_mode0() {

	for (unsigned int j=0;j<numcorrelators;++j) {
		for (unsigned int i=0;i<p;++i) {
			shift[j][i] = -2E10;
			correlation[j][i] = 0;
			ncorrelation[j][i] = 0;
			shift2[j][i] = -2E10;
			correlation2[j][i] = 0;
			ncorrelation2[j][i] = 0;
			shift3[j][i] = -2E10;
			correlation3[j][i] = 0;
			ncorrelation3[j][i] = 0;
		}
		accumulator[j] = 0.0;
		naccumulator[j] = 0;
		insertindex[j] = 0;
		accumulator2[j] = 0.0;
		naccumulator2[j] = 0;
		insertindex2[j] = 0;
		accumulator3[j] = 0.0;
		naccumulator3[j] = 0;
		insertindex3[j] = 0;
	}

	for (unsigned int i=0;i<length;++i) {
		t[i] = 0;
		f1_mode0[i] = 0;
		f2_mode0[i] = 0;
		f3_mode0[i] = 0;
	}

	npcorr =0;
	kmax=0;
	accval=0;
	accval2=0;
	accval3=0;
}

void Correlator::add_mode0(const double w, const double w2, const double w3, const unsigned int k) {

	// If we exceed the correlator side, the value is discarded
	if (k == numcorrelators) return;
	if (k > kmax) kmax=k;

	// Insert new value in shift array
	shift[k][insertindex[k]] = w;
	shift2[k][insertindex2[k]] = w2;
	shift3[k][insertindex3[k]] = w3;

	// Add to average value
	if (k==0){
		accval += w;
		accval2 += w2;
		accval3 += w3;
	}

	// Add to accumulator and, if needed, add to next correlator
	accumulator[k] += w;
	++naccumulator[k];
	accumulator2[k] += w2;
	++naccumulator2[k];
	accumulator3[k] += w3;
	++naccumulator3[k];
	if (naccumulator[k]==m && naccumulator2[k]==m && naccumulator3[k]==m) {
		add_mode0(accumulator[k]/m, accumulator2[k]/m, accumulator3[k]/m, k+1);
		accumulator[k]=0;
		naccumulator[k]=0;
		accumulator2[k]=0;
		naccumulator2[k]=0;
		accumulator3[k]=0;
		naccumulator3[k]=0;
	}

	// Calculate correlation function
	unsigned int ind1=insertindex[k];
	unsigned int ind12=insertindex2[k];
	unsigned int ind13=insertindex3[k];
	if (k==0) { // First correlator is different
		int ind2=ind1;
		int ind22=ind12;
		int ind23=ind13;
		for (unsigned int j=0;j<p;++j) {
			if (shift[k][ind2] > -1e10) {
				correlation[k][j]+= shift[k][ind1]*shift[k][ind2];
				++ncorrelation[k][j];
			}
			--ind2;
			if (ind2<0) ind2+=p;
			
			if (shift2[k][ind22] > -1e10) {
				correlation2[k][j]+= shift2[k][ind12]*shift2[k][ind22];
				++ncorrelation2[k][j];
			}
			--ind22;
			if (ind22<0) ind22+=p;	

			if (shift3[k][ind23] > -1e10) {
				correlation3[k][j]+= shift3[k][ind13]*shift3[k][ind23];
				++ncorrelation3[k][j];
			}
			--ind23;
			if (ind23<0) ind23+=p;				
		}
	}
	else {
		int ind2=ind1-dmin;
		int ind22=ind12-dmin;
		int ind23=ind13-dmin;
		for (unsigned int j=dmin;j<p;++j) {
			if (ind2<0) ind2+=p;				
			if (shift[k][ind2] > -1e10) {
				correlation[k][j]+= shift[k][ind1]*shift[k][ind2];
				++ncorrelation[k][j];				
			}
			--ind2;
			
			if (ind22<0) ind22+=p;				
			if (shift2[k][ind22] > -1e10) {
				correlation2[k][j]+= shift2[k][ind12]*shift2[k][ind22];
				++ncorrelation2[k][j];				
			}
			--ind22;
			
			if (ind23<0) ind23+=p;				
			if (shift3[k][ind23] > -1e10) {
				correlation3[k][j]+= shift3[k][ind13]*shift3[k][ind23];
				++ncorrelation3[k][j];				
			}
			--ind23;
		}
	}

	++insertindex[k];
	if (insertindex[k]==p) insertindex[k]=0;
	++insertindex2[k];
	if (insertindex2[k]==p) insertindex2[k]=0;
	++insertindex3[k];
	if (insertindex3[k]==p) insertindex3[k]=0;
}

void Correlator::evaluate_mode0(const bool norm) {
	unsigned int im=0;

	double aux=0;
	double aux2=0;
	double aux3=0;
	if (norm){
		aux = (accval/ncorrelation[0][0])*(accval/ncorrelation[0][0]);
		aux2 = (accval2/ncorrelation2[0][0])*(accval2/ncorrelation2[0][0]);
		aux3 = (accval3/ncorrelation3[0][0])*(accval3/ncorrelation3[0][0]);
	}

	// First correlator
	for (unsigned int i=0;i<p;++i) {
		if (ncorrelation[0][i]>0 && ncorrelation2[0][i]>0 && ncorrelation3[0][i]>0) {
			t[im] = i;
			f1_mode0[im] = correlation[0][i]/ncorrelation[0][i] - aux;
			f2_mode0[im] = correlation2[0][i]/ncorrelation2[0][i] - aux2;
			f3_mode0[im] = correlation3[0][i]/ncorrelation3[0][i] - aux3;
			++im;
		}
	}

	// Subsequent correlators
	for (int k=1;k<kmax;++k) {
		for (int i=dmin;i<p;++i) {
			if (ncorrelation[k][i]>0 && ncorrelation2[k][i]>0 && ncorrelation3[k][i]>0) {
				t[im] = i * pow((double)m, k);
				f1_mode0[im] = correlation[k][i] / ncorrelation[k][i] - aux;
				f2_mode0[im] = correlation2[k][i] / ncorrelation2[k][i] - aux2;
				f3_mode0[im] = correlation3[k][i] / ncorrelation3[k][i] - aux3;
				++im;
			}
		}
	}
	

	npcorr = im;
}

