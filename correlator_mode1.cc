#include "correlator.h"
#include <math.h>

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

void Correlator::initialize_mode1() {

	for (unsigned int j=0;j<numcorrelators;++j) {
		for (unsigned int i=0;i<p;++i) {
			shift_aa[j][i] = -2E10;
			correlation_aa[j][i] = 0;
			ncorrelation_aa[j][i] = 0;
			
			shift_bb[j][i] = -2E10;
			correlation_bb[j][i] = 0;
			ncorrelation_bb[j][i] = 0;
			
			shift_cc[j][i] = -2E10;
			correlation_cc[j][i] = 0;
			ncorrelation_cc[j][i] = 0;
			
			shift_ab[j][i] = -2E10;
			correlation_ab[j][i] = 0;
			ncorrelation_ab[j][i] = 0;
			
			shift_ac[j][i] = -2E10;
			correlation_ac[j][i] = 0;
			ncorrelation_ac[j][i] = 0;
			
			shift_bc[j][i] = -2E10;
			correlation_bc[j][i] = 0;
			ncorrelation_bc[j][i] = 0;
		}
		accumulator_aa[j] = 0.0;
		naccumulator_aa[j] = 0;
		insertindex_aa[j] = 0;
		
		accumulator_bb[j] = 0.0;
		naccumulator_bb[j] = 0;
		insertindex_bb[j] = 0;
		
		accumulator_cc[j] = 0.0;
		naccumulator_cc[j] = 0;
		insertindex_cc[j] = 0;
		
		accumulator_ab[j] = 0.0;
		naccumulator_ab[j] = 0;
		insertindex_ab[j] = 0;
		
		accumulator_ac[j] = 0.0;
		naccumulator_ac[j] = 0;
		insertindex_ac[j] = 0;
		
		accumulator_bc[j] = 0.0;
		naccumulator_bc[j] = 0;
		insertindex_bc[j] = 0;
	}

	for (unsigned int i=0;i<length;++i) {
		t[i] = 0;
		f1_mode1[i] = 0;
		f2_mode1[i] = 0;
		f3_mode1[i] = 0;
		f4_mode1[i] = 0;
		f5_mode1[i] = 0;
		f6_mode1[i] = 0;
	}

	npcorr =0;
	kmax=0;
	accval_aa=0;
	accval_bb=0;
	accval_cc=0;
	accval_ab=0;
	accval_ac=0;
	accval_bc=0;
}

void Correlator::add_mode1(const double w1, const double w2, const double w3, const double w4, const double w5, const double w6, const unsigned int k) {

	// If we exceed the correlator side, the value is discarded
	if (k == numcorrelators) return;
	if (k > kmax) kmax=k;

	// Insert new value in shift array
	shift_aa[k][insertindex_aa[k]] = w1;
	shift_bb[k][insertindex_bb[k]] = w2;
	shift_cc[k][insertindex_cc[k]] = w3;
	shift_ab[k][insertindex_ab[k]] = w4;
	shift_ac[k][insertindex_ac[k]] = w5;
	shift_bc[k][insertindex_bc[k]] = w6;

	// Add to average value
	if (k==0){
		accval_aa += w1;
		accval_bb += w2;
		accval_cc += w3;
		accval_ab += w4;
		accval_ac += w5;
		accval_bc += w6;
	}

	// Add to accumulator and, if needed, add to next correlator
	accumulator_aa[k] += w1;
	++naccumulator_aa[k];
	accumulator_bb[k] += w2;
	++naccumulator_bb[k];
	accumulator_cc[k] += w3;
	++naccumulator_cc[k];
	accumulator_ab[k] += w4;
	++naccumulator_ab[k];
	accumulator_ac[k] += w5;
	++naccumulator_ac[k];
	accumulator_bc[k] += w6;
	++naccumulator_bc[k];
	if (naccumulator_aa[k]==m && naccumulator_bb[k]==m && naccumulator_cc[k]==m && naccumulator_ab[k]==m && naccumulator_ac[k]==m && naccumulator_bc[k]==m) {
		add_mode1(accumulator_aa[k]/m, accumulator_bb[k]/m, accumulator_cc[k]/m, accumulator_ab[k]/m, accumulator_ac[k]/m, accumulator_bc[k]/m, k+1);
		accumulator_aa[k]=0;
		naccumulator_aa[k]=0;
		accumulator_bb[k]=0;
		naccumulator_bb[k]=0;
		accumulator_cc[k]=0;
		naccumulator_cc[k]=0;
		accumulator_ab[k]=0;
		naccumulator_ab[k]=0;
		accumulator_ac[k]=0;
		naccumulator_ac[k]=0;
		accumulator_bc[k]=0;
		naccumulator_bc[k]=0;
	}

	// Calculate correlation function
	unsigned int ind1_aa=insertindex_aa[k];
	unsigned int ind1_bb=insertindex_bb[k];
	unsigned int ind1_cc=insertindex_cc[k];
	unsigned int ind1_ab=insertindex_ab[k];
	unsigned int ind1_ac=insertindex_ac[k];
	unsigned int ind1_bc=insertindex_bc[k];
	if (k==0) { // First correlator is different
		int ind2_aa=ind1_aa;
		int ind2_bb=ind1_bb;
		int ind2_cc=ind1_cc;
		int ind2_ab=ind1_ab;
		int ind2_ac=ind1_ac;
		int ind2_bc=ind1_bc;
		for (unsigned int j=0;j<p;++j) {
			if (shift_aa[k][ind2_aa] > -1e10) {
				correlation_aa[k][j]+= shift_aa[k][ind1_aa]*shift_aa[k][ind2_aa];
				++ncorrelation_aa[k][j];
			}
			--ind2_aa;
			if (ind2_aa<0) ind2_aa+=p;
			
			if (shift_bb[k][ind2_bb] > -1e10) {
				correlation_bb[k][j]+= shift_bb[k][ind1_bb]*shift_bb[k][ind2_bb];
				++ncorrelation_bb[k][j];
			}
			--ind2_bb;
			if (ind2_bb<0) ind2_bb+=p;
				
			if (shift_cc[k][ind2_cc] > -1e10) {
				correlation_cc[k][j]+= shift_cc[k][ind1_cc]*shift_cc[k][ind2_cc];
				++ncorrelation_cc[k][j];
			}
			--ind2_cc;
			if (ind2_cc<0) ind2_cc+=p;
			
			if (shift_ab[k][ind2_ab] > -1e10) {
				correlation_ab[k][j]+= shift_ab[k][ind1_ab]*shift_ab[k][ind2_ab];
				++ncorrelation_ab[k][j];
			}
			--ind2_ab;
			if (ind2_ab<0) ind2_ab+=p;
			
			if (shift_ac[k][ind2_ac] > -1e10) {
				correlation_ac[k][j]+= shift_ac[k][ind1_ac]*shift_ac[k][ind2_ac];
				++ncorrelation_ac[k][j];
			}
			--ind2_ac;
			if (ind2_ac<0) ind2_ac+=p;
				
			if (shift_bc[k][ind2_bc] > -1e10) {
				correlation_bc[k][j]+= shift_bc[k][ind1_bc]*shift_bc[k][ind2_bc];
				++ncorrelation_bc[k][j];
			}
			--ind2_bc;
			if (ind2_bc<0) ind2_bc+=p;				
		}
	}
	else {
		int ind2_aa=ind1_aa-dmin;
		int ind2_bb=ind1_bb-dmin;
		int ind2_cc=ind1_cc-dmin;
		int ind2_ab=ind1_ab-dmin;
		int ind2_ac=ind1_ac-dmin;
		int ind2_bc=ind1_bc-dmin;
		for (unsigned int j=dmin;j<p;++j) {
			if (ind2_aa<0) ind2_aa+=p;				
			if (shift_aa[k][ind2_aa] > -1e10) {
				correlation_aa[k][j]+= shift_aa[k][ind1_aa]*shift_aa[k][ind2_aa];
				++ncorrelation_aa[k][j];				
			}
			--ind2_aa;
			
			if (ind2_bb<0) ind2_bb+=p;				
			if (shift_bb[k][ind2_bb] > -1e10) {
				correlation_bb[k][j]+= shift_bb[k][ind1_bb]*shift_bb[k][ind2_bb];
				++ncorrelation_bb[k][j];				
			}
			--ind2_bb;
			
			if (ind2_cc<0) ind2_cc+=p;				
			if (shift_cc[k][ind2_cc] > -1e10) {
				correlation_cc[k][j]+= shift_cc[k][ind1_cc]*shift_cc[k][ind2_cc];
				++ncorrelation_cc[k][j];				
			}
			--ind2_cc;
			
			if (ind2_ab<0) ind2_ab+=p;				
			if (shift_ab[k][ind2_ab] > -1e10) {
				correlation_ab[k][j]+= shift_ab[k][ind1_ab]*shift_ab[k][ind2_ab];
				++ncorrelation_ab[k][j];				
			}
			--ind2_ab;
			
			if (ind2_ac<0) ind2_ac+=p;				
			if (shift_ac[k][ind2_ac] > -1e10) {
				correlation_ac[k][j]+= shift_ac[k][ind1_ac]*shift_ac[k][ind2_ac];
				++ncorrelation_ac[k][j];				
			}
			--ind2_ac;
			
			if (ind2_bc<0) ind2_bc+=p;				
			if (shift_bc[k][ind2_bc] > -1e10) {
				correlation_bc[k][j]+= shift_bc[k][ind1_bc]*shift_bc[k][ind2_bc];
				++ncorrelation_bc[k][j];				
			}
			--ind2_bc;
		}
	}

	++insertindex_aa[k];
	if (insertindex_aa[k]==p) insertindex_aa[k]=0;
	++insertindex_bb[k];
	if (insertindex_bb[k]==p) insertindex_bb[k]=0;
	++insertindex_cc[k];
	if (insertindex_cc[k]==p) insertindex_cc[k]=0;
	++insertindex_ab[k];
	if (insertindex_ab[k]==p) insertindex_ab[k]=0;
	++insertindex_ac[k];
	if (insertindex_ac[k]==p) insertindex_ac[k]=0;
	++insertindex_bc[k];
	if (insertindex_bc[k]==p) insertindex_bc[k]=0;
}

void Correlator::evaluate_mode1(const bool norm) {
	unsigned int im=0;

	double aux_aa=0;
	double aux_bb=0;
	double aux_cc=0;
	double aux_ab=0;
	double aux_ac=0;
	double aux_bc=0;
	if (norm){
		aux_aa = (accval_aa/ncorrelation_aa[0][0])*(accval_aa/ncorrelation_aa[0][0]);
		aux_bb = (accval_bb/ncorrelation_bb[0][0])*(accval_bb/ncorrelation_bb[0][0]);
		aux_cc = (accval_cc/ncorrelation_cc[0][0])*(accval_cc/ncorrelation_cc[0][0]);
		aux_ab = (accval_ab/ncorrelation_ab[0][0])*(accval_ab/ncorrelation_ab[0][0]);
		aux_ac = (accval_ac/ncorrelation_ac[0][0])*(accval_ac/ncorrelation_ac[0][0]);
		aux_bc = (accval_bc/ncorrelation_bc[0][0])*(accval_bc/ncorrelation_bc[0][0]);
	}

	// First correlator
	for (unsigned int i=0;i<p;++i) {
		if (ncorrelation_aa[0][i]>0 && ncorrelation_bb[0][i]>0 && ncorrelation_cc[0][i]>0 && ncorrelation_ab[0][i]>0 && ncorrelation_ac[0][i]>0 && ncorrelation_bc[0][i]>0) {
			t[im] = i;
			f1_mode1[im] = correlation_aa[0][i]/ncorrelation_aa[0][i] - aux_aa;
			f2_mode1[im] = correlation_bb[0][i]/ncorrelation_bb[0][i] - aux_bb;
			f3_mode1[im] = correlation_cc[0][i]/ncorrelation_cc[0][i] - aux_cc;
			f4_mode1[im] = correlation_ab[0][i]/ncorrelation_ab[0][i] - aux_ab;
			f5_mode1[im] = correlation_ac[0][i]/ncorrelation_ac[0][i] - aux_ac;
			f6_mode1[im] = correlation_bc[0][i]/ncorrelation_bc[0][i] - aux_bc;
			++im;
		}
	}

	// Subsequent correlators
	for (int k=1;k<kmax;++k) {
		for (int i=dmin;i<p;++i) {
			if (ncorrelation_aa[k][i]>0 && ncorrelation_bb[k][i]>0 && ncorrelation_cc[k][i]>0 && ncorrelation_ab[k][i]>0 && ncorrelation_ac[k][i]>0 && ncorrelation_bc[k][i]>0) {
				t[im] = i * pow((double)m, k);
				f1_mode1[im] = correlation_aa[k][i] / ncorrelation_aa[k][i] - aux_aa;
				f2_mode1[im] = correlation_bb[k][i] / ncorrelation_bb[k][i] - aux_bb;
				f3_mode1[im] = correlation_cc[k][i] / ncorrelation_cc[k][i] - aux_cc;
				f4_mode1[im] = correlation_ab[k][i] / ncorrelation_ab[k][i] - aux_ab;
				f5_mode1[im] = correlation_ac[k][i] / ncorrelation_ac[k][i] - aux_ac;
				f6_mode1[im] = correlation_bc[k][i] / ncorrelation_bc[k][i] - aux_bc;
				++im;
			}
		}
	}
	

	npcorr = im;
}

