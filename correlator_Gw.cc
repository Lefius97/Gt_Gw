#include "correlator.h"
#include <math.h>

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

void Correlator::initialize_Gw() {

	for (unsigned int i=0;i<length;++i) {
		wt[i] = 0;
		Gw[i] = 0;
		Gw_storage[i] = 0;
		Gw_loss[i] = 0;
	}

}

void Correlator::calculate_Gw(const double dtime, const int time, const int Gtcut) {
	
	double X1=0.0,X2=0.0,costk_1=0.0,sintk_1=0.0,coswdtk=0.0,sinwdtk=0.0,dtk,w2;
	double dt=dtime;
	int tt=time;
	double* wtt=new double [tt];
	double* Gw_1=new double [tt];
	double* Gw_2=new double [tt];
	unsigned int im=0;
	int judge=Gtcut;

	for (unsigned int i=1;i<tt;++i) {
		wtt[i]=1/(i*dt);
		w2=wtt[i]*wtt[i];
		for (unsigned int j=1;j<npcorr;++j){
			
			if(i==1){
				if(j>int(0.5*npcorr) && Gt[j]<0){
					judge=j;
				}
				if(j>=judge){
					Gt[j]=0;
				}
				//std::cout << i<<" "<<j << " " << Gt[j] << " " << judge << std::endl;
			}
			
			costk_1 = cos(wtt[i]*t[j-1]*dt);
			sintk_1 = sin(wtt[i]*t[j-1]*dt);
			dtk = (t[j]-t[j-1])*dt;
			coswdtk = cos(wtt[i]*dtk);
			sinwdtk = sin(wtt[i]*dtk);
			
			X1 = Gt[j]*( (coswdtk-1)/(w2*dtk) + sinwdtk/wtt[i] ) - Gt[j-1]*( (coswdtk-1)/(w2*dtk) );
			X2 = Gt[j]*( coswdtk/wtt[i] - sinwdtk/(w2*dtk) ) - Gt[j-1]*( 1.0/wtt[i] - sinwdtk/(w2*dtk) );
			
			Gw_1[i] += sintk_1*X1 - costk_1*X2;
			Gw_2[i] += sintk_1*X2 + costk_1*X1;
		}
		Gw_1[i] = Gw_1[i]*wtt[i];
		Gw_2[i] = Gw_2[i]*wtt[i];
		//std::cout << wtt[i] << " " << Gw_1[i] << " " << Gw_2[i] << std::endl;
	}

	for(int i=0;i<p;++i){
		wt[im]=i;
		Gw_storage[im]=Gw_1[im];
		Gw_loss[im]=Gw_2[im];
		Gw[im] = sqrt( Gw_storage[im]*Gw_storage[im] + Gw_loss[im]*Gw_loss[im] );
		++im;
	}
	int ss=im;
	for (int k=1;k<kmax;++k) {
		for (int i=dmin;i<p;++i) {
			wt[im]=i * pow((double)m, k);
			++im;
		}
	}

	for (int k=1;k<kmax;++k) {
		for (int i=dmin;i<p;++i) {
			int lgth=0;
			if(ss+1>=npcorr) break;
			for(int s=wt[ss-1];s<=wt[ss+1];++s){
				Gw_storage[ss]+=Gw_1[s];
				Gw_loss[ss]+=Gw_2[s];
				lgth++;
			}
			Gw_storage[ss]=Gw_storage[ss]/lgth;
			Gw_loss[ss]=Gw_loss[ss]/lgth;
			Gw[ss] = sqrt( Gw_storage[ss]*Gw_storage[ss] + Gw_loss[ss]*Gw_loss[ss] );
			ss++;
		}
	}

}

