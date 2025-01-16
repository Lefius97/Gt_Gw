#include "correlator.h"
#include <math.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>

int main(int argc, char *argv[]) {

	if (argc<2) {
		std::cout << "Usage:" << std::endl 
			<< "Program INPUT_FILE OUTPUT_FILE_G(T) OUTPUT_FILE_G(W) DT(YOUR_SIMALLESTDATA=DT*TIME_UNIT) MODE G(T)CUT" << std::endl
			<< "Note: You need to multiply the result of the modulus by V/kBT" << std::endl
			<< "      MODE==0:G(t) used stress in xy,xz,yz" << std::endl
			<< "      MODE==1:G(t) used stress in xx,yy,zz,xy,xz,yz; and your input file should put xx,yy,zz in first three columns" << std::endl
			<< "      G(t)cut could make your G(t)==0 after t>=t_cut; t_cut can be selected according to the first column in the G(t) file and this t_cut defaults to the time maximum" << std::endl;
		return 1;
	}

	Correlator c;
	c.setsize(32, 16, 2);    

	std::ifstream fin;
	double valab,valac,valbc,valaa,valbb,valcc;
	
	fin.open(argv[1]);
	std::ofstream fout;
	fout.open(argv[2]);
	std::ofstream fout2;
	fout2.open(argv[3]);
	double dt=atof(argv[4]);
	int mode=atoi(argv[5]);
	int Gtcut=-1;
	if(argv[6]!='\0'){
		Gtcut=atoi(argv[6]);
		if(Gtcut<0){
			std::cout << "ERROR G(T)CUT!" << std::endl;
			return 1;
		}
	}
	if(!fin.is_open() || fin.bad()) {
		std::cout << "ERROR INPUT FILE!" << std::endl;
		return 1;
	}
	if(mode!=0 && mode!=1) {
		std::cout << "ERROR MODE SET!" << std::endl;
		return 1;
	}
	
	int tt;
	if(mode==0){
		c.initialize_mode0();
		tt=0;
		while(!fin.eof()) {
			fin >> valab >> valac >> valbc;
			c.add_mode0 (valab,valac,valbc);
			tt++;
		}
		fin.close();
		c.evaluate_mode0();
		fout << "#Time " << "G_t " << std::endl;
		for (unsigned int i=0;i<c.npcorr;++i){
			c.Gt[i] = (c.f1_mode0[i] + c.f2_mode0[i] + c.f3_mode0[i]) / 3.0;
			fout << c.t[i]*dt << " " << c.Gt[i] << std::endl;
		}
		fout.close();
	}else if(mode==1){
		c.initialize_mode1();
		tt=0;
		while(!fin.eof()) {
			fin >> valaa >> valbb >> valcc >> valab >> valac >> valbc;
			//c.add_mode1 (valaa-valbb,valaa-valcc,valbb-valcc,valab,valac,valbc);
			c.add_mode1 (valaa,valaa,valbb,valab,valac,valbc);
			tt++;
		}
		fin.close();
		c.evaluate_mode1();
		fout << "#Time " << "G_t " << std::endl;
		for (unsigned int i=0;i<c.npcorr;++i){
			c.Gt[i] = (c.f1_mode1[i] + c.f2_mode1[i] + c.f3_mode1[i]) / 30.0 + (c.f4_mode1[i] + c.f5_mode1[i] + c.f6_mode1[i]) / 5.0;
			fout << i << " " <<c.t[i]*dt << " " << c.Gt[i] << std::endl;
		}
		fout.close();
	}
	
	c.initialize_Gw();
	if(Gtcut==-1) Gtcut=c.npcorr;
	c.calculate_Gw(dt,tt,Gtcut);
	fout2 << "#Frequency " << "Gw_storage " << "Gw_loss " << "Gw " << std::endl;
	for (unsigned int i=c.npcorr-1;i>0;--i){
		fout2 << 1.0/(c.wt[i]*dt) << " " << c.Gw_storage[i] << " " << c.Gw_loss[i] << " " << c.Gw[i] << std::endl;
	}
	fout2.close();

	return 0;
}
