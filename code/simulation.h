#ifndef _SIMULATION_H
#define	_SIMULATION_H

#define IM1 2147483563 
#define IM2 2147483399 
#define AM (1.0/IM1) 
#define IMM1 (IM1-1) 
#define IA1 40014 
#define IA2 40692 
#define IQ1 53668 
#define IQ2 52774 
#define IR1 12211 
#define IR2 3791 
#define NTAB 32
#define NDIV (1+IMM1/NTAB) 
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#define PI (3.141592653589793)
#include "izhi.h"
#include "alphasynapse.h"

using namespace std;
using namespace Cells;
using namespace synapse;

class simulation {
    public:
		simulation();
		~simulation();
		void run();
		void setSeed(long seed);
		void setGin(double gin);
		void setGexc(double gexc);	
		void setH(double h);
		void setTmax(double tmax);
		void setNNeuron(int num);
		void setM(int m);
		void setProbConn(double prob);
		void setProbExc(double probEx);
		void set_percInh(double perc_inh);
		void set_percex1(int perc_ex1);
		void set_percex2(int perc_ex2);
		void set_sn(int sn);
		void set_tStim(double tStim);
		void set_InpAmp(double InpAmp);
		void set_n_stim(int n_stim_aux);
		void rasterdata(char * pchar);
		void LifeTime(char * pchar);
		void set_NameFile(char * NameFile);
		void set_Type_Inh(int type_inh);
		void set_Type_Ex1(int type_ex1);
		void set_Type_Ex2(int type_ex2);
		
		void set_changeNneurons(int changeNneurons);
		void set_numIC(int numIC);
		void set_tStimAfterChange(double tStimAfterChange);
		void set_InpAmpAfterChange(int InpAmpAfterChange);
		
		void set_TimeFreeFalltoTrajectory(double TimeFreeFall);
		void set_dtTrajectory(double dtTraj);
		void set_dtStepsTrajectory(int dtStepsTraj);
		
		void printTime(string name);
		void setNumSim(int NumSim);
		void createNet();
		double LastSpike ();
		//double * avgVolt, * volt, * lfp, *I0;
		double noiseAmp;
    protected:
    	float randNext();
		float a_inh, b_inh, c_inh, d_inh;
		float a_ex1, b_ex1, c_ex1, d_ex1;
		float a_ex2, b_ex2, c_ex2, d_ex2;
    	vector <izhi*> * listn;
    	vector <alphasynapse*> * lists;
    	long idum;
    	int nneuron, m, nMax, changeNneurons, numIC, InpAmpAfterChange;
		double tStimAfterChange, TimeFreeFall, dtTraj;
		int dtStepsTraj;
		double diffStim;
    	double h, tmax, InpAmp;
    	double tStim, gin, gexc;
    	double prob, probEx;
		int * ex12in;
		char * NameFile;
		int perc_inh, perc_ex1, perc_ex2, type_inh, type_ex1, type_ex2, n_stim, n_stim_aux, NNeuron, sn, NumSim;
    //regi[]
};
    
#endif
