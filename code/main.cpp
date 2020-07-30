#include <sstream>
#include <iostream>
#include <string>
#include "simulation.h"



using namespace std;

int main(int argc, char * argv[]) {

	FILE * pfile;
	simulation * sim;
	ostringstream nameFile;	
	double  h = 0.01;	//integration precision
	double prob = 0.01;	//connection probability
	double probEx = 0.9;	//rewiring prob. for ex neurons 
		
	double perc_inh=20;	// % of inhib neurons in the netw
		
	//types of neurons
	//int type_inh=2;  // 1 = FS ; 2 = LTS; 3 = RS
	//int type_ex1=1;  // 1 = RS ; 2 = IB ; 3 = CH
	//int type_ex2=3;  // 1 = RS ; 2 = IB ; 3 = CH
		
	long seed = -184503872;	
			
	sim = new simulation();
	sim->setH(h);
	sim->setSeed(seed+atof(argv[16])*1000000);
	sim->setProbConn(prob);
	sim->setProbExc(probEx);
	sim->set_percInh(perc_inh);
	sim->setM(atof(argv[1]));
	sim->set_percex1(atof(argv[2]));
	sim->set_percex2(atof(argv[3]));
	sim->set_tStim(atof(argv[4]));
	sim->set_InpAmp(atof(argv[5]));
	sim->set_n_stim(atof(argv[6]));
	sim->setGin(atof(argv[7]));
	sim->setGexc(atof(argv[8]));
	sim->set_sn(atof(argv[9]));
	sim->set_Type_Inh(atof(argv[11]));
	sim->set_Type_Ex1(atof(argv[12]));
	sim->set_Type_Ex2(atof(argv[13]));
	sim->setNNeuron(atof(argv[14]));
	sim->setTmax(atof(argv[15]));
	sim->setNumSim(atof(argv[16]));
	sim->set_numIC(atof(argv[17]));
	sim->set_changeNneurons(atof(argv[18]));
	sim->set_tStimAfterChange(atof(argv[19]));
	sim->set_InpAmpAfterChange(atof(argv[20]));
	sim->set_TimeFreeFalltoTrajectory(atof(argv[21]));
	sim->set_dtTrajectory(atof(argv[22]));
	sim->set_dtStepsTrajectory(atof(argv[23]));
	sim->noiseAmp = atof(argv[24]);
	
	sim->createNet();	
	sim->run();
	
	sim->LifeTime(argv[10]);
	sim->rasterdata(argv[10]); 
	
	return 0;
}
