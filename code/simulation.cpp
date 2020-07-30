#include "simulation.h"

simulation::simulation () {
	listn = new vector <izhi*>;
	lists = new vector <alphasynapse*>;		
}

simulation::~simulation () {
	izhi * ptrn = NULL;
	alphasynapse * ptrs = NULL;
	for(int i = 0; i < listn->size(); i++){
		ptrn = listn->at(i);
		delete ptrn;
		ptrn = NULL;
	}
	listn->clear();
	delete listn;
	for(int i = 0; i < lists->size(); i++){
		ptrs = lists->at(i);
		delete ptrs;
		ptrs = NULL;		
	}
	lists->clear();
	delete lists;
}

void simulation::run () {
	
	for (int i = 0; i < listn->size(); i++) {	
		(listn->at(i))->updateNeNi();  
	}
	FILE * pfile;

	ostringstream nameFile;  

	int nMax;
	double t = 0.0;	
	double hn;

	int stim, stim1;
	
	int StimMtr[listn->size()];	
	int StimMtr1[listn->size()];
	
	//the stimulation region
	int nreg=2;
	int regi[6]= {0, n_stim-1}; 
	int si=0;
	for (int i = 0; i <= listn->size(); i++){
	  StimMtr1[i]=0;
	if (i >= regi[si] && i<= regi[si+1] && si<nreg){
	  StimMtr[i]=1;
	  if (i == regi[si+1]){
		si=si+2;
	  }
	}else{
	  StimMtr[i]=0;	
	}	
	}
	
	if (numIC>0){
	//for (int i=0; i<=nneuron*(numIC); i++ ) randNext();	
	int ss=0;
	int ii=0;

	while (ss< changeNneurons){
	  ii++;
	  if (ii==listn->size()) ii=0;
	  if (StimMtr1[ii]==0 && randNext()<= ((double) changeNneurons/(double) nneuron)){
	    StimMtr1[ii]=1;
	    ss++;
	    }	
	}	  
	}
				
	nMax = (int) ((double)tmax/h) + 1;	
	stim = 1;
	stim1 = 0;
	
	
	for (int k = 0; k < nMax; k++) {
        t += h;
	
	if (t > tStim+TimeFreeFall+dtStepsTraj*dtTraj+tStimAfterChange){
		stim = 0; 
		stim1= 0;		
	}else if(t > tStim && 
	         t <= tStim+TimeFreeFall+dtStepsTraj*dtTraj){	  
		stim = 0; 
		stim1= 0;  	
        }else if (t > tStim+TimeFreeFall+dtStepsTraj*dtTraj
	       && t < tStim+TimeFreeFall+dtStepsTraj*dtTraj+tStimAfterChange){
		stim = 0; 
		stim1= 1;		
        }
        
		double vmean = 0, umean=0;
	        for (int i = 0; i < listn->size(); i++) {	  

				int f = 0;
				diffStim = (StimMtr[i] + InpAmpAfterChange*stim1*StimMtr1[i]); 

				
			if (t<tStim){
				(listn->at(i))->evaluate(InpAmp,t);} // (diffStim, t);}
			else{
				(listn->at(i))->evaluate(0, t);}
	
	        }

    }
    

    if ( 1 ){ 		//output of the parameters 
		cout <<endl;
		cout <<NumSim<<" % REPETITION NR"<<endl;
		cout <<endl;
		cout << "%****************"<<endl;
		cout <<sn<<" % SIMULATION NR"<<endl;
		cout <<numIC<<" % IC"<<endl;
		cout << "%****************"<<endl;
		cout <<endl;
		cout <<h<<" % integr precision"<<endl;
		cout <<tmax << "  % tmax" <<endl;
		cout <<nneuron << "  % nneuron" <<endl;
		cout <<prob  << "  % prob" <<endl;
		cout <<probEx << "  % Probexc" <<endl;
		cout <<perc_inh << "  % perc_inh" <<endl;
		cout <<type_inh << "  % type_inh" <<endl;
		cout <<type_ex1 << "  % type_ex1" <<endl;
		cout <<type_ex2 << "  % type_ex2" <<endl;
		cout << "%================"<<endl;
		cout <<m << "  % m" <<endl;
		cout <<perc_ex1 << "  % perc_ex1" <<endl;
		cout <<perc_ex2 << "  % perc_ex2" <<endl;
		cout <<n_stim << "  % n_stim" <<endl;
		cout << "%------------"<<endl;
		cout <<tStim << "  % tStim" <<endl;
		cout <<InpAmp  << "  % InpAmp" <<endl;
		cout << "%------------"<<endl;
		cout <<TimeFreeFall <<"% Time of Free Fall after initial stimulation " <<endl;
		cout << dtTraj <<"% dT along trajectory after free fall " <<endl ;
		cout << dtStepsTraj<<"% # of steps along trajectory after free fall " <<endl;
		cout << dtStepsTraj*dtTraj<<"% => Time along the trajectory after the free fall " <<endl;
		cout << "%------------"<<endl;
		cout <<tStimAfterChange<< "  % tStimAfterChange" <<endl;
		cout <<InpAmpAfterChange  << "  % InpAmpAfterChange" <<endl;
		cout <<changeNneurons  << "  % changeNneurons" <<endl;
		cout << "%------------"<<endl;
		cout <<gin<< "  % gin" <<endl;
		cout <<gexc<< "  % gex" <<endl;
		cout <<endl;
		cout << "%=========================================="<<endl;
		cout << "%=========================================="<<endl;
    }
	
}

double simulation::LastSpike () {
        FILE * pfile;
        ostringstream nameFile;  	
	izhi * nrn;
	double aux = 0.0, last = 0.0;
	for (int i = 0; i < listn->size(); i++){
		nrn = listn->at(i);
		if ((nrn->getevents())->empty()) continue;
		aux = (nrn->getevents())->back();
		if (aux > last) last = aux;
	}
	return last;
}


void simulation::setSeed (long seed) {
	this->idum = seed;
}

float simulation::randNext() {
	int j; 
	long k; 
	static long idum2=123456789; 
	static long iy=0; 
	static long iv[NTAB]; 
	float temp;
	if (idum <= 0) { 
		if (-(idum) < 1) idum=1; 
		else idum = -(idum); 
		idum2=(idum); 
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1; 
			idum=IA1*(idum-k*IQ1)-k*IR1; 
			if (idum < 0) idum += IM1; 
			if (j < NTAB) iv[j] = idum;
		} 
		iy=iv[0];
	} 
	k=(idum)/IQ1; 
	idum=IA1*(idum-k*IQ1)-k*IR1; 
	if (idum < 0) idum += IM1; 
	k=idum2/IQ2; 
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2; 
	j=iy/NDIV; iy=iv[j]-idum2; iv[j] = idum; 
	if (iy < 1) iy += IMM1; 
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp;
}  


void simulation::set_Type_Inh(int type_inh){
	this->type_inh = type_inh;
}

void simulation::set_Type_Ex1(int type_ex1){
	this->type_ex1 = type_ex1;
}

void simulation::set_Type_Ex2(int type_ex2){
	this->type_ex2 = type_ex2;
}

void simulation::setGin(double gin){
	this->gin = gin;
}

void simulation::setGexc(double gexc){
	this->gexc = gexc;	
}

void simulation::setH(double h){
	this->h = h;	
}

void simulation::setTmax(double tmax){
	this->tmax = tmax;	
}

void simulation::setNNeuron(int num){
	this->nneuron = num;	
}

void simulation::setM(int m){
	this->m = m;	
}

void simulation::setProbConn(double prob){
	this->prob = prob;	
}

void simulation::setProbExc(double probEx){
	this->probEx = probEx;	
}

void simulation::set_percex1(int perc_ex1){
	this->perc_ex1 = perc_ex1;	
}

void simulation::set_percex2(int perc_ex2){
	this->perc_ex2 = perc_ex2;	
}

void simulation::set_tStim(double tStim){
	this->tStim = tStim;	
}

void simulation::set_changeNneurons(int changeNneurons){
	this->changeNneurons = changeNneurons;	
}

void simulation::set_numIC(int numIC){
	this->numIC = numIC;	
}

void simulation::set_tStimAfterChange(double tStimAfterChange){
	this->tStimAfterChange = tStimAfterChange;	
}

void simulation::set_InpAmp(double InpAmp){
	this->InpAmp = InpAmp;	
}

void simulation::set_InpAmpAfterChange(int InpAmpAfterChange){
	this->InpAmpAfterChange = InpAmpAfterChange;	
}

void simulation::set_TimeFreeFalltoTrajectory(double TimeFreeFall){
	this->TimeFreeFall = TimeFreeFall;
}

void simulation::set_dtTrajectory(double dtTraj){
	this->dtTraj = dtTraj;  
}

void simulation::set_dtStepsTrajectory(int dtStepsTraj){
	this->dtStepsTraj = dtStepsTraj;  
}


void simulation::set_n_stim(int n_stim_aux){
	if (n_stim_aux==0){
	this->n_stim = 1024;}
	else if (n_stim_aux==1){
	this->n_stim = 512;}
	else if (n_stim_aux==2){
	this->n_stim = 256;	  
	}else if (n_stim_aux==3){
	this->n_stim = 128;}	
	else if (n_stim_aux==4){
	this->n_stim = 64;}
}


void simulation::setNumSim(int NumSim){
	this->NumSim = NumSim;	
}

void simulation::set_sn(int sn){
	this->sn = sn;	
}

void simulation::set_percInh(double perc_inh){
	this->perc_inh = perc_inh;	
}

void simulation::createNet () {
    
    FILE * pfile;
    ostringstream nameFile;
      
    //type of inhibitory neurons 
    if (type_inh == 1){		//FS: fast spiking
	a_inh= 0.1; b_inh= 0.2; c_inh= -65; d_inh= 2;
    }
    else if (type_inh == 2){	//LTS: low-treshold spiking
	a_inh= 0.02; b_inh= 0.25; c_inh= -65; d_inh= 2;  
    }
    else if (type_inh == 3){	//RS: regular spiking
	a_inh= 0.02; b_inh= 0.2; c_inh= -65; d_inh= 8;
    }
    
    //type I of excitatory neurons  
    if (type_ex1==1){	 //RS: regular spiking
      	a_ex1= 0.02; b_ex1= 0.2; c_ex1= -65; d_ex1= 8;
    }
    else if (type_ex1==2){	//IB: intrinsiclly bursting
        a_ex1= 0.02; b_ex1= 0.2; c_ex1= -55; d_ex1= 4;
    }
    else if (type_ex1==3){	//CH: chattering
        a_ex1= 0.02; b_ex1= 0.2; c_ex1= -50; d_ex1= 2;
    }
    
    //type II of excitatory neurons 
    if (type_ex2==1){		 //RS: regular spiking
        a_ex2= 0.02; b_ex2= 0.2; c_ex2= -65; d_ex2= 8;
    }
    else if (type_ex2==2){	//IB: intrinsiclly bursting
        a_ex2= 0.02; b_ex2= 0.2; c_ex2= -55; d_ex2= 4;      
    }
    else if (type_ex2==3){	//CH: chattering    
        a_ex2= 0.02; b_ex2= 0.2; c_ex2= -50; d_ex2= 2;        
    }

    ex12in = new int[nneuron];
  
    long idum = -1089532789; 
    alphasynapse * syn;
    int i, j, count, sizem = nneuron, num, num1, nmod = 1, mod1, mod2, n1, n2, nn;
    double sort;
    int ** listsyn = new int *[1000000];
	for (i=0; i<1000000; i++)
		listsyn[i] = new int[3];		
    
    izhi * cell;
    double rnd;
    for (i = 0; i < nneuron; i++) {
    	cell = new izhi;
        cell->seth(h);
        listn->push_back(cell);    

		num =  (int) (randNext()*100);
		num1 =  (int) (randNext()*100);
        if (num < perc_inh){
        	rnd = randNext();
			cell->setpar(a_inh,b_inh,c_inh,d_inh);
        	cell->setInhibitory();
			ex12in[i]=1;		
        }else{
	  		if (num1 < perc_ex1) {
	        	rnd = randNext();
	       		cell->setpar(a_ex1,b_ex1,c_ex1,d_ex1);
	       		cell->setExcitatory();
				ex12in[i]=0;
	  		}else{
		        rnd = randNext();
	       		cell->setpar(a_ex2,b_ex2,c_ex2,d_ex2);
	       		cell->setExcitatory();
				ex12in[i]=-1;
	  		}
        }
      	cell->set_noiseAmplitude(noiseAmp, numIC*1000);
    }
    
    nameFile.str("");
	nameFile << "ex12in_M" << m<< "_In_"<<type_inh<<"_Ex2_"<< type_ex2<<"_Sim_"<< sn <<"_NumSim_" << NumSim <<  ".dat";
	pfile = fopen((nameFile.str()).c_str(),"w");
	for (int j = 0; j <nneuron; j++){
		fprintf(pfile,"%d\n",ex12in[j]);
	}
	fclose(pfile); 	
    
    count = 0;
	for (i = 0; i < nneuron; i++) {  
		for (j = 0; j < nneuron; j++) {
			sort = randNext();			
			if (j!=i & sort < prob) {
				listsyn[count][0] = i;
				listsyn[count][1] = j;	
				listsyn[count][2] = 1;		
				count++;
			}
		}		
	}	
	listsyn[count][0] = -1;	
	
    for (i = 1; i <= m; i++){
    	sizem = sizem/2;
    	nmod = nmod*2;
    	count = 0;
		while(listsyn[count][0] != -1){  									
			n1 = 1; n2 = 1; mod1 = 0; mod2 = 0;			
			while(listsyn[count][1] > mod2 + sizem){
				mod2 += sizem;
				n2++;	
			}												
			while(listsyn[count][0] > mod1 + sizem){
				mod1 += sizem;
				n1++;	
			}				
			if (n1!=n2 & listn->at(listsyn[count][0])->getTypeSyn() == 0){	
				nn = mod1 + (int) (randNext()*sizem);
				listsyn[count][1] = nn;				
			}
			if (n1!=n2 &  listn->at(listsyn[count][0])->getTypeSyn() == 1){	
				sort = randNext();
				if (sort < probEx) {
					nn = mod1 + (int) (randNext()*sizem);
					listsyn[count][1] = nn;	
					listsyn[count][2] = 1;
				}
				else {
					listsyn[count][2]++;
				}											
			}							  	
	    	count++;
    	}
    }
    
    count = 0;
    while(listsyn[count][0]!=-1){
    	syn = new alphasynapse;
    	lists->push_back(syn);
    	if( listn->at(listsyn[count][0])->getTypeSyn()==0)
    		syn->setpar(gin, 1, false);
    	else{
	    	syn->setpar(gexc, 1, true);
    	}	
    		    	
	     listn->at(listsyn[count][0])->makeconnection( listn->at(listsyn[count][1]), syn);		    	
    	count++;
    }
    
    nameFile.str("");
	nameFile << "listsyn_M" << m<<"_In_"<<type_inh<<"_Ex2_"<< type_ex2<<"_Sim_" << sn <<"_NumSim_" << NumSim <<  ".dat";
	pfile = fopen((nameFile.str()).c_str(),"w");
	for (i = 0; i <= count-1; i++){
	  fprintf(pfile,"%d\t%d\n",listsyn[i][0],listsyn[i][1]);
	} 
	fclose(pfile);
       
	for (i=0; i<1000000; i++)
		 delete[] listsyn[i];
	delete[] listsyn;
	delete[] ex12in;
}   

void simulation::rasterdata(char * pchar) {
    ofstream fout, fplot;
    fout.open(pchar);
    int i, j,k, max_size = 0;
    
    for (i = 0; i < listn->size(); i++)
        if (!listn->at(i)->getevents()->empty())
        	if (listn->at(i)->getevents()->size() > max_size)
        		max_size = listn->at(i)->getevents()->size();    
	
    for (i = 0; i < listn->size(); i++) {
        if (!listn->at(i)->getevents()->empty()) {
        	fout << i << "\t";
            for (j = 0; j < listn->at(i)->getevents()->size(); j++) {
                fout << listn->at(i)->getevents()->at(j) << "\t";
            }
            if (listn->at(i)->getevents()->size() < max_size){
            	k = max_size - listn->at(i)->getevents()->size();
            	for (j = 0;j < k; j++)
            		fout << 0.0 << "\t";	
            }
            fout << "\n";
        }
    }    
}

void simulation::LifeTime(char * pchar) {
    ofstream fout;
    fout.open(pchar);
    fout<<LastSpike()-tStim-TimeFreeFall-dtStepsTraj*dtTraj-tStimAfterChange<<"\t%%LIFETIME after initial External Stimulation AND Free Fall AND dT anlong Traj AND dT Perturbation\n\n";
    
    fout<<"\n%Time of External Stimulation [ms]: "<< tStim;
    fout<<"\n%Amplitude of the External Current (integer): "<< InpAmp;
    fout<<"\n%Number of neurons receiving external stimulus starting from the first neuron upwards: "<< n_stim;
    fout<<"\n%Inhibitory Conductance: "<< gin;
    fout<<"\n%Excitatory Conductance: "<< gexc;
    
    fout<<"\n\n% Time of Free Fall after initial stimulation "<<TimeFreeFall ;
    fout<<"\n% dT along trajectory after free fall "<< dtTraj;
    fout<<"\n% # of steps along trajectory after free fall "<< dtStepsTraj;
    fout<<"\n%=> Time along the trajectory after the free fall "<< dtStepsTraj*dtTraj;
    
    
    fout<<"\n\n%DeltaT Perturbation "<< tStimAfterChange;
    fout<<"\n%Amplitude of Perturbation: "<< InpAmpAfterChange;
    fout<<"\n%Secondary Perturbed Neurons: "<< changeNneurons;
    
    fout<<"\n\n%SIMULATION NUMBER: "<< sn;
    fout<<"\n%IntegrPrecisioin: "<< h;
    fout<<"\n%Total time [ms]: "<< tmax;
    fout<<"\n%Total number of neurons: "<< nneuron;
    fout<<"\n%Connection Probability: "<< prob;
    fout<<"\n%Rewiring prob for exc neurons: "<< probEx;
    fout<<"\n%Perc of inh neurons: "<< perc_inh;
    fout<<"\n%Type of inh neurons: FS  LTS RS: "<< type_inh;
    fout<<"\n%Type I of exc neurons: RS  IB  CH: "<< type_ex1;
    fout<<"\n%Type II of exc neurons: RS  IB  CH: "<< type_ex2;
    fout<<"\n%Hierarchical levels: "<< m;
    fout<<"\n%Perc of type I exc neurons: "<< perc_ex1;
    fout<<"\n%Perc of type II exc neurons: "<< perc_ex2;
}

void simulation::printTime(string name){
	FILE * pfile;
	pfile = fopen(name.c_str(), "a");
	izhi * nrn;
	double tmp = 0.0, tmp2 = 0.0, aux = 0.0, last = 0.0, std = 0.0;
	for (int i = 0; i < nneuron; i++){
		nrn = listn->at(i);
		if ((nrn->getevents())->empty()) continue;
		aux = (nrn->getevents())->back();
		tmp += aux;
		tmp2 += aux*aux;
		if (aux > last) last = aux;
	}
	tmp = tmp/(listn->size());
	tmp2 = tmp2/(listn->size());
	std = sqrt(tmp2 - tmp*tmp);
	fprintf(pfile,"%lf\t%lf\t%lf\t%lf\t%lf\n",gin, gexc, last, tmp, std);
	fclose(pfile);
}  	

