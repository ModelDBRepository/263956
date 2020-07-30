#include "izhi.h"
namespace Cells {
	
    izhi::izhi() {
	fbuf = new double[2];
        w = new double[4];
        waux = new double[4];
        setpar(0.02, 0.2, -65.0, 8.0);
        seth(0.1);
        setequilibrium();
        saxon = new vector <alphasynapse*>;
        sdend = new vector <alphasynapse*>;
        events = new vector <double>;
        kv = new double[2]; 
        ku = new double[2]; 
        kgex = new double[2]; 
        kgin = new double[2]; 
        startKs();
        /* A pointer to a gls random number generator object */
	  	//r_aloc = gsl_rng_alloc(gsl_rng_mt19937);
	  	noiseAmplitude=0;
	  	NexCon=0; NinCon=0;	
	    incrementEx=0; incrementIn=0;
	    tauex=5; tauin=6; Eex=0; Ein=-80;
	    D=0;
	    neuronCurr=0;
    }
	    
    izhi::~izhi(){
    	delete[] fbuf;
    	delete[] w;
    	delete[] waux;
    	delete[] kv;
    	delete[] ku;
    	delete[] kgex;
    	delete[] kgin;
    	delete saxon;
    	delete sdend;
    	delete events; 	
    }     
	
    void izhi::setw0(double w0, double w1) {
        w[0] = w0;
        w[1] = w1;
        w[2] = 0;
        w[3] = 0;
    }
	
    void izhi::seth(double hArg) {
        h = hArg;
    }
	
    double izhi::getw(int ind) {
        return w[ind - 1];
    }
	     
    void izhi::setpar(double aArg, double bArg, double cArg, double dArg) {
        a = aArg;
        b = bArg;
        c = cArg;
        d = dArg;
        setequilibrium();
    }
	
    void izhi::setequilibrium() {
        complex <double> ac(a, 0.0), bc(b, 0.0), cc(c, 0.0), dc(d, 0.0);
        complex <double> eq1 = (b - 5.0 + sqrt((5.0 - b)*(5.0 - b) - 22.4)) / 0.08;
        complex <double> eq2 = (b - 5.0 - sqrt((5.0 - b)*(5.0 - b) - 22.4)) / 0.08;
        complex <double> tr = 0.08 * eq1 + 5.0 - a;
        complex <double> det = -a * (0.08 * eq1 + 5.0) + a*b;
        complex <double> lambda1 = 0.5 * (tr + sqrt(tr * tr - 4.0 * det));
        complex <double> lambda2 = 0.5 * (tr - sqrt(tr * tr - 4.0 * det));
        if ((double) lambda1.real() < 0.0 && (double) lambda2.real() < 0.0)
            setw0((double) eq1.real(), b * (double) eq1.real());
        else
            setw0((double) eq2.real(), b * (double) eq2.real());
    }
	
    void izhi::fx(double inj, double time) {
        double isyn = 0.0;
        isyn += waux[2]*(waux[0] - Eex) + waux[3]*(waux[0] - Ein);
        inj -= isyn;
        fbuf[0] = (0.04 * waux[0] * waux[0] + 5 * waux[0] + 140 - waux[1] + inj);
        fbuf[1] = a * (b * waux[0] - waux[1]);
        neuronCurr = inj;
    }
    
    double izhi::getCurr(){
    	return neuronCurr;		
    }
    
    void izhi::set_noiseAmplitude(double amp, int IC){
		noiseAmplitude = sqrt(amp*h);
		//gsl_rng_set(r_aloc,IC); //seed
		D=amp;
	}
	
    void izhi::updateNeNi(){
   		for(int i = 0; i < sdend->size(); i++) {
			if((sdend->at(i))->getType()) //excitatory
				NexCon++;
			else
				NinCon++;
		}
	}
	
    void izhi::evaluate(double inj, double time) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator (seed);
  		  std::normal_distribution<double> distribution1(0,sqrt(NexCon));
                  std::normal_distribution<double> distribution2(0,sqrt(NinCon));	
		waux[0] = w[0]; //v
		waux[1] = w[1];	//u
		waux[2] = w[2]; //gex
		waux[3] = w[3]; //gin   
		if(noiseAmplitude!=0){	    
		    double n1 = distribution1(generator); 
		    double n2 = distribution2(generator);
	 	    incrementEx = noiseAmplitude*n1; 
		    incrementIn = noiseAmplitude*n2; 
            //std::normal_distribution<double> distribution (0.0,NexCon);   //gsl_ran_gaussian(r_aloc,NexCon);
            //std::normal_distribution<double> distribution (0.0,NinCon); //gsl_ran_gaussian(r_aloc,NinCon); 
        }
	
		kgex[0] = - waux[2]/tauex; //f0 exc
		kgin[0] = - waux[3]/tauin; //f0 inh			            
        fx(inj, time);     		
        kv[0] = fbuf[0]; //f0
        ku[0] = fbuf[1]; //f0
        
		waux[0] = w[0] + h*kv[0]; //x1
		waux[1] = w[1] + h*ku[0];  
		kgex[1] = w[2] + incrementEx + h*kgex[0]; //x1
		kgin[1] = w[3] + incrementIn + h*kgin[0];
		
		waux[2] = kgex[1];
		waux[3] = kgin[1];
        fx(inj, time);
        kv[1] = fbuf[0]; //f(x1)
        ku[1] = fbuf[1]; //f(x1)
	        
        w[0] += h/2*(kv[0]+kv[1]);
        w[1] += h/2*(ku[0]+ku[1]);  
        w[2] += incrementEx + h/2*(kgex[0] - kgex[1]/tauex);
        w[3] += incrementIn + h/2*(kgin[0] - kgin[1]/tauin);  
        
           
        for(int i = 0; i < sdend->size(); i++) {
        	(sdend->at(i))->evaluate(time);
        	
			if((sdend->at(i))->getType()) //excitatory
				w[2] += (sdend->at(i))->getGsyn();
			else
				w[3] += (sdend->at(i))->getGsyn();
		}
		
		if (w[2]<0) w[2]=0;
		if (w[3]<0) w[3]=0;              
		
		if (w[0] >= 30.0) {
		  w[0] = c;
		  w[1] = w[1] + d;
		  sendevent(time);		  
		}
    }
    	
    void izhi::makeconnection(izhi * dend, alphasynapse * syn) {
        saxon->push_back(syn);
        dend->addsyndend(syn);
    }
	
    void izhi::addsyndend(alphasynapse * syn) {
        sdend->push_back(syn);
    }
        	
    void izhi::sendevent(double time) {
        int i = 0;
        events->push_back(time);
        while (i < saxon->size()) {
            (saxon->at(i))->addevent(time);
            i++;
        }
    }
	
    vector<double> *  izhi::getevents(){
        return events;
    }
     
    void izhi::setExcitatory(){
		typesyn = 1;
    }
    
    void izhi::setInhibitory(){
		typesyn = 0;
    }
	
    int izhi::getTypeSyn(){
		return typesyn;
    }    
	
    void izhi::startKs(){
    	for(int i = 0; i < 2; i++){
    		kv[i] = 0.0; 	
    		ku[i] = 0.0; 
    		kgex[i] = 0.0;
    		kgin[i] = 0.0;
    	}
    }
    
}

