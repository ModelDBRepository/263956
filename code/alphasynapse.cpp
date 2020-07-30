#include "alphasynapse.h"
namespace synapse {
	
    alphasynapse::alphasynapse() {
        gmax = 1.0;
        delay = 1;
        spikes = new vector <double>;
        gsyn = 0.0;
        type = true; //0 = inh, 1 = exc      
    }

	alphasynapse::~alphasynapse(){
     	delete spikes;
	}
		
    void alphasynapse::setpar(double gmaxArg, unsigned char delayArg, bool typeArg) {
        gmax = gmaxArg;
        delay = delayArg;
        type = typeArg;
    }
	
    void alphasynapse::addevent(double spk) {
		spikes->push_back(spk);
    }
    
    void alphasynapse::evaluate(double time) {	 
    	gsyn = 0;   	
	    if (!spikes->empty()){
	       	double s = time - spikes->at(0) - delay;          
	       	if (s>=0){
	       		gsyn += gmax;
	       		spikes->erase(spikes->begin());
	       	}
	   	}    
    }

    double alphasynapse::getGsyn() {
        return gsyn;
    }
	
    bool alphasynapse::getType() {
        return type;
    }
    
}

