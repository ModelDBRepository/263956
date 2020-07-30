#ifndef _ALPHASYNAPSE_H
#define	_ALPHASYNAPSE_H

#include <vector>
#include <stdlib.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "izhi.h"
#include <random>
#include <iostream>
#include <string>

namespace Cells{
	class izhi;}
using namespace std;
using namespace Cells;

namespace synapse {
	
    class alphasynapse {
    public:
        alphasynapse();
        ~alphasynapse();
        void setpar(double gmaxArg, unsigned char delayArg, bool typeArg);
        void evaluate(double time);
        void addevent(double spk);
        bool getType();
        double getGsyn();
     protected:    	
        vector <double> * spikes;       
        unsigned char delay;
        double gmax;
        double gsyn;
        bool type;
   };
}

#endif	/* _ALPHASYNAPSE_H */

