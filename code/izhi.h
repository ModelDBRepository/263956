#ifndef _IZHI_
#define _IZHI_
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <complex>
#include <vector>
#include <algorithm>
#include "alphasynapse.h"
#include <chrono>

#define a_21 (1/5)
#define a_31 (3/40)       
#define a_32 (9/40)
#define a_41 (44/45)
#define a_42 (-56/15)
#define a_43 (32/9)
#define a_51 (19372/6561)
#define a_52 (-25360/2187)
#define a_53 (64448/6561)
#define a_54 (-212/729)
#define a_61 (9017/3168)
#define a_62 (-355/33)
#define a_63 (46732/5247)
#define a_64 (49/176)
#define a_65 (-5103/18656)
	     
#define b_1 (35/384)
#define b_2 (0)
#define b_3 (500/1113)
#define b_4 (125/192)
#define b_5 (-2187/6784)
#define b_6 (11/84)
#define bp_1 (5179/57600)
#define bp_2 (0)
#define bp_3 (77571/16695)
#define bp_4 (393/640)
#define bp_5 (-92097/339200)
#define bp_6 (187/2100)

#define atol (10^(-2))
#define rtol (10^(-4))

namespace synapse{
	class alphasynapse;}
using namespace synapse;
using namespace std;
namespace Cells {
	
    class izhi {
    public:
        izhi();  
        ~izhi();     
        void evaluate (double inj,double time);
        void setw0(double w0, double w1);
        double getw(int ind);
        void seth(double hArg);
        void setpar(double aArg, double bArg, double cArg, double dArg);
        void setequilibrium();
        void makeconnection(izhi * dend, alphasynapse * syn);
        void addsyndend(alphasynapse * syn);
        vector<double> * getevents();
        void setExcitatory();
        void setInhibitory();
        int getTypeSyn(); 
	    vector <alphasynapse*> * sdend;
	    void set_noiseAmplitude(double amp, int IC);
		void updateNeNi();
		double getCurr();
    protected:
    	void startKs();
        void fx(double inj,double time);
        void sendevent(double time);
        vector <alphasynapse*> * saxon;
        vector <double> * events;
        double * kv, * ku, * kgin, * kgex;
        double * fbuf;
        double * w, * waux;
        double a, b, c, d;
        double h;
        int typesyn;
        //gsl_rng *r_aloc;
        double noiseAmplitude, D;
        double incrementEx, incrementIn;
        double NexCon, NinCon;
        double tauex, tauin, Eex, Ein;
        double neuronCurr;
    };
}
#endif
