/* Header File containing parameters for the discriminator code. 
   Required by netinmain.cpp, discrimnet.cpp and structure.cpp */

using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

string createString(string, int, int);
const int NEsub = 400; // cells per excitatory (E-)population in memory network
const int NIsub = 200; // cells per inhibitory (I-)population
const int NEpops = 12; // number of graded bistable E-populations
const int NIpops0 = 1; // number of I-populations
const int Nreadout = 400; // number of cells in Readout Population
const int NEbins = 12; // number of populations bins for memory cell output
const int NIbins = 1; // number of bins for interneuron output
const int NE = NEsub*NEpops;
const int NI = NIsub*NIpops0;
const int Ngates = 3; // detail for facilitating synapses
const int Ncues = 3; // initial cue of zero followed by f1 then f2

const double pi = 4.0*atan(1.0); 
  
double Ran_Gaussian();

class Params{
 public:
  /* Single neuron parameters */
  double cm;       
  double vreset;
  double vl;
  double vth0;
  double gl0;
  double tref;
  double gExt;
  double gExtI;
  double gApp;
  double tauExt;
  double rateExt;
  double gSynE;
  double gSynI;
  double vSynE;
  double vSynI;

  /* Synaptic parameters */
  double tau_f[Ngates]; // Set of time constants for facilitation
    double C_f[Ngates]; // Sets the scale of facilitation (Matveev and Wang)
  double tau_d; //time constant for synaptic depression 
  int N0;  // number of vesicles in synapse

  double vary; // variable to allow heterogeneity in synapses

};

/* Each neuron in the network is a member of the class cell 
   with the following set of parameters and variables */
class CELL{
  double Fired();
  static long seed; 
  static long seed2;
  double vold;
  double vnew;
 public:
  Params *pc;
  CELL() {  }
  ~CELL(){  }
  void datain();
  void Init(Params*,double,double);
  double findV(double,double);
  double vnow;
  double tm;
  double vth;
  double gl;
  double tref;
  double cellInh;
  double lftime;
  double llftime;
  double sSynE;
  double sSynI;
  double sExt;
  double sExtI;
  double sApp;
  double sAppFact;
  double gtilde;
  double gtildebar;
  double v0;
  double v0bar;
  double vshadow;
  double vshadowbar;
  double ISIbar;
  double binrate;
  int numfired;
  double pv;

  /* Fired is a subroutine to interpolate precise time of spike 
     when the membrane potential crosses threshold */

  double Fired(double v1, double v2, double time, double dt){
    double fdt = dt*(v2-vth)/(v2-v1);
    return time - fdt;
  }


  /* SExt gives the time-varying external conductance 
     (when multiplied by gExt) due to Poisson noise */
  double SExt(double, double, double, double);

  /* Integrate updates the membrane potential for eahc neuron */
  void Integrate(double, double);

  /* FacDep calculates the synaptic strength which can vary 
     because of depression and facilitation */
  double FacDep();

  /* Ogate is the state of each of 3 facilitating channels */
  double Ogate[Ngates];
  double Nves;   // number of vesicles

  /* EGetParams() and IGetParams() help initialize excitatory 
     and inhibitory cells respectively */
  void IGetParams();
  void EGetParams();
  
  /* cellfire is TRUE if the cell fired in the time increment */
  bool cellfire;
};


/* OK -- bad C++ style, but this is a set of parameters needed 
   throughout the code  */

class NET{
  Params Epar;
  Params Ipar;
 public:
  /*  Totalnumfired counts total spikes and acts as a safety 
      valve, as simulation stops if Totalnumfired > maxfired */
  int Totalnumfired;
  double dt;
  double gApp;
  double gInit;
  double gEE;
  double gEI;
  double gIE;
  double gII;
  double gNMDAfact;
  double gAMPAfact;
  double gNMDAfactback;
  double gAMPAfactback;
  double gNMDAfactfwd;
  double gAMPAfactfwd;
  double Ealpha;
  double Ialpha;
  double WEE[NE][NEpops];
  double WEI[NE][NIpops0];
  double WIE[NI][NEpops];
  double WII[NI][NIpops0];
  double WEEpop[NEpops][NEpops];
  double WEIpop[NEpops][NIpops0];
  double WIEpop[NIpops0][NEpops];
  double WIIpop[NIpops0][NIpops0];
  double Eaxonfact[NE];
  double Iaxonfact[NI];
  double Eaxonfact2[NE];
  double Iaxonfact2[NI];
  double sigma_EE;
  double sigma_EI;
  double sigma_IE;
  double sigma_II;
  double WER[NEpops];
  double WIback[NIpops0];
  double WEback[NEpops];
  double WEforward[NEpops];
  double WIforward[NEpops];
  double WIRO;
  double WEEcross;
  double WEIcross;
  double WIEcross;
  double WIIcross;  
  double WEEin;
  double WEEon;
  double WEon;
  double ECm;
  double EVth;
  double EVreset;
  double EVl;
  double Egl;
  double Etref;
  double Etaus;
  double EtauAMPA;
  double EgExt;
  double EgExtI;
  double EgApp;
  double rApp0;
  double AppFact;
  double EtauExt;
  double ErateExt;
  double ICm;
  double IVth;
  double IVreset;
  double IVl;
  double Igl;
  double Itref;
  double Itaus;
  double IgExt;
  double IgExtI;
  double ItauExt;
  double IrateExt;
  double VSynE;
  double VSynI;
  double Vthchange;
  double glchange;
  double VRshift;
  double glRshift;
  double Vinshift;
  double glinshift;
  double Vonshift;
  double glonshift;
  double gonApp;
  double Ethresh_shift;
  double Ithresh_shift;
  double Evary;
  double Ivary;
  double maxfired;
  double maxrate;
  double bintime;
  CELL *ECell;
  CELL *ICell;
  CELL *EinCell;
  CELL *EonCell;
  double *sE;
  double *sI;
  double *sEin;
  double *sEon;
  double *sEPlus;
  double *sIPlus;
  double *sEinPlus;
  double *sEonPlus;
  double *sAMPA;
  double *sAMPAPlus;
  double *sAMPAin;
  double *sAMPAinPlus;
  double *sAMPAon;
  double *sAMPAonPlus;
  NET(){};
  ~NET(){};
  void IGetParams();
  void EGetParams();
  void Datin(int,int);
  void WDatin();
  void Init();
  void GatherInputs();
  void Connect();
  double sigmaDistrib;
  double Distrib(int,int,int,int,double);
  double cuetime[Ncues];
  double cuelength[Ncues];
  double cuestrength1;
  double cuestrength2; 
  double Einresponse[Ncues];
  double Eresponse[NEbins+1][Ncues];
  double Iresponse[NIbins][Ncues];
  int rcount[Ncues];
  double range;
  double sApp0;
  ofstream logOut;
  ofstream Efiredout;
  ofstream Ifiredout;
  ofstream Eout[NEbins+1];
  ofstream Iout[NIbins];
  ofstream Einout;
  ofstream Eonout;
};



