/* MAIN DISCRIMINATION C++ CODE 
To run the code: 

When the code is compiled, run the executable followed by a number which 
sets the frequencies of the two stimuli and initializes the random number 
generator. 
The last 4 digits give a number, N, which encodes the frequency, 
any digits before the last 4 are the initial seed to the random 
number generator. 

N/2 is the value of f1.
An odd number for N equates to f2>f1, while an even number 
for f2<f1. 

That is, set N= 2 times f1 +1 for f2 > f1 and 
set N= 2 times f1 for f2 < f1.
See netinmain.cpp. 
 
On a parallel machine use the script file which runs 14 trials in parallel, 
each trial with a different pair of f1,f2. 
The first input runs from 10 to 34 in steps of 4, 
corresponding to f1=10Hz to 34Hz in steps of 4Hz. 
f2 is set as f1+8Hz or f1-8Hz.


Each neuron is a member of class CELL, (within class NET). 

The 5 types of cell are: 
EinCell[i] (i=0 to 399) the C-neurons or comparison cells 
                        in the model which receive input.
EonCell[i] (i=0 to 399) the `ON' neurons which are a single bistable group, 
                        switched on but untuned during the task.
ECell[i] (i=0 to 4799) the 12x400 bistable subpopulations of the memory network.
ECell[i] (i=4800 to 5199) the 400 Readout cells of the memory network.
ICell[i] (i=0 to 199) the 200 interneurons, which feedback the inhibition.

Input files: 

connections_in.dat 
   Contains the sets of weights for all connections in the network.
   Used by code in structure.cpp to form the complete weight-matrix.
discrimnet_in.dat
   Contains all other inputs. 
   First a set of single neuron parameters.
   Finally the set of cue times and durations for the simulation.

Output files:
   Efiredout << goes to "Efired.dat" 
      spike times and ISIs of selected E-neurons
   Ifiredout  << goes to "Efired.dat"
      spike times and ISIs of selected I-neurons
   Eout[i] << goes to "Eout"i+1".dat" 
      where i = 0 to 11 for the 12 bistable memory populations 
      and i = 12 for the readout cells.
      Outputs time-binned population average firing rate.
   Eout[i] << goes to "Iout"i+1".dat" 
      where i = 0 for the single interneuron population.
      Outputs time-binned population average firing rate.
   Einout << goes to "Einout.dat" 
      Outputs average time-binned rates of C-neurons.
   Eonout << goes to "Eonout.dat" 
      Outputs average time-binned rates of `ON' neurons.
*/

#include "MersenneTwister.h" //random number generator
#include "dnet.h"

extern MTRand rand1; 
extern MTRand rand2;


void CELL::Init(Params* extpc, double vthadd=0.0, double gladd =0.0){
  /* Initialization of each neuron extpc points to the values for the 
     particular cell type (excitatory/inhibitory) */

  vth = extpc->vth0 + vthadd; // voltage threshold
  gl = extpc->gl0 + gladd; // leak conductance


  /* The following section sets a random time since the last spike 
     and an initial membrane potential accordingly, based on an 
     assumed spontaneous rate of avrate */
  double avrate = 1.0;
  lftime = -rand1()/avrate;
  if(lftime > -extpc->tref)
    vnow = extpc->vreset;
  else
    vnow = (vth - extpc->vreset)*(-lftime)*avrate + extpc->vreset;
  vold = vnow;
  vnew = vnow;

  numfired = 0;

  /* This section sets initial value of external noise input based 
     one estimated time since last spike. */
  double wait = -log(rand1())/extpc->rateExt;
  sExt = exp(-wait/extpc->tauExt);

  /* These variables ending in -bar are to calculate means across time */
  gtildebar = 0.0;
  v0bar = 0.0;
  vshadowbar = 0.0;
  ISIbar = 0.0;

  /* pc = present cell */
  pc = extpc;
  binrate = 0.0;

  Nves = extpc->N0; // number of vesicles in axon terminals for this cell type

  /* cellInh adds leak conductance, since the leak potential 
     is equalt to the inhibitory reversal potential in this model. 
     Default for parameter "vary" is zero, but when non-zero gives 
     heterogeneity in leak conductances */

  cellInh = Ran_Gaussian() * extpc->vary; 

}


double CELL::findV(double tau, double deltat){
  /* Basic temporal integration */

  /* RUNGE KUTTA-4
  double k1 = deltat*(-(vnow-v0)/tau);
  double k2 = deltat*(-(vnow+0.5*k1-v0)/tau);
  double k3 = deltat*(-(vnow+0.5*k2-v0)/tau);
  double k4 = deltat*(-(vnow+k3-v0)/tau);

  return vnow + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;
  */
  /* RUNGE KUTTA-2 */
  double k1 = deltat*(-(vnow-v0)/tau);
  double k2 = deltat*(-(vnow+0.5*k1-v0)/tau);

  return vnow + k2;

}

double CELL::SExt(double sext,double rate,double tauext, double dt){
  /* Updates the background noise input */

  bool extfire = false;
  do
    if ( 1.0 - exp(-dt*rate) > rand1() ){
      double delay = dt*rand1();
      sext *= exp(-delay/tauext);
      sext +=1.0;
      extfire = true;
      dt = dt - delay;
    }
    else{
      sext *= exp(-dt/tauext);
      extfire = false;
    }
  while (extfire == true);

  return sext;
}

void CELL::Integrate(double dt,double t){
  /* Calculates currents needed for time integration of 
     the membrane potential, then updates the membrane potential 
     and checks if cell has fired. 
     The shadow voltage is used (a membrane potential with no 
     reset) as an approximation to the dendritic membrane 
     potential used when calculating NMDA conductance 
  */

  sExt = SExt(sExt,pc->rateExt,pc->tauExt,dt);
  sExtI = SExt(sExtI,pc->rateExt,pc->tauExt,dt);

  gtilde = gl + pc->gExt*sExt + pc->gExtI*sExtI 
    + pc->gSynE*sSynE + pc->gSynI*sSynI + pc->gApp*sApp*sAppFact;
  gtildebar += gtilde;
  
  v0 = (pc->vl*gl + (pc->gExtI*sExtI + pc->gSynI*sSynI)*pc->vSynI) /gtilde;
  v0bar += v0;
  tm = pc->cm / gtilde;

  vold = vnow;
  if ( t > lftime+pc->tref ){
    if(t-dt >= lftime+pc->tref){
      vnew = findV(tm,dt);
    }
    else{
      vnew = findV(tm,t-lftime-pc->tref);
    }
  }
  
  vshadow = v0; // shadow voltage for NMDA conductance is not reset
  vshadowbar += vshadow;

  /* If threshold is crossed, calculate exact crossing time 
     and reset for a time given by refractory period, tref */

  if (vnew >= vth ){
    cellfire = true;
    llftime = lftime;
    lftime = Fired(vold,vnew,t,dt);
    
    numfired++;
    binrate++;
    ISIbar += lftime-llftime;
    vnew = pc->vreset;

    if ( t-lftime > pc->tref ) {
      vnow = pc->vreset;
      vnew = findV(tm,t-lftime-pc->tref);
    }
  }
  
  vnow = vnew;
  
}

double CELL::FacDep(){
  pv = 1.0;
  double deltat = lftime-llftime;
  for (int i=0; i<Ngates; i++){
    Ogate[i] *= exp(-deltat/pc->tau_f[i]);
    Ogate[i] = 1.0 - (1.0-Ogate[i])*pc->C_f[i];
    pv *= Ogate[i];
  }
  double p = 0.0;
  if ( pc->tau_d > 0.0 ) 
    p = exp(-deltat/pc->tau_d);

  Nves += (pc->N0-Nves)*(1.0-p);

  return Nves*pv/pc->N0;
}

void NET::Datin(int istr1, int istr2){
  /* Input all the single cell properties 
     and simulation protocol. 
  */

  cuestrength1 = double(istr1);
  cuestrength2 = double(istr2);

  cout << " cuestrength1 " << cuestrength1 
       << " cuestrength2 " << cuestrength2 << endl;

  ifstream input;
  input.open("discrimnet_in.dat");
  string dummy;
  input >> ECm >> dummy;
  input >> EVth >> dummy;
  input >> EVreset >> dummy;
  input >> EVl >> dummy;
  input >> Egl >> dummy;
  input >> Etref >> dummy;
  input >> Etaus >> dummy;
  input >> EtauAMPA >> dummy;
  input >> EgExt >> dummy;
  input >> EgExtI >> dummy;
  input >> ErateExt >> dummy;
  input >> EtauExt >> dummy;
  input >> ICm >> dummy;
  input >> IVth >> dummy;
  input >> IVreset >> dummy;
  input >> IVl >> dummy;
  input >> Igl >> dummy;
  input >> Itref >> dummy;
  input >> Itaus >> dummy;
  input >> IgExt >> dummy;
  input >> IgExtI >> dummy;
  input >> IrateExt >> dummy;
  input >> ItauExt >> dummy;
  input >> gEE >> dummy;
  input >> gEI >> dummy;
  input >> gIE >> dummy;
  input >> gII >> dummy;
  input >> gNMDAfact >> dummy;
  input >> gAMPAfact >> dummy;
  input >> gNMDAfactback >> dummy;
  input >> gAMPAfactback >> dummy;
  input >> gNMDAfactfwd >> dummy;
  input >> gAMPAfactfwd >> dummy;
  input >> VSynE >> dummy;
  input >> VSynI >> dummy;
  input >> Ealpha >> dummy;
  input >> Ialpha >> dummy;
  input >> dt >> dummy;
  input >> sigmaDistrib >> dummy;
  input >> Vthchange >> dummy;
  input >> glchange >> dummy;
  input >> VRshift >> dummy;
  input >> glRshift >> dummy;
  input >> Vinshift >> dummy;
  input >> glinshift >> dummy;
  input >> Vonshift >> dummy;
  input >> glonshift >> dummy;
  input >> gonApp >> dummy;
  input >> maxrate >> dummy;
  input >> bintime >> dummy;
  input >> rApp0 >> dummy;
  input >> EgApp >> dummy;
  input >> AppFact >> dummy;
  input >> range >> dummy;
  input >> Ethresh_shift >> dummy;
  input >> Ithresh_shift >> dummy;
  input >> Evary >> dummy;
  input >> Ivary >> dummy;

  for (int j = 0; j < Ncues; j++) 
    input >> cuetime[j] >> dummy;
  cuetime[Ncues] = 999;
  for (int j = 0; j < Ncues; j++) 
    input >> cuelength[j] >> dummy;

  cout << " inputs done " << endl;

  string logString = createString( "net1log",NEsub,3);
  logString += ".dat";
  cout << logString << endl;

  logOut.open(logString.c_str());
 
  logOut << " NE " << NE << endl;
  logOut << " NI " << NI << endl;
  logOut << " Nreadout " << Nreadout << endl;
  logOut << " ECm " << ECm << endl;
  logOut << " EVth " << EVth << endl;
  logOut << " EVreset " << EVreset << endl;
  logOut << " EVl " << EVl << endl;
  logOut << " Egl " << Egl << endl;
  logOut << " Etref " << Etref << endl;
  logOut << " Etaus " << Etaus << endl;
  logOut << " EtauAMPA " << EtauAMPA << endl;
  logOut << " EgExt " << EgExt << endl;
  logOut << " EgExtI " << EgExtI << endl;
  logOut << " ErateExt " << ErateExt << endl;
  logOut << " EtauExt " << EtauExt  << endl;
  logOut << " ICm " << ICm << endl;
  logOut << " IVth " << IVth << endl;
  logOut << " IVreset " << IVreset << endl;
  logOut << " IVl " << IVl << endl;
  logOut << " Igl " << Igl << endl;
  logOut << " Itref " << Itref << endl;
  logOut << " Itaus " << Itaus << endl;
  logOut << " IgExt " << IgExt << endl;
  logOut << " IgExtI " << IgExtI << endl;
  logOut << " IrateExt " << IrateExt << endl;
  logOut << " ItauExt " << ItauExt  << endl;
  logOut << " gEE " << gEE << endl;
  logOut << " gEI " << gEI << endl;
  logOut << " gIE " << gIE << endl;
  logOut << " gII " << gII << endl;
  logOut << " gNMDAfact " << gNMDAfact << endl;
  logOut << " gAMPAfact " << gAMPAfact << endl;
  logOut << " gNMDAfactback " << gNMDAfactback << endl;
  logOut << " gAMPAfactback " << gAMPAfactback << endl;
  logOut << " gNMDAfactfwd " << gNMDAfactfwd << endl;
  logOut << " gAMPAfactfwd " << gAMPAfactfwd << endl;
  logOut << " VSynE " << VSynE << endl;
  logOut << " VSynI " << VSynI << endl;
  logOut << " Ealpha " << Ealpha << endl;
  logOut << " Ialpha " << Ialpha << endl;
  logOut << " dt " << dt << endl;
  logOut << " sigmaDistrib " << sigmaDistrib << endl;
  logOut << " Vthchange " << Vthchange << endl;
  logOut << " glchange " << glchange << endl;
  logOut << " VRshift " << VRshift << endl;
  logOut << " glRshift " << glRshift << endl;
  logOut << " Vinshift " << Vinshift << endl;
  logOut << " glinshift " << glinshift << endl;
  logOut << " Vonshift " << Vonshift << endl;
  logOut << " glonshift " << glonshift << endl;
  logOut << " gonApp " << gonApp << endl;
  logOut << " maxrate " << maxrate << endl;
  maxfired = maxrate*float(NE+NI)*(cuetime[Ncues-1]+5.0);
  logOut << " maxfired " << maxfired << endl;
  logOut << " bintime " << bintime << endl;
  logOut << " Ethresh_shift " << Ethresh_shift << endl;
  logOut << " Ithresh_shift " << Ithresh_shift << endl;
  logOut << " Evary " << Evary << endl;
  logOut << " Ivary " << Ivary << endl;
  logOut << " gApp " << EgApp << endl;
  logOut << " rApp0 " << rApp0 << endl;
  logOut << " AppFact " << AppFact << endl;
  for ( int j = 0; j < Ncues ; j++)
    logOut << " cuetime " << j << " " << cuetime[j] << endl;
  for ( int j = 0; j < Ncues ; j++)
    logOut << " cuelength " << j << " " << cuelength[j] << endl;
  
}

void NET::Init(){
  ECell = new CELL[NE+Nreadout];
  EinCell = new CELL[NEsub];
  EonCell = new CELL[NEsub];
  ICell = new CELL[NI];
  sE = new double[NE+Nreadout];
  sAMPA = new double[NE+Nreadout];
  sEin = new double[NEsub];
  sAMPAin = new double[NEsub];
  sEon = new double[NEsub];
  sAMPAon = new double[NEsub];
  sI = new double[NI];
  sEPlus = new double[NE+Nreadout];
  sAMPAPlus = new double[NE+Nreadout];
  sEinPlus = new double[NEsub];
  sAMPAinPlus = new double[NEsub];
  sEonPlus = new double[NEsub];
  sAMPAonPlus = new double[NEsub];
  sIPlus = new double[NI];
  cout << " IN Init" << endl;
  EGetParams();
  IGetParams();

  for ( int i = 0; i<NE+Nreadout ; i++){
    sEPlus[i] = (1.0 - (1.0 - 0.5*rand1())*exp(-Ealpha))*
      Ealpha*Epar.N0*Etaus/(Ealpha*Epar.N0*Etaus + Epar.tau_d);
    sAMPAPlus[i] = 0.0;
  }

  for ( int i = 0; i<NEsub ; i++) {
    sEinPlus[i] = (1.0 - (1.0 - 0.5*rand1())*exp(-Ealpha))*
      Ealpha*Epar.N0*Etaus/(Ealpha*Epar.N0*Etaus + Epar.tau_d);
    sEonPlus[i] = (1.0 - (1.0 - 0.5*rand1())*exp(-Ealpha))*
      Ealpha*Epar.N0*Etaus/(Ealpha*Epar.N0*Etaus + Epar.tau_d);
  }

  for ( int i = 0; i<NI ; i++)
    sIPlus[i] = 1.0 - (1.0 - 0.5*rand1())*exp(-Ialpha);

  for (int j = 0; j < NEpops; j++){
    double fac;
    if ( NEpops > 1 ) 
      fac = -0.5 + float(j)/float(NEpops-1);
    else 
      fac = 0.0;
    
    for ( int i = 0; i < NEsub ; i++){
      ECell[i+NEsub*j].Init(&Epar,Vthchange*fac,glchange*fac);
    }
  }
  
  for ( int i = 0; i < NE ; i++){
    ECell[i].sAppFact = 0.0;
  }
  /* Now initialize readout cells */
  for ( int i = 0; i < Nreadout ; i++){
    ECell[i+NE].Init(&Epar,VRshift,glRshift);
    ECell[i+NE].sAppFact = 0.0;
  }  
  /* Now initialize input cells */
  for ( int i = 0; i < NEsub ; i++){
    EinCell[i].Init(&Epar,Vinshift,glinshift);
    EinCell[i].sAppFact = 1.0;
    EonCell[i].Init(&Epar,Vonshift,glonshift);
    EonCell[i].sAppFact = 1.0;
  }  

  /* Now initialize interneurons */
  for ( int i = 0; i < NI ; i++){
    ICell[i].Init(&Ipar);
  }


  Efiredout.open("Efired.dat");
  Ifiredout.open("Ifired.dat");
  for ( int i = 0; i<NEbins+1; i++){
    string EString = createString( "Eout",i+1,2);
    EString += ".dat";
    Eout[i].open(EString.c_str());
  }
  Einout.open("Einout.dat");
  Eonout.open("Eonout.dat");

  for ( int i = 0; i<NIbins; i++){
    string IString = createString( "Iout",i+1,2);
    IString += ".dat";
    Iout[i].open(IString.c_str());
  }
      
}

void NET::EGetParams(){
  /* Parameters for excitatory cells */
  Epar.cm = ECm;
  Epar.vth0 = EVth;
  Epar.vreset = EVreset;
  Epar.vl = EVl;
  Epar.gl0 = Egl;
  Epar.tref = Etref;
  Epar.gExt = EgExt;
  Epar.gExtI = EgExtI;
  Epar.gApp = EgApp;
  Epar.tauExt = EtauExt;
  Epar.rateExt = ErateExt;
  Epar.gSynE = gEE;
  Epar.gSynI = gEI;
  Epar.vSynE = VSynE;
  Epar.vSynI = VSynI;
  Epar.vary = Evary;

  Epar.tau_d = 0.2;
  Epar.N0 = 8;
 
  Epar.tau_f[0] = 0.05;
  Epar.tau_f[1] = 0.2;
  Epar.tau_f[2] = 2.0;
  Epar.C_f[0] = 0.1;
  Epar.C_f[1] = 0.2;
  Epar.C_f[2] = 0.4;
  

}

void NET::IGetParams(){
  /* Parameters for inhibitory cells */
  Ipar.cm = ICm;
  Ipar.vth0 = IVth;
  Ipar.vreset = IVreset;
  Ipar.vl = IVl;
  Ipar.gl0 = Igl;
  Ipar.tref = Itref;
  Ipar.gExt = IgExt;
  Ipar.gExtI = IgExtI;
  Ipar.tauExt = ItauExt;
  Ipar.rateExt = IrateExt;
  Ipar.gSynE = gIE;
  Ipar.gSynI = gII;
  Ipar.vSynE = VSynE;
  Ipar.vSynI = VSynI;
  Ipar.vary = Ivary;
}

void NET::GatherInputs(){

  int NIpops = NIpops0;

  double sum = 0.0;
  double sum2 = 0.0;
  double t = 0.0;
  double bint;
  double rApp = 0.0;
  double rApp2 = 0.0;

  double cellbin[11];
  double cell2bin[11];
  for (int i=0; i<11; i++) cellbin[i] = 0.0;
  for (int i=0; i<11; i++) cell2bin[i] = 0.0;
  double EpopInh[NEpops+1];
  double EpopinInh;

  if (NEpops > 1 ){
    for (int i=0; i<NEpops; i++)
      EpopInh[i] = Ethresh_shift + range*float(i)/float(NEpops-1);
  }
  else{
    for (int i=0; i<NEpops; i++)
      EpopInh[i] = Ethresh_shift;
  }
  
  /* Now do readout population */
  EpopInh[NEpops] = Ethresh_shift;
  EpopinInh = Ethresh_shift;

  ofstream Eglout;
  Eglout.open("Ethresh.dat");
  double Eglbar[NEpops];
  for (int i=0; i<NEpops; i++){
    Eglbar[i] = 0.0; 
    for (int k=0; k<NEsub; k++){
      ECell[i*NEsub+k].cellInh += EpopInh[i];
      Eglout << i*NEsub+k << " " << 
	Egl + ECell[i*NEsub+k].cellInh * gEI << endl;
      Eglbar[i] += Egl + ECell[i*NEsub+k].cellInh * gEI;
    }
    Eglbar[i] /= float(NEsub);
  }

  for (int i=0; i<NEpops; i++)
    Eglout << i << " " <<  Eglbar[i] << endl;

  ofstream Iglout;
  Iglout.open("Ithresh.dat");
  double Iglbar[NIpops0];
  for (int i=0; i<NIpops; i++){
    for (int k=0; k<NIsub; k++){
      ICell[i*NIsub+k].cellInh += Ithresh_shift;
      Iglout << i << " " << 
	Igl + ICell[i*NIsub+k].cellInh * gII << endl;
      Iglbar[i] += Igl + ICell[i*NIsub+k].cellInh * gII;
    }
    Iglbar[i] /= float(NIsub);
  }
  for (int i=0; i<NIpops; i++)
    Iglout << i << " " << Iglbar[i] << endl;

  cout << " Gather Inputs " << endl;
  int cue = 0;
  ofstream cellout;
  cellout.open("cellout.dat");
  ofstream cell2out;
  cell2out.open("cell2out.dat");
  
  double sEpop[NEpops];
  double sAMPApop[NEpops];
  double sIpop[NIpops];
  double sEinpop = 0.0;
  double sAMPAinpop = 0.0;
  double sEonpop = 0.0;
  double sAMPAonpop = 0.0;
  double sEbin[NEpops];
  double sIbin[NIpops];
  double sEinbin = 0.0;
  double sEonbin = 0.0;
  double sROpop = 0.0;
  double sROAMPApop = 0.0;

  for (int i = 0; i<NEpops; i++){
    sEpop[i] = 0.0;
    sEbin[i] = 0.0;
  }

  for (int i = 0; i<NIpops; i++){
    sIpop[i] = 0.0;
    sIbin[i] = 0.0;
  }

  for (int i = 0; i < Ncues; i++) rcount[i]= 0;
  int cellcount = 0;

  double Epopreadin = 0.0;
  
  while ( (t< cuetime[2] + 1.0 ) && 
	  (Totalnumfired < maxfired) ){
    
    bint = 0.0;

    while ( bint < bintime ){
     
      t += dt;
      bint += dt;
      
      if ( ( t>cuetime[cue] ) && 
	   (t < cuetime[cue] + cuelength[cue] ) ) {
	rApp = 0.0;
	if ( cue == 1 )
	  rApp = rApp0*cuestrength1;
	if ( cue == 2 ) 
	  rApp = rApp0*cuestrength2;

	rApp2 = 0.0;
      }
      else{
	rApp = 0.0;
	rApp2 = 0.0;
      }
      
      if ( cue < Ncues-1 ){  
	if ( t > cuetime[cue+1] - 1  ) cue++;
      }

      /* First do input cells */

      for ( int i = 0; i<NEsub; i++ ){
	if (EinCell[i].cellfire == true ){
	  EinCell[i].cellfire = false;
	  double kicksize = EinCell[i].FacDep();
	  
	  EinCell[i].Nves -= kicksize;
	  
	  sEin[i] = sEinPlus[i]
	    *exp(-(EinCell[i].lftime-EinCell[i].llftime) / Etaus );
	  sEinPlus[i] = sEin[i] + kicksize*(1.0-exp(-Ealpha))*(1.0-sEin[i]);
	  sAMPAin[i] = sAMPAinPlus[i]
	    *exp(-(EinCell[i].lftime-EinCell[i].llftime) / EtauAMPA );
	  sAMPAinPlus[i] = sAMPAin[i] + kicksize;
	}
  	sEin[i] = sEinPlus[i]*exp(-(t-EinCell[i].lftime) / Etaus );	
	sAMPAin[i] = sAMPAinPlus[i]*exp(-(t-EinCell[i].lftime)/EtauAMPA );
      }
      sEinpop = 0.0;
      sAMPAinpop = 0.0;
      for ( int i=0; i<NEsub; i++){
	sEinpop += sEin[i];
	sAMPAinpop += sAMPAin[i];
      }
      sEinpop /= double(NEsub);
      sAMPAinpop /= double(NEsub);

      for ( int i = 0; i<NEsub; i++ ){
	if (EonCell[i].cellfire == true ){
	  EonCell[i].cellfire = false;
	  double kicksize = EonCell[i].FacDep();
	  
	  EonCell[i].Nves -= kicksize;
	  
	  sEon[i] = sEonPlus[i]
	    *exp(-(EonCell[i].lftime-EonCell[i].llftime) / Etaus );
	  sEonPlus[i] = sEon[i] + kicksize*(1.0-exp(-Ealpha))*(1.0-sEon[i]);
	  sAMPAon[i] = sAMPAonPlus[i]
	    *exp(-(EonCell[i].lftime-EonCell[i].llftime) / EtauAMPA );
	  sAMPAonPlus[i] = sAMPAon[i] + kicksize;
	}
  	sEon[i] = sEonPlus[i]*exp(-(t-EonCell[i].lftime) / Etaus );	
	sAMPAon[i] = sAMPAonPlus[i]*exp(-(t-EonCell[i].lftime)/EtauAMPA );
      }
      sEonpop = 0.0;
      sAMPAonpop = 0.0;
      for ( int i=0; i<NEsub; i++){
	sEonpop += sEon[i];
	sAMPAonpop += sAMPAon[i];
      }
      sEonpop /= double(NEsub);
      sAMPAonpop /= double(NEsub);

      for(int i = 0; i<NE+Nreadout; i++){
	
	if (ECell[i].cellfire == true ){
	  ECell[i].cellfire = false;
	  double kicksize = ECell[i].FacDep();
	  
	  ECell[i].Nves -= kicksize;
	  
	  sE[i] = sEPlus[i]
	    *exp(-(ECell[i].lftime-ECell[i].llftime) / Etaus );
	  sEPlus[i] = sE[i] + kicksize*(1.0-exp(-Ealpha))*(1.0-sE[i]);
	  sAMPA[i] = sAMPAPlus[i]
	    *exp(-(ECell[i].lftime-ECell[i].llftime) / EtauAMPA );
	  sAMPAPlus[i] = sAMPA[i] + kicksize;

	  //	  if ( i == NE+1 ) 
	  //	    cout << " Fire sAMPA " << i << " " << sAMPA[i] << " " << sAMPAPlus[i] << 
	  //	      " sE " << sE[i] << " kicksize " << kicksize << endl;

	  if(i == 1) {
	    cellbin[9] += kicksize;
	    cellbin[10] += ECell[i].pv;
	  }
	  if(i == 4*NEsub+1) {
	    cell2bin[9] += kicksize;
	    cell2bin[10] += ECell[i].pv;
	  }
	}
  	sE[i] = sEPlus[i]*exp(-(t-ECell[i].lftime) / Etaus );
	sAMPA[i] = sAMPAPlus[i]*exp(-(t-ECell[i].lftime)/EtauAMPA );
	//	if ( i == NE+1 ) 
	//	  cout << " sAMPA " << i << " " << sAMPA[i] << " " << sAMPAPlus[i] << 
	//	    " sE " << sE[i] << " sEPlus " << sEPlus[i] << 
	//	    " lftime " << ECell[i].lftime << " llftime " << ECell[i].llftime << endl;
      }	

      for(int i = 0; i<NI; i++){
	if (ICell[i].cellfire == true ){
	  ICell[i].cellfire = false;
	  sI[i] = sIPlus[i]
	    *exp(-(ICell[i].lftime-ICell[i].llftime) / Itaus );
	  sIPlus[i] = sI[i]+1.0;
	  }
	sI[i] = sIPlus[i]*exp(-(t-ICell[i].lftime) / Itaus );
      }
      
      for (int j = 0; j<NEpops; j++){
	sEpop[j] = 0.0;
	sAMPApop[j] = 0.0;
	for (int i = 0; i<NEsub; i++){
	    sEpop[j] += sE[i+j*NEsub]*Eaxonfact[i+j*NEsub];	
	  sAMPApop[j] += sAMPA[i+j*NEsub]*Eaxonfact[i+j*NEsub];
	}
	sEpop[j] /= NEsub;
	sAMPApop[j] /= NEsub;
      }

      for ( int i=0; i< NEpops; i++){
      	sEbin[i] += sEpop[i];
      }      

      sROpop = 0.0;
      sROAMPApop = 0.0;
      for (int i = 0; i<Nreadout; i++){
	sROpop += sE[i+NE];	
	sROAMPApop += sAMPA[i+NE];
      }
      sROpop /= Nreadout;
      sROAMPApop /= Nreadout;

      for (int j = 0; j<NIpops; j++){
	sIpop[j] = 0.0;
	for (int i = 0; i<NIsub; i++)
	  sIpop[j] += sI[i+j*NIsub]*Iaxonfact[i+j*NIsub];
	sIpop[j] /= NIsub;
      }

      /* Input cells first */

      //Excitatory input first

      for (int i=0; i<NEsub; i++){
	sum = 0.0;
	sum2 = 0.0;
	for (int j = 0; j<NEpops; j++){
	  sum += WEback[j]*sEpop[j];
	  sum2 += WEback[j]*sAMPApop[j];
	}

	EinCell[i].sSynE = (sum*gNMDAfactback 
			    + sEinpop*WEEin*gNMDAfact 
			    + sEonpop*WEon*gNMDAfactfwd)/
	  ( 1.0 + exp(-62*EinCell[i].vshadow)/3.57 ) 
	  + sum2*gAMPAfactback 
	  + sAMPAinpop*WEEin*gAMPAfact 
	  + sAMPAonpop*WEon*gAMPAfactfwd;

	
	// Then inhibitory input 
	
	sum = 0.0;
	for (int j = 0; j<NIpops; j++)
	  sum += WIback[j]*sIpop[j];

	EinCell[i].sSynI = sum + EinCell[i].cellInh;
      }
      
      /* Then find out if input cell has fired. */
      
      for (int i = 0; i<NEsub; i++){
	//EinCell[i].sApp = ECell[i].SExt(EinCell[i].sApp, rApp, Epar.tauExt, dt);
	EinCell[i].sApp = rApp*Epar.tauExt;
	EinCell[i].Integrate(dt,t);
	
	if ( (EinCell[i].cellfire == true ) ){
	  if ( i%(NEsub/2) == 0 ) {
	    Efiredout << i+NE+Nreadout << " " 
		      << setprecision(8) << EinCell[i].lftime;
	    Efiredout << " " << EinCell[i].lftime-EinCell[i].llftime; 
	    Efiredout << endl;
	  }
	  Totalnumfired ++;
	}

      }

      //Now ON cells //

      for (int i=0; i<NEsub; i++){
	EonCell[i].sSynE = (sEonpop*WEEon*gNMDAfact)/
	  (1.0 + exp(-62*EonCell[i].vshadow)/3.57 ) 
	  + sAMPAonpop*WEEon*gAMPAfact ;
	
	// Then inhibitory input 
	sum = 0.0;
	EonCell[i].sSynI = sum + EonCell[i].cellInh;
      }
      
      /* Then find out if ON cell has fired. */
      
      for (int i = 0; i<NEsub; i++){
	EonCell[i].sApp = gonApp*rApp*Epar.tauExt;
	EonCell[i].Integrate(dt,t);
      }
      
      /* The Memory E-cells */

      //Excitatory input first

      for (int i=0; i<NE; i++){
	double sumNMDA = 0.0;
	double sumAMPA = 0.0;
	for (int j = 0; j<NEpops; j++){
	  sumNMDA += WEE[i][j]*sEpop[j];
	  sumAMPA += WEE[i][j]*sAMPApop[j];
	}
	double sumfwd = WEforward[i/NEsub]*sEinpop;

	ECell[i].sSynE = sumNMDA*gNMDAfact/
	      (1.0 + exp(-62*ECell[i].vshadow)/3.57 )  
	  + sumAMPA*gAMPAfact + sumfwd;

	// Then inhibitory input 

	sum = 0.0;
	for (int j = 0; j<NIpops; j++)
	  sum += WEI[i][j]*sIpop[j];
	ECell[i].sSynI = sum + ECell[i].cellInh;
	
      }

      /* Then find out if excitatory cell has fired. */
      
      for (int i = 0; i<NE; i++){
	
	ECell[i].sApp = 0.0;
	ECell[i].Integrate(dt,t);
	  
	if ( (ECell[i].cellfire == true ) && ( i < NEpops*NEsub ) ){
	  if ( ( i%(NEsub/2) == 0 ) ) 
	    {
	      Efiredout << i << " " << setprecision(8) << ECell[i].lftime;
	      Efiredout << " " << ECell[i].lftime-ECell[i].llftime; 
	      Efiredout << endl;
	    }
	  Totalnumfired ++;
	}
	
      }
      
      /* Now find inputs to inhibitory cells. */

      for (int i=0; i<NI; i++){

	// Excitatory input first

	double sumNMDA = 0.0;
	double sumAMPA = 0.0;
	for (int j = 0; j<NEpops; j++){
	  sumNMDA += WIE[i][j]*sEpop[j];
	  sumAMPA += WIE[i][j]*sAMPApop[j];
	}
	double sumfwd = WIforward[i/NIsub]*(sEinpop*gNMDAfactfwd/
	      (1.0 + exp(-62*ICell[i].vshadow)/3.57 )  + sAMPAinpop*gAMPAfactfwd);
	double sumRO = WIRO*(sROpop*gNMDAfactfwd/
	      (1.0 + exp(-62*ICell[i].vshadow)/3.57 )  + sROAMPApop*gAMPAfactfwd);

	ICell[i].sSynE = sumNMDA*gNMDAfact/
	      (1.0 + exp(-62*ICell[i].vshadow)/3.57 ) + sumAMPA *gAMPAfact
	  + sumfwd + sumRO;

	// Then inhibitory input.
	
	sum = 0.0;
	for (int j = 0; j<NIpops; j++)
	  sum += WII[i][j]*sIpop[j];
	ICell[i].sSynI = sum + ICell[i].cellInh;
	
      }
      for ( int i=0; i < NIpops; i++ ){
	sIbin[i] += sIpop[i];
      }

      /* Now update membrane potential and see if 
	 inhibitory cell has fired */

      for (int i = 0; i<NI; i++){
	
	ICell[i].Integrate(dt,t);
	
	if (ICell[i].cellfire == true ){
	  if ( (i%(NIsub/2)) == 0 ) {
	    Ifiredout << i << " " << setprecision(8) << ICell[i].lftime;
	    Ifiredout  << " " << ICell[i].lftime-ICell[i].llftime;
	    Ifiredout << endl;
	  }
	  Totalnumfired ++;
	}
	
      }

      /* Calculate inputs to readout cells */
      
      //Excitatory input only

      double sumNMDA = 0.0;
      double sumAMPA = 0.0;
      for (int j = 0; j<NEpops; j++){
	sumNMDA += WER[j]*sEpop[j];
	sumAMPA += WER[j]*sAMPApop[j];
      }

      Epopreadin += sumNMDA*gNMDAfact/
	      (1.0 + exp(-62*ECell[NE].vshadow)/3.57 ) + sumAMPA *gAMPAfact;

      for (int k=0; k<Nreadout; k++){
	ECell[NE+k].sSynE = sumNMDA*gNMDAfact/
	      (1.0 + exp(-62*ECell[NE+k].vshadow)/3.57 ) + sumAMPA *gAMPAfact;;
	ECell[NE+k].sSynI = EpopInh[NEpops];
      }

      /* Then find out if readout cell has fired. */

      for (int i = NE; i<NE+Nreadout; i++){
	ECell[i].sApp = 0.0; 
	ECell[i].Integrate(dt,t);	  
	if ( ECell[i].cellfire == true ){
	  if ( ( (i-NE)%(Nreadout/2) == 0 ) ) 
	    {
	      Efiredout << i << " " << setprecision(8) << ECell[i].lftime;
	      Efiredout << " " << ECell[i].lftime-ECell[i].llftime; 
	      Efiredout << endl;
	    }
	}     
      }

      cellbin[0] += ECell[1].sSynE;
      cellbin[1] += ECell[1].sSynI;
      cellbin[2] += ECell[1].sExt;
      cellbin[3] += ECell[1].v0;
      cellbin[4] += ECell[1].vnow;
      cellbin[5] += ECell[1].tm;
      cellbin[6] += ECell[1].Nves;
      cellbin[7] += sEpop[0]*gNMDAfact/
	      (1.0 + exp(-62*ECell[1].vshadow)/3.57 );
      cellbin[8] += sAMPApop[0]*gAMPAfact;

      cell2bin[0] += ECell[4*NEsub+1].sSynE;
      cell2bin[1] += ECell[4*NEsub+1].sSynI;
      cell2bin[2] += ECell[4*NEsub+1].sExt;
      cell2bin[3] += ECell[4*NEsub+1].v0;
      cell2bin[4] += ECell[4*NEsub+1].vnow;
      cell2bin[5] += ECell[4*NEsub+1].tm;
      cell2bin[6] += ECell[4*NEsub+1].Nves;
      cell2bin[7] += sEpop[4]*gNMDAfact/
	      (1.0 + exp(-62*ECell[4*NEsub+1].vshadow)/3.57 );
      cell2bin[8] += sAMPApop[4]*gAMPAfact;
      cellcount++;

    } // End of bintime loop
    
    cellout << t;
    for (int i=0;i<9;i++){
      cellbin[i] *=dt/bint;
      cellout  << " " << setw(6) << cellbin[i];
    }
    if ( ECell[1].binrate > 0.0 ){
      cellbin[9] /= ECell[1].binrate;
      cellbin[10] /= ECell[1].binrate;
    }
    else{
      cellbin[9] = 0.0;
      cellbin[10] = 0.0;
    }

    cellout  << " " << setw(6) << cellbin[9];
    cellout  << " " << setw(6) << cellbin[10];
    cellout <<  " " << rApp << endl;
    for (int i=0;i<11;i++){
      cellbin[i] = 0.0;
    }

    cell2out << t;
    for (int i=0;i<9;i++){
      cell2bin[i] *=dt/bint;
      cell2out  << " " << setw(6) << cell2bin[i];
    }
    if ( ECell[4*NEsub+1].binrate > 0.0 ){
      cell2bin[9] /= ECell[4*NEsub+1].binrate;
      cell2bin[10] /= ECell[4*NEsub+1].binrate;
    }
    else{
      cell2bin[9] = 0.0;
      cell2bin[10] = 0.0;
    }

    cell2out  << " " << setw(6) << cell2bin[9];
    cell2out  << " " << setw(6) << cell2bin[10];
    cell2out <<  " " << rApp << endl;
    for (int i=0;i<11;i++){
      cell2bin[i] = 0.0;
    }
  
    double Eratebin[NEbins+1];
    double Einratebin;
    double Eonratebin;
    for (int i=0; i<NEbins; i++){
      Eratebin[i] = 0.0;
      for ( int j=0;j<NE/NEbins;j++){
	Eratebin[i] += ECell[j + i*(NE/NEbins)].binrate;
      }
      Eratebin[i] /= (bintime*float(NE/NEbins));
      sEbin[i] *= dt/bintime;
      Eout[i] << t << " " << Eratebin[i] << " " << sEbin[i] << endl;
      sEbin[i] = 0.0;
    }
    
    /* now input cells */
    Einratebin = 0.0;
    for ( int i =0; i<NEsub; i++)
      Einratebin += EinCell[i].binrate;
    Einratebin /= (bintime*float(NEsub));
    sEinbin *= dt/bintime;
    Einout << t << " " << Einratebin << " " << sEinbin << endl;
    sEinbin = 0.0;

    /* now on cells */
    Eonratebin = 0.0;
    for ( int i =0; i<NEsub; i++)
      Eonratebin += EonCell[i].binrate;
    Eonratebin /= (bintime*float(NEsub));
    sEonbin *= dt/bintime;
    Eonout << t << " " << Eonratebin << " " << sEonbin << endl;
    sEonbin = 0.0;

    /* Add bin updates for readout cells */

    Eratebin[NEbins] = 0.0;
    for ( int j=0;j<Nreadout;j++){
      Eratebin[NEbins] += ECell[j + NE].binrate;
    }
    Eratebin[NEbins] /= (bintime*float(Nreadout));
    Eout[NEbins] << t << " " << Eratebin[NEbins] << " " 
		 << Epopreadin << endl;

    Epopreadin = 0.0;

    double Iratebin[NIbins];
    for (int i=0; i<NIbins; i++){
      Iratebin[i] = 0.0;
      for ( int j=0;j<NI/NIbins;j++){
	Iratebin[i] += ICell[j + i*(NI/NIbins)].binrate;
      }    
      Iratebin[i] /= (bintime*float(NI/NIbins));
      sIbin[i] *= dt/bintime;
      Iout[i] << t << " " << Iratebin[i] << " " << sIbin[i] << endl;
      sIbin[i] = 0.0;
    }

    if ( ( t > cuetime[cue] ) && 
	 ( t < cuetime[cue] + cuelength[cue] ) ){
      for ( int i = 0; i<NEbins+1; i++){
	Eresponse[i][cue] += Eratebin[i];
      }
      for ( int i = 0; i<NIbins; i++){
	Iresponse[i][cue] += Iratebin[i];
      }
      Einresponse[cue] += Einratebin;
      rcount[cue] ++;
    }
        
    for ( int i = 0; i< NEsub; i++){
      EinCell[i].binrate = 0.0;
      EonCell[i].binrate = 0.0;
    }
    for ( int i = 0; i< NE+Nreadout; i++){
      ECell[i].binrate = 0.0;
    }
    for ( int i = 0; i< NI; i++){
      ICell[i].binrate = 0.0;
    }

  } // End of time loop

  if ( Totalnumfired >= maxfired ) cout << " MAX. FIRED REACHED! " << endl;

  for ( int j=0; j < Ncues; j++){
    for ( int i = 0; i<NEbins+1; i++){
      Eresponse[i][j] /= float(rcount[j]);
    }
    for ( int i = 0; i<NIbins; i++){
      Iresponse[i][j] /= float(rcount[j]);
    }
    Einresponse[j] /= float(rcount[j]);
  }

     cout << 1 << " " << cuestrength1 << " " << Einresponse[1]  << 
      " " << Eresponse[NEbins][1] << endl;
     cout << 2 << " " << cuestrength2 << " " << Einresponse[2]  << 
      " " << Eresponse[NEbins][2] << endl;

}




