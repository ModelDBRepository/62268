/* structure.cpp uses data input from connections_in.dat to 
   generate the full set of connection strengths between 
   subpopulations of neurons */

#include "MersenneTwister.h"
#include "dnet.h"

extern MTRand rand1;
extern MTRand rand2;

void NET::WDatin(){
  string dummy;
  int NIpops = NIpops0;

  ifstream connectin;
  connectin.open("connections_in.dat");
  double state;
  double W_amp;
  double W_diag;
  double W_sigma;
  double WEE_amp[NEpops];
  double WEE_sigma[NEpops];
  double Eaxon_sig;
  double Iaxon_sig;

  /* The parameters sigma_EE etc. are standard deviations 
     for connection weights (set to zero as default) */
  connectin >> sigma_EE >> dummy;
  connectin >> sigma_EI >> dummy;
  connectin >> sigma_IE >> dummy;
  connectin >> sigma_II >> dummy;
  /* Eaxon_sig and Iaxon_sig can be set to give standard deviation 
     for heterogeneity in synaptic weights (zero by default) */
  connectin >> Eaxon_sig >> dummy;
  connectin >> Iaxon_sig >> dummy;
  logOut << " Eaxon_sig " << Eaxon_sig << endl;
  logOut << " Iaxon_sig " << Iaxon_sig << endl;

  for ( int i = 0; i < NE ; i++){
    Eaxonfact[i] = 1.0 + Eaxon_sig*Ran_Gaussian();
    if ( Eaxonfact[i] < 0.0 ) Eaxonfact[i] = 0.0;
    if ( Eaxonfact[i] > 2.0 ) Eaxonfact[i] = 2.0;
  }
  for ( int i = 0; i < NE ; i++){
    Eaxonfact2[i] = 1.0 + Eaxon_sig*Ran_Gaussian();
    if ( Eaxonfact2[i] < 0.0 ) Eaxonfact2[i] = 0.0;
    if ( Eaxonfact2[i] > 2.0 ) Eaxonfact2[i] = 2.0;
  }
  for ( int i = 0; i < NI ; i++){
    Iaxonfact[i] = 1.0 + Iaxon_sig*Ran_Gaussian();
    if ( Iaxonfact[i] < 0.0 ) Iaxonfact[i] = 0.0;
    if ( Iaxonfact[i] > 2.0 ) Iaxonfact[i] = 2.0;
  }
  for ( int i = 0; i < NI ; i++){
    Iaxonfact2[i] = 1.0 + Iaxon_sig*Ran_Gaussian();
    if ( Iaxonfact2[i] < 0.0 ) Iaxonfact2[i] = 0.0;
    if ( Iaxonfact2[i] > 2.0 ) Iaxonfact2[i] = 2.0;
  }
    
  /*   Now the main connection strengths */

  double WEEij;
  connectin >> WEEij >> dummy;

  for ( int j = 0; j < NEpops; j++){
    connectin >> state >> dummy;
    WEE_amp[j] = state*(1.0+sigma_EE*Ran_Gaussian());
  }

  for ( int j = 0; j < NEpops; j++){
    connectin >> state >> dummy;
    WEE_sigma[j] = state;
  }

  /* W_asym is an asymmetry factor, set to 1 for symmetric 
     connections between numbered memory subpopulations. 
     W_asym > 1 means slower decay, hence stronger connections 
     from high threshold to low threshold.  */

  connectin >> state >> dummy;
  double W_asym = state;

  ofstream WEEout;
  WEEout.open("WEEout.dat");
  for ( int i = 0; i < NEpops ; i++){
    float frac_i = float(i)/float(NEpops-1);
    for ( int j = 0; j < i ; j++){
      float frac_j = float(j)/float(NEpops-1);
      double d_ij = -frac_j + frac_i;
      double sig_ij = WEE_sigma[i];

      /* Reduce s.d. by a factor W_asym and weaken (strengthen) 
	 connections if j < i when W_asym is greater (less) than one */
      sig_ij /= W_asym;
      /* WEEpop[i][j] is strength from j to i */
      WEEpop[i][j] = WEEij*exp(-d_ij/sig_ij)/float(NEpops);

      if ( WEEpop[i][j] < 0.0 ) 
	WEEpop[i][j] = 0.0;

    }
    WEEpop[i][i] = WEE_amp[i];
    logOut << " WEE " << i+1 << "," << i+1 << "  " << WEEpop[i][i] << endl;
    WEEout << i+1 << " " << i+1 << " " << WEEpop[i][i] << endl;

    for ( int j = i+1; j < NEpops ; j++){
      float frac_j = float(j)/float(NEpops-1);
      double d_ij = frac_j - frac_i;
      double sig_ij = WEE_sigma[i];

      /* Increase s.d. by a factor W_asym and strengthen (weaken) 
	 connections if j > i when W_asym is greater (less) than one */
      sig_ij *= W_asym;

      /* WEEpop[i][j] is strength from j to i */
      WEEpop[i][j] = WEEij*exp(-d_ij/sig_ij)/float(NEpops);

      if ( WEEpop[i][j] < 0.0 ) 
	WEEpop[i][j] = 0.0;

      //      logOut << " WEE " << i+1 << "," << j+1 << "  " << WEEpop[i][j] << endl;
      //      WEEout << i+1 << " " << j+1 << " " << WEEpop[i][j] << endl;

    }
    
  }

  connectin >> state >> dummy;
  W_amp = state*(1.0+sigma_EI*Ran_Gaussian());
  connectin >> state >> dummy;
  W_diag = state;
  connectin >> state >> dummy;
  W_asym = state;
  connectin >> state >> dummy;
  W_sigma = state;
  logOut << " WEI_amp " << W_amp << "  WEI_diag " << W_diag << endl;
  logOut << " WEI_asym " << W_asym << "  WEI_sigma " << W_sigma << endl;

  for ( int i = 0; i < NEpops ; i++){
    float frac_i = float(i)/float(NEpops-1);
    for ( int j = 0; j < NIpops ; j++){
      float frac_j;
      if ( NIpops > 1 ) 
	frac_j = float(j)/float(NIpops-1);
      else 
	frac_j = 0.5;

      double d_ij = frac_j - frac_i;
      if ( d_ij >= 0.0 ) 
	//	WEIpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*d_ij)
	WEIpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*double(2*i-NEpops)/double(NEpops))
	    * exp(-d_ij/W_sigma) / float(NIpops);
      else
	//	WEIpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*d_ij)	
	WEIpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*double(2*i-NEpops)/double(NEpops))
	  * exp(d_ij/W_sigma) / float(NIpops);

      if ( WEIpop[i][j] < 0.0 ) 
	WEIpop[i][j] = 0.0;
      //      logOut << " WEI " << i+1 << "," << j+1 << "  " << WEIpop[i][j] << endl;
    }
  }

  connectin >> state >> dummy;
  W_amp = state*(1.0+sigma_IE*Ran_Gaussian());
  connectin >> state >> dummy;
  W_diag = state;
  connectin >> state >> dummy;
  W_asym = state;
  connectin >> state >> dummy;
  W_sigma = state;
  logOut << " WIE_amp " << W_amp << "  WIE_diag " << W_diag << endl;
  logOut << " WIE_asym " << W_asym << "  WIE_sigma " << W_sigma << endl;

  for ( int i = 0; i < NIpops ; i++){
    float frac_i;
    if ( NIpops > 1 ) 
      frac_i = float(i)/float(NIpops-1);
    else 
      frac_i = 0.5;

    for ( int j = 0; j < NEpops ; j++){
      float frac_j = float(j)/float(NEpops-1);
      double d_ij = frac_j - frac_i;
      if ( d_ij >= 0.0 ) 
	//	WIEpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*d_ij)
	WIEpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*double(2*j-NEpops)/double(NEpops) )
	    * exp(-d_ij/W_sigma) / float(NEpops);
      else
	//	WIEpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*d_ij)	
	WIEpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*double(2*j-NEpops)/double(NEpops) )
	  * exp(d_ij/W_sigma) / float(NEpops);

      if ( WIEpop[i][j] < 0.0 ) 
	WIEpop[i][j] = 0.0;
      //      logOut << " WIE " << i+1 << "," << j+1 << "  " << WIEpop[i][j] << endl;
    }
  }

  connectin >> state >> dummy;
  W_amp = state*(1.0+sigma_II*Ran_Gaussian());
  connectin >> state >> dummy;
  W_diag = state;
  connectin >> state >> dummy;
  W_asym = state;
  connectin >> state >> dummy;
  W_sigma = state;
  logOut << " WII_amp " << W_amp << "  WII_diag " << W_diag << endl;
  logOut << " WII_asym " << W_asym << "  WII_sigma " << W_sigma << endl;

  for ( int i = 0; i < NIpops ; i++){
    float frac_i;
    if ( NIpops > 1 ) 
      frac_i = float(i)/float(NIpops-1);
    else 
      frac_i = 0.5;

    for ( int j = 0; j < NIpops ; j++){
      float frac_j;
      if ( NIpops > 1 ) 
	frac_j = float(j)/float(NIpops-1);
      else 
	frac_j = 0.5;

      double d_ij = frac_j - frac_i;
      if ( d_ij >= 0.0 ) 
	WIIpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*d_ij)
	    * exp(-d_ij/W_sigma) / float(NIpops);
      else
	WIIpop[i][j] = W_amp*(1.0 + W_diag*frac_j + W_asym*d_ij)	
	  * exp(d_ij/W_sigma) / float(NIpops);

      if ( WIIpop[i][j] < 0.0 ) 
	WIIpop[i][j] = 0.0;
      //      logOut << " WII " << i+1 << "," << j+1 << "  " << WIIpop[i][j] << endl;
    }
  }

  for (int j = 0; j < NEpops; j++) {
    for ( int i = 0; i < NEpops ; i++){
      for ( int k = 0; k< NEsub; k++ ){
	double fact = 1.0 + sigma_EE*Ran_Gaussian();
	if ( fact < 0.0 ) fact = 0.0;
	if ( fact > 2.0 ) fact = 2.0;
	WEE[i*NEsub + k][j] = WEEpop[i][j]*fact;
	//	WEEout << i*NEsub + k << " " << j << " " <<  WEE[i*NEsub + k][j] << endl;
      }
      WEEout << i << " " << j << " " <<  WEE[i*NEsub][j] << endl;
    }
  }


  for (int j = 0; j < NIpops; j++) {
    for ( int i = 0; i < NEpops ; i++){
      for ( int k = 0; k< NEsub; k++ ){
	double fact = 1.0 + sigma_EI*Ran_Gaussian();
	if ( fact < 0.0 ) fact = 0.0;
	if ( fact > 2,0 ) fact = 2.0;
	WEI[i*NEsub + k][j] = WEIpop[i][j]*fact;
      }
    }
  }


  for (int j = 0; j < NEpops; j++) {
    for ( int i = 0; i < NIpops ; i++){
      for ( int k = 0; k< NIsub; k++ ){
	double fact = 1.0 + sigma_IE*Ran_Gaussian();
	if ( fact < 0.0 ) fact = 0.0;
	if ( fact > 2,0 ) fact = 2.0;
	WIE[i*NIsub + k][j] = WIEpop[i][j]*fact;
      }
    }
  }

  for (int j = 0; j < NIpops; j++) {
    for ( int i = 0; i < NIpops ; i++){
      for ( int k = 0; k< NIsub; k++ ){
	double fact = 1.0 + sigma_II*Ran_Gaussian();
	if ( fact < 0.0 ) fact = 0.0;
	if ( fact > 2,0 ) fact = 2.0;
	WII[i*NIsub + k][j] = WIIpop[i][j]*fact;
      }
    }
  }

    
  /* Now add in connections to readout group */

  for (int i = 0; i < NEpops; i++){
    connectin >> state >> dummy;
    WER[i] = state / float(NEpops);;
    logOut << " WER " << i+1  << " " << WER[i] << endl;
  }

  /* Now add in input-integrator connections */
  double WEf; // Average feedforward excitation
  double WEfv; // If non-zero gives a range of excitation (= 0)  
  double WEb; // Excitatory feedback (= 0)
  double WIf; // Feedforward inhibition (= 0)
  double WIb; // Feedback inhibition strength from interneurons

  connectin >> WEf >> dummy;
  connectin >> WEfv >> dummy;
  connectin >> WEb >> dummy;
  connectin >> WIf >> dummy;
  connectin >> WIb >> dummy;
  connectin >> WEEin >> dummy;

  logOut << " WEf " << WEf << endl;
  logOut << " WEfvary " << WEfv << endl;
  logOut << " WEb " << WEb << endl;
  logOut << " WIf " << WIf << endl;
  logOut << " WIb " << WIb << endl;
  logOut << " WEEin " << WEEin << endl;

  for (int i = 0; i< NEpops; i++ ){
    WEforward[i] = WEf + WEfv*(double(NEpops-1)/2.0-i)/double(NEpops-1);
    WEback[i] = WEb/double(NEpops);
  }
  for (int i = 0; i< NIpops; i++ ){
    WIforward[i] = WIf;
    WIback[i] = WIb/double(NIpops);
  }

  /* Now add in cross-connections to oppositely tuned groups.
     Note that in the discrete integrator used in the basic 
     discrimination model, there is no negatively tuned set 
     of memory neurons, so all cross-connections are zero */ 
  
  connectin >> WEEcross >> dummy;
  connectin >> WEIcross >> dummy;
  connectin >> WIEcross >> dummy;
  connectin >> WIIcross >> dummy;
  logOut << " WEEcross " << WEEcross << endl;
  logOut << " WEIcross " << WEIcross << endl;
  logOut << " WIEcross " << WIEcross << endl;
  logOut << " WIIcross " << WIIcross << endl;


  /* WIRO is connection strength from Readout to Inhibitory neurons*/
  connectin >> WIRO >> dummy;
  logOut << " WIRO " << WIRO << endl;

  /* WEon is strength of excitation from ON cells to Comparison cells.
     WEEon is recurrent excitation within ON cells to make the
     ON populations bistable */  
  connectin >> WEon >> dummy;
  connectin >> WEEon >> dummy;
  logOut << " WEon " << WEon << endl;
  logOut << " WEEon " << WEEon << endl;


}






