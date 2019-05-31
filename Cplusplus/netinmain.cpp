/* main code for running the discriminator model 
   should be compiled with discrimnet.cpp and structure.cpp */

#include "MersenneTwister.h"
#include "dnet.h"

MTRand rand1;
MTRand rand2;
int seed1 = 1;
int seed2 = 1;

int main(int argc, char **argv){

  int numin;

  if(argc > 1)   //if there IS an argument
    numin = atoi(argv[1]);
  
  cout << " numin " << numin << endl;

  if(numin >= 10000) {                    //code for 2 arguments
    seed1 = atoi(argv[1])/10000;          // first digits
    numin = atoi(argv[1]) - 10000*seed1;  // last 4 digits
  }
  seed2 = seed1+100;
  rand1.seed(seed1);
  rand2.seed(seed2);
  
  cout << " numin " << numin << endl;
  int istrength1 = numin/2;
  cout << " istrength 1 " << istrength1 << endl;
  int iud = numin - 2*istrength1;

  cout << " iud " << iud << endl;

  /* The following lines set f2 in relation to f1*/
  int istrength2 = istrength1;
  if ( iud > 0 ) 
    istrength2 += 8;
  else
    istrength2 -= 8;

  cout << " seed1 " << seed1 << " seed2 " << seed2 << 
    " cue 1 " << istrength1 << " cue 2 " << istrength2 << endl;
  double a = rand1();
  double b = rand2();
  cout << " rand1 " << a << " rand2 " << b << endl;
  NET net1;

  cout << " Main" << endl;

  /* Input single cell and simulation parameters */
  net1.Datin(istrength1,istrength2);
  cout << " done Datin " << endl;

  /* Input cpnnection strengths */
  net1.WDatin();

  /* Initialize all cells */
  net1.Init();

  /* Run the trial*/
  net1.GatherInputs();

}

string createString(string prefix, int num, int digits){

  string filenumber = "0123456789";
  string RUNNING;
  int thousands;
  int hundreds;
  int tens;
  int units;
  string U = "";
  string T = "";
  string H = "";
  string TH = "";

  switch(digits){
  case 4:
    thousands = num/1000;
    TH = filenumber.at(thousands);
  case 3:
    cout << " 3 " << endl;
    hundreds = (num%1000)/100;
    H = filenumber.at(hundreds);
  case 2:
    tens = (num%100)/10;
    T = filenumber.at(tens);
  case 1:
    units = num%10;    
    U = filenumber.at(units);
    break;
  default:
    TH = "0";
    H = "0";
    T="0";
    U = "0";
    break;
  }
  cout << prefix+TH+H+T+U << endl;
  return prefix+TH+H+T+U;
}

/* Ran_Gaussian requires 2 random number calls, which it 
   uses to produce a Gaussian random variable with mean of 
   zero and s.d. of 1 */
double Ran_Gaussian(){

  static int jj=0;
  static double kk;
  double root_factor;
  double sum_square;
  double ran_num_1;
  double ran_num_2;
  
  if (jj == 0) {
    do {
      ran_num_1 = 2.0*rand2()-1.0;
      ran_num_2 = 2.0*rand2()-1.0;
      sum_square = ran_num_1*ran_num_1+ran_num_2*ran_num_2;
    } while (sum_square >= 1.0 || sum_square == 0);
    
    root_factor = sqrt(-2.0*log(sum_square)/sum_square);
    kk=ran_num_1*root_factor;
    jj=1;
    return ran_num_2*root_factor;
  } 
  else {
    jj = 0;
    return kk;
  }
  
}
