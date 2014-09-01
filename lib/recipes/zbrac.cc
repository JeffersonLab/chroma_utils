#include <cmath>
#include "nr.h"
using namespace std;

bool NR::zbrac(const ScalFunc &func, DP &x1, DP &x2)
{
  //bool NR::zbrac(DP func(const DP), DP &x1, DP &x2)
  //{
  const int NTRY=10;     //50; //jjd made this give up much quicker
  const DP FACTOR=1.6;
  int j;
  DP f1,f2;
  
  if (x1 == x2) nrerror("Bad initial range in zbrac");
  f1=func(x1);
  f2=func(x2);
  for (j=0;j<NTRY;j++) {

    // std::cout << " zbrac j= " <<j << ", f(x1= " << x1 <<")=" << f1 
    //	      << ", f(x2= " << x2 <<")=" << f2 << endl;


    if (f1*f2 < 0.0) return true;
    if (fabs(f1) < fabs(f2))
      f1=func(x1 += FACTOR*(x1-x2));
    else
      f2=func(x2 += FACTOR*(x2-x1));
  }
  return false;
}

//and to use the original form of the call:
bool NR::zbrac(DP func(const DP), DP &x1, DP &x2)
{
  PtrScalFunc p(func);
  return NR::zbrac( p, x1, x2 );
};


//DJW: Version that limits x1 and x2 to positive values only
bool NR::zbrac_pos(const ScalFunc &func, DP &x1, DP &x2)
{
  //bool NR::zbrac(DP func(const DP), DP &x1, DP &x2)
  //{
  const int NTRY=10;     //50; //jjd made this give up much quicker
  const DP FACTOR=1.6;
  int j;
  DP f1,f2;
  
  if (x1 == x2) nrerror("Bad initial range in zbrac");
  f1=func(x1);
  f2=func(x2);
  for (j=0;j<NTRY;j++) {

    // std::cout << " zbrac j= " <<j << ", f(x1= " << x1 <<")=" << f1 
    //	      << ", f(x2= " << x2 <<")=" << f2 << endl;

    if (f1*f2 < 0.0) return true;
    if (fabs(f1) < fabs(f2)){
      x1 = x1 + FACTOR*(x1-x2);
      if (x1<0.0){x1=0.0;};
      f1 = func(x1);
    } else {
      x2 = abs(x2 + FACTOR*(x2-x1));
      f2 = func(x2);
    };
  };
  return false;
}

//and to use the original form of the call:
bool NR::zbrac_pos(DP func(const DP), DP &x1, DP &x2)
{
  PtrScalFunc p(func);
  return NR::zbrac_pos( p, x1, x2 );
};


//JJD: Version that limits x1 and x2 to negative values only
bool NR::zbrac_neg(const ScalFunc &func, DP &x1, DP &x2)
{
  const int NTRY=10;     //50; //jjd made this give up much quicker
  const DP FACTOR=1.6;
  int j;
  DP f1,f2;

  if (x1 == x2) nrerror("Bad initial range in zbrac");
  f1=func(x1);
  f2=func(x2);
  for (j=0;j<NTRY;j++) {

    // std::cout << " zbrac j= " <<j << ", f(x1= " << x1 <<")=" << f1 
    //	      << ", f(x2= " << x2 <<")=" << f2 << endl;

    if (f1*f2 < 0.0) return true;
    if (fabs(f1) < fabs(f2)){
      x1 =  x1 + FACTOR*(x1-x2) ;
      f1 = func(x1);
    } 
    else {
      x2 = x2 + FACTOR*(x2-x1);
      if (x2>0.0){x2=0.0;};
      f2 = func(x2);
    };
  };
  return false;
}

//and to use the original form of the call:
bool NR::zbrac_neg(DP func(const DP), DP &x1, DP &x2)
{
  PtrScalFunc p(func);
  return NR::zbrac_neg( p, x1, x2 );
};
