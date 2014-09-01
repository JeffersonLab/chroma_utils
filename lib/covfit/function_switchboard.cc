
#include "covfit/fitters/fitfunctions.h"

using namespace std ;

namespace CovFit {

  Function* CreateFunction(std::string& name, int N){
    // The N should be the number of parametes   
    Function *f ;
    
    if(name == "AltExponentials"){
      f = new AltExponentials(N/4);
    }
    else if(name == "Exponentials"){
      f = new Exponentials(N/2);
    }
    else if(name == "Polynomium"){
      f = new Polynomium(N);
    }
    else if(name == "BoundMom"){
       f = new BoundMom();
    }
    else{
      string foo = "Unknown function: "+name;
      throw  foo ;
    }
  
    return f ;
  }

  Function* CreateFunction(std::string& name, double p, int N){
    // The N should be the number of parametes   
  Function *f ;

  if(name == "AltExponentialsPeriodic"){
    f = new AltExponentialsPeriodic(p,N/4);
  }
  if(name == "ExponentialsAPeriodicBar"){
    f = new ExponentialsAPeriodicBar(p,N/4);
  }
  else if(name == "ExponentialsAPeriodic"){
    f = new ExponentialsAPeriodic(p,N/2);
  }
  else if(name == "ExponentialsPeriodic"){
    f = new ExponentialsPeriodic(p,N/2);
  }
  else if(name == "ExponentialsPeriodicSubtr"){
    f = new ExponentialsPeriodicSubtr(p,N/2);
  }
  else if(name == "AltExponentialsPeriodicOneMass"){
    f = new AltExponentialsPeriodicOneMass(p,N/3);
  }
  else{
    string foo = "Unknown periodic function: "+name;
    throw  foo ;
  }
  
  return f ;
}




}// namespace CovFit

