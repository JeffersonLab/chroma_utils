//$Id: Fitter.cc,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: Fitter.cc,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.2  2008/11/30 02:49:34  edwards
//Put in some judicious const char* in place of char* to make the compiler
//happy.
//
//Revision 1.1  2005/04/27 19:53:52  edwards
//Completely overhauled include and lib structure. Put all
//of Kostas' fit stuff into their own covfit subdir in include
//and lib. All in its own namespace. Similarly, all of
//ensem is in it's own dirs and namespaces. Unified use of
//XML Array as base array type. I've not touched TTT which
//uses someother classes.
//The new convention should be is people put there little mini
//packages in their own dirs and namespaces. One can compose
//a bigger namespace with namespace composition.
//
//Revision 1.4  2004/05/20 16:32:23  kostas
//improved the  interface to the fitters and added xml parsing stuff
//for spectrum calculations
//
//Revision 1.3  2004/02/26 20:39:42  kostas
//added functionality to the fitters so that they can be called
//fro c++ programs
//
//Revision 1.2  2004/02/13 15:37:58  kostas
//Just comments
//

//
// This is just a c++ wrap around the c fitting codes provided by
// Doug Toussaint. A description of the algorithms used can be
// found in Doug's TASI lectures.
//
#include "covfit/fitter.h"

namespace CovFit {

void Fitter::ReadData(const char* dataf)
{
  FILE  *fp ;    
  fp = fopen(dataf,"r") ;
  if(fp == NULL){
    fprintf(stderr,"Cannot open data file %s\n",dataf);
    exit(1);
  }
  
  double xx,yy,ww;
  while(fscanf(fp,"%le %le  %le",&xx,&yy,&ww) == 3 ){
    if(xx < xlow || xx > xhigh ) continue;
    xdata  = (double *)realloc(xdata, (ndata+1)*sizeof(double));
    ydata  = (double *)realloc(ydata, (ndata+1)*sizeof(double));
    sigma  = (double *)realloc(sigma, (ndata+1)*sizeof(double));
    weight = (double *)realloc(weight,(ndata+1)*sizeof(double));
    xdata[ndata] = xx;   ydata[ndata] = yy; sigma[ndata] = ww; 
    weight[ndata] = 1.0/(ww*ww) ;
    ndata++;
  }
  printf("\n\t read in %d data points\n",ndata);
  fclose(fp);
}


void Fitter::setData(int N, double *x, double *y, double *w)
{
  ndata = 0 ;
  for(int i(0);i<N;i++)
    if( x[i] >= xlow && x[i] <= xhigh ) 
      ndata++ ;

  weight = (double *)realloc(weight,ndata*sizeof(double));
  xdata  = (double *)realloc(xdata ,ndata*sizeof(double));
  ydata  = (double *)realloc(ydata ,ndata*sizeof(double));
  sigma  = (double *)realloc(sigma ,ndata*sizeof(double));

  int d(0);
  for(int i(0);i<N;i++)
    {
      if( x[i] >= xlow && x[i] <= xhigh ) 
	{
	  xdata [d] = x[i];
	  ydata [d] = y[i];
	  sigma [d] = w[i];
	  weight[d] = 1.0/(w[i]*w[i]);
	  d++ ;
	}
    }
  ndof = ndata-npar+fixpar;
}

void Fitter::setData(Array<double>& x, Array<double>& y, Array<double>& w)
{
  int N(x.size());

  ndata = 0 ;
  for(int i(0);i<N;i++)
    if( x[i] >= xlow && x[i] <= xhigh ) 
      ndata++ ;
  
  weight = (double *)realloc(weight,ndata*sizeof(double));
  xdata  = (double *)realloc(xdata ,ndata*sizeof(double));
  ydata  = (double *)realloc(ydata ,ndata*sizeof(double));
  sigma  = (double *)realloc(sigma ,ndata*sizeof(double));

  int d(0);
  for(int i(0);i<N;i++)
    {
      if( x[i] >= xlow && x[i] <= xhigh ) 
	{
	  xdata [d] = x[i];
	  ydata [d] = y[i];
	  sigma [d] = w[i];
	  weight[d] = 1.0/(w[i]*w[i]);
	  d++ ;
	}
    }
  ndof = ndata-npar+fixpar;
}


double Fitter::phi(double *b)
{
	double sum,diff;
	int i;
	sum=0;
	for(i=0,sum=0 ;i<ndata;i++){
		diff= F->f( xdata[i] , b ) - ydata[i];
		sum += weight[i] * diff * diff ;
	}
	return(sum);
}

void Fitter::dphi(double *b, double *grad)
{
	double diff;
	int i,j,k;
	for(j=0;j<npar;j++) {
		grad[j]=0;
		for(k=0;k<fixpar;k++)if(j==parlist[k])goto bottom;
		for(i=0;i<ndata;i++){
			diff= F->f( xdata[i] , b ) - ydata[i];
			grad[j] += 2*weight[i] * diff * F->df(xdata[i],b,j); 
		}
bottom:	
		;
	}
}

void Fitter::ddphi(double *b, double **a)
{
	int i,j,k,l;
	for(j=0;j<npar;j++)for(k=0;k<=j;k++){
		a[j][k]=0;
		for(l=0;l<fixpar;l++)if(j==parlist[l] || k ==parlist[l] ){
			if(j==k)a[j][k]=1;
			goto bottom;
		}
		for(i=0;i<ndata;i++) a[j][k] += 2*weight[i]*
		    ( F->df(xdata[i],b,j) * F->df(xdata[i],b,k) + 
		    (F->f(xdata[i],b)-ydata[i]) * F->ddf(xdata[i],b,j,k) );
		a[k][j] = a[j][k] ;
bottom:	
		;
	}
}


/*! 
  Does uncorrelated fit
*/
void Fitter::Fit()
{
  int i,j,k ;
  double res ;
  double *x=par;
  double *dx=epar;
  double y ;
  double *yt,*b ;
  double **a ;
  int *p,*q;

  yt = new double[npar] ;
  b = new double[npar] ;
  p = new int[npar];
  q = new int[npar];
  //a=(double**)malloc(npar*sizeof(double*));
  //a[0]=(double*)malloc(npar*npar*sizeof(double));
  //for(i=1;i<npar;i++)a[i]=a[0]+i*npar ;
  a=(double**)malloc(npar*sizeof(double*));
  for(i=0;i<npar;i++)a[i]=(double*)malloc(npar*sizeof(double));

  if(display){
    printf("data to be fitted: x,y,error\n");
    for(i=0;i<ndata;i++)
      printf("%e\t%e\t+/- %e\n",xdata[i],ydata[i],1.0/sqrt(weight[i]));
  }
  
  chi_square = min() ;
  if(display){
    /*	print differences		*/
    printf("\n\n\n");
    printf("x\t\tobserved\tpredicted\tdifference/sigma\n");
    for(i=0;i<ndata;i++){
      y=F->f(xdata[i],x);
      res=ydata[i]-y;
      printf("%6e\t%6e\t",xdata[i],ydata[i]);
      printf("%6e\t%f\n",y,res*sqrt(weight[i]));
    }
  }
  printf("\nchi-square with %d degrees of freedom is %e\n",ndof,chi_square);
  confidence = gammq(0.5*ndof,0.5*chi_square) ;
  printf("confidence of fit is\t%e\n\n",confidence);
  ddphi(x,a);
  for(i=0;i<npar;i++)dx[i]=0;
  if( 0 != factor(a,p,q,npar) ){
    for(i=0;i<ndata;i++){
      for(j=0;j<npar;j++){
	b[j]=0;
	for(k=0;k<fixpar;k++)if(j==parlist[k])goto bottom;
	b[j] = 2*weight[i] * F->df(xdata[i],x,j); 
      bottom:	
	;
      }
      subst(a,b,yt,p,q,npar);
      for(j=0;j<npar;j++) dx[j] += yt[j]*(yt[j]/weight[i]);
    }
    for(i=0;i<npar;i++) epar[i]=dx[i]=sqrt(dx[i]);
    singular_Hessian = false ;
  }
  else{
    printf("Hessian is singular\tno error estimates\n");
    singular_Hessian = true ;
  }
  printf("\n\nFINAL PARAMETERS:\n\n");
  for(i = 0 ; i < npar ; i++){
    for(j=0;j<fixpar;j++)if(parlist[j]==i)break;
    if(j==fixpar)printf("x[%d]\t=\t%e\t+-\t%e\n",i,x[i],dx[i]);
    else       printf("x[%d]\t=\t%e\t+-\t(fixed)\n",i,x[i]);
  }
  //graph("fit_graph.ax");

  delete [] yt ;
  delete [] b ;
  delete [] p ;
  delete [] q ;
  for(i=0;i<npar;i++) free(a[i]) ;
  free(a);
}

} // namespace CovFit
