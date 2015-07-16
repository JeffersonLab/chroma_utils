//$Id: coFitter.cc,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: coFitter.cc,v $
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
//Revision 1.4  2004/05/20 17:55:57  kostas
//fixed bug in stdout printing
//
//Revision 1.3  2004/05/20 16:32:23  kostas
//improved the  interface to the fitters and added xml parsing stuff
//for spectrum calculations
//
//Revision 1.2  2004/02/13 15:37:59  kostas
//Just comments
//

// Implements the covariant fitter class
#include "covfit/fitter.h"

namespace CovFit {

  use namespace std ;

/*! sets the covariance matrix. checks the limits xmin,xmax */
void coFitter::SetCovarMat(CovarMat& CM){

  ndata=0;
  for(int i(0);i<CM.Npoints();i++){
    if( CM.Xdata(i)<xlow || CM.Xdata(i)>xhigh )continue;/*reject datapoint*/
    xdata = (double *)realloc(xdata,(ndata+1)*sizeof(double));
    ydata = (double *)realloc(ydata,(ndata+1)*sizeof(double));
    xdata[ndata] = CM.Xdata(i);   ydata[ndata] = CM.Ydata(i);
    ndata++;
  }
  ndof = ndata-npar+fixpar;
  
  nblocks=CM.Nblocks();
  printf("\n\t using %d data points\n",ndata);
  if( nblocks <= ndata ){
    fprintf(stderr,"Not enough blocks in covariance matrix\n");
    exit(0);
  }

  /* initialize matrix */
  covar = (double*)malloc(ndata*ndata*sizeof(double));
  covarinv = (double*)malloc(ndata*ndata*sizeof(double));
  for(int i=0;i<ndata;i++){
    for(int j=0;j<ndata;j++)covar[i*ndata+j] = 0.0;
    covar[i*ndata+i] = -1.0;
    /*diagonal can not be < 0, use to make sure it is specified */
  }
  for(int i(0);i<CM.Npoints();i++)
    for(int j(0);j<CM.Npoints();j++){
      /* find indices of this element */
      int ii,jj ;
      for(ii=0;ii<ndata;ii++)if(CM.Xdata(i) == xdata[ii])break;
      for(jj=0;jj<ndata;jj++)if(CM.Xdata(j) == xdata[jj])break;
      if(ii < ndata && jj < ndata)covar[ii*ndata+jj] = CM[i*CM.Npoints()+j];
      /****if(i != j )covar[i*ndata+j] = 0.0;	 temporary test */
    }
  int flag(0);
  /* check to make sure diagonal elements were specified */
  for(int i=0;i<ndata;i++){
    if(covar[i*ndata+i] < 0.0){
      fprintf(stderr,
	      "Diagonal covariance missing or negative: x = %e\n",xdata[i]);
      flag=1;
    }
  }
  if(flag)exit(1);

  /*invert covariance matrix */
  matinv(covar,covarinv,ndata);
  
  sigma = (double *) realloc(sigma,ndata*sizeof(double));
  
  /* set up the errors on the data points */
  for(int i(0);i<ndata;i++)
    sigma[i]=sqrt(covar[i*ndata+i]) ;

  /* allocate work space */
  wdata1 =  new double[ndata] ;
  wdata2 =  new double[ndata] ;
  wdata3 =  new double[ndata] ;
  wdata4 =  new double[ndata] ;
  wd_alloc = true ;
}

/* dump a matrix in parameter space, skip fixed parameters */
void coFitter::dumpmat(const char *label,const char *format,double *matrix)
{
  int i,j,k;
  printf("\n%s\n",label);
  for(i=0;i<npar;i++){
    for(k=0;k<fixpar;k++)if(parlist[k]==i)break;if(k!=fixpar)continue;
    for(j=0;j<npar;j++){
      for(k=0;k<fixpar;k++)if(parlist[k]==j)break;if(k!=fixpar)continue;
      printf(format,matrix[i*npar+j]);printf(" \t");
    }
    printf("\n");
  }
}


double coFitter::phi(double *b)
{
  double sum(0.0);

  for(int i(0);i<ndata;i++) wdata1[i] = F->f(xdata[i],b) - ydata[i] ;

  for(int i=0;i<ndata;i++)
    for(int j=0;j<ndata;j++)
      sum += wdata1[i]*covarinv[i*ndata+j] * wdata1[j] ;

  return(sum);
}

void coFitter::dphi(double *b, double *grad)
{
  double deriv;
  int i,j,k;

  for(i=0;i<ndata;i++) wdata1[i] = F->f(xdata[i],b) - ydata[i];
  /* compute difference outside double loop */
  for(j=0;j<npar;j++) {
    grad[j]=0.0;
    for(k=0;k<fixpar;k++)if(j==parlist[k])goto bottom;
    for(i=0;i<ndata;i++){
      deriv = F->df(xdata[i],b,j);
      for(k=0;k<ndata;k++){
	grad[j] += 2.0 * deriv * covarinv[i*ndata+k] * wdata1[k];
      }
    }
  bottom:	
    ;
  }

}

void coFitter::ddphi(double *b, double **a)
{
  int i,j,k,l;
  
  for(i=0;i<ndata;i++)wdata4[i] = F->f(xdata[i],b) - ydata[i];
  for(j=0;j<npar;j++){
    for(i=0;i<ndata;i++)wdata1[i] = F->df(xdata[i],b,j);
    for(k=0;k<=j;k++){
      a[j][k]=0.0;
      for(l=0;l<fixpar;l++)if(j==parlist[l] || k ==parlist[l] ){
	if(j==k)a[j][k]=1.0;
	goto bottom;
      }
      for(i=0;i<ndata;i++){
	wdata2[i] = F->df(xdata[i],b,k);
	wdata3[i] = F->ddf(xdata[i],b,j,k);
      }
      for(i=0;i<ndata;i++){
	for(l=0;l<ndata;l++){
	  a[j][k] += 2.0*(wdata1[i]*covarinv[i*ndata+l]*wdata2[l] + 
			  wdata4[i] * covarinv[i*ndata+l] * wdata3[l] );
	}
      }
    bottom:     ;
      a[k][j] = a[j][k] ;
    }
  }
}



/*! 
  correlated least squares fit (designed for QCD hadron propagators).	
									
 Using average function and covariance matrix			
 it minimize quadratic form:						
  (y[i]-f(x[i])) * (inverse of covariance matrix[ij]) * (y[j]-f(x[j]))	
	
*/
void coFitter::Fit(){
  
  singular_Hessian = false ;

  double *delpar   = new double[npar*npar];
  double *wparmat1 = new double[npar*npar];
  double *wparmat2 = new double[npar*npar];

  double *wdata1 = new double[ndata*ndata];
  double *wdata2 = new double[ndata*ndata];

  double **a=(double**)malloc(npar*sizeof(double*));
  for(int i(0);i<npar;i++)a[i]=(double*)malloc(npar*sizeof(double));


  if(display){
    printf("data to be fitted: x,y,error\n");
    for(int i(0);i<ndata;i++)
      printf("%e\t%e\t+/- %e\n",xdata[i],ydata[i],sigma[i]);
  }

  chi_square = min();
  
  if(display){
    /*	print differences		*/
    printf("\n\n\n");
    printf("x\t\tobserved\tpredicted\tdifference/naive_sigma\n");
    for(int i(0);i<ndata;i++){
      double y(F->f(xdata[i],par));	
      double res(ydata[i]-y);
      printf("%6e\t%6e\t",xdata[i],ydata[i]);
      printf("%6e\t%f\n",y,res/sigma[i]) ;
    }
  }
  
  /* make error estimates on parameters */
  
  /* Note the correction factor of sqrt( nblocks/(nblocks-ndata))
     This is analogous to the sqrt(n/(n-1)) in computing the error on
     the mean, which comes from the difference between the sample
     mean and the real mean.  For the fitting program, "n-1"
     is replaced by "n-ndata".  Actually, I should confess that
     I worked this out to lowest order in 1/nblocks, and since
     it blows up at the right point, nblocks=ndata, I use it for
     all nblocks.  DT,  7/22/92 */
  
  /* first look at derivatives of chi-squared */
  /* essentially ask how they vary as chi-squared goes to one more
     than its minimum value */
  printf("\nERROR ANALYSIS FROM CHI-SQUARED CONTOURS\n");
  ddphi(par,a);
  /* second derivatives of chi-square. Unit matrix for fixed params */
  for(int i=0;i<npar;i++)for(int j=0;j<npar;j++)wparmat1[i*npar+j]=a[i][j];
  matinv(wparmat1,delpar,npar);
  dumpmat("Hessian matrix","%e",wparmat1);
  dumpmat("Inverse Hessian matrix","%e",delpar);
  for(int i=0;i<npar;i++){
    for(int j=0;j<npar;j++)
      wparmat2[i*npar+j]=
	delpar[i*npar+j]/sqrt(delpar[i*npar+i]*delpar[j*npar+j]);
  }
  dumpmat("Scaled inverse Hessian matrix","%.3f",wparmat2);

  /* now look at dependence of fit parameters on data points */
  printf("\nERROR ANALYSIS FROM FULL EXPRESSION\n");
  /* wparmat1 = df/dpar * Cinv * df/dpar */

  for(int i=0;i<npar;i++)
    for(int j=0;j<npar;j++){
      for(int l=0;l<fixpar;l++)
	if(i==parlist[l] || j ==parlist[l] ){
	  if(i==j)
	    wparmat1[i*npar+j] = 1.0; 
	  else 
	    wparmat1[i*npar+j]=0.0;
	  goto bottom2;
	}
      wparmat1[i*npar+j] = 0.0;
      for(int k=0;k<ndata;k++){
	wdata1[k] = F->df(xdata[k],par,i);
	wdata2[k] = F->df(xdata[k],par,j);
      }
      for(int k=0;k<ndata;k++)for(int l=0;l<ndata;l++)
	wparmat1[i*npar+j] += wdata1[k]*covarinv[k*ndata+l]*wdata2[l];
    bottom2:    ;
    }
  /* wparmat2 = wparmat1 * delpar (delpar = inverse Hessian) */
  for(int i=0;i<npar;i++)
    for(int j=0;j<npar;j++)
      {
	wparmat2[i*npar+j]=0.0 ;
	for(int k(0);k<npar;k++)
	  wparmat2[i*npar+j] += wparmat1[i*npar+k]*delpar[k*npar+j];
      }
  /* wparmat1 = delpar * wparmat1 (delpar = inverse Hessian) */
  for(int i=0;i<npar;i++)
    for(int j=0;j<npar;j++)
      {
	wparmat1[i*npar+j]=0.0;
	for(int k(0);k<npar;k++)
	  wparmat1[i*npar+j] += 4.0*delpar[i*npar+k]*wparmat2[k*npar+j];
      }
  dumpmat("Parameter variance matrix","%e",wparmat1);
  for(int i=0;i<npar;i++){
    for(int j=0;j<npar;j++)
      wparmat2[i*npar+j] = 
	wparmat1[i*npar+j]/sqrt(wparmat1[i*npar+i]*wparmat1[j*npar+j]);
  }
  dumpmat("Scaled parameter variance matrix","%.3f",wparmat2);
  
  printf("\nFINAL PARAMETERS, CHI_SQUARE ERRORS, FULL ERRORS\n");
  printf("chi-square with %d degrees of freedom is %e\n",ndof,chi_square);
  confidence = gammq(0.5*ndof,0.5*chi_square) ;
  printf("confidence of fit is\t%e\n\n", confidence);
  for(int i(0);i<npar;i++){
    epar[i]=0.0 ;
    int j;
    for(j=0;j<fixpar;j++)if(parlist[j]==i)break;
    if(j!=fixpar) printf("par[%d] = %e  +-  (fixed)\n",i,par[i]);
    else{
      epar[i] = sqrt( (double)nblocks/(nblocks-ndata) * wparmat1[i*npar+i]) ;
      double err = sqrt((double)nblocks/(nblocks-ndata)*2.0*delpar[i*npar+i]) ;
      printf("par[%d] = %e  +-  %e \t%e\n",i,par[i],err,epar[i]);
      if(isnan(err)||isnan(epar[i]))
	fit_failed=true ;
      else
	fit_failed=false ;
    }
  }

  delete [] delpar   ;
  delete [] wparmat1 ;
  delete [] wparmat2 ;

  delete [] wdata1 ;
  delete [] wdata2 ;

  for(int i=0;i<npar;i++) free(a[i]) ;
  free(a);
}

} // namespace CovFit
