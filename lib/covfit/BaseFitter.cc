//$Id: BaseFitter.cc,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: BaseFitter.cc,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.2  2008/11/30 02:49:33  edwards
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

// Implements the Fitter base class

#include "covfit/fitter.h"

namespace CovFit {

  use namespace std; 

/* From Numerical Recipes:  */
/*! To get confidence: gammq( 0.5*dof , 0.5*chisq ) */
double BaseFitter::gammq(double a, double x)
{
  double gamser,gammcf,gln;
  
  if( a==0.0 )return(0.0);
  if (x < 0.0 || a <= 0.0)
    fprintf(stderr,"Invalid arguments in routine GAMMQ %e %e\n", a,x);
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}

#define ITMAX 200
#define EPS 1.0e-10

void BaseFitter::gser(double *gamser,double a,double x,double *gln)
{
  int n;
  double sum,del,ap;
  
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) fprintf(stderr,"x less than 0 in routine GSER");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    fprintf(stderr,"a too large, ITMAX too small in routine GSER");
    return;
  }
}

void BaseFitter::gcf(double *gammcf,double a, double x,double *gln)
{
  int n;
  double gold=0.0,g,fac=1.0,b1=1.0;
  double b0=0.0,anf,ana,an,a1,a0=1.0;
  
  *gln=gammln(a);
  a1=x;
  for (n=1;n<=ITMAX;n++) {
    an=(double) n;
    ana=an-a;
    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;
    anf=an*fac;
    a1=x*a0+anf*a1;
    b1=x*b0+anf*b1;
    if (a1) {
      fac=1.0/a1;
      g=b1*fac;
      if (fabs((g-gold)/g) < EPS) {
	*gammcf=exp(-x+a*log(x)-(*gln))*g;
	return;
      }
      gold=g;
    }
  }
  fprintf(stderr,"a too large, ITMAX too small in routine GCF");
}

#undef ITMAX
#undef EPS

double BaseFitter::gammln(double xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

void BaseFitter::graph(char *fitgraph)
{
  /* Routine to draw experimental data and compare to fit	*/
  int i,j;
  double tx,z,xmin=1e30,xmax= -1e30;
  FILE *fp;
  fp=fopen(fitgraph,"w");
  for(i=0;i<ndata;i++){
    z = xdata[i];
    if( z < xmin )xmin = z;
    if( z > xmax )xmax = z;
  }
  
  fprintf(fp,"%s","%");
  fprintf(fp,"chi-square with %d degrees of freedom is %e\n",ndof,chi_square);
  fprintf(fp,"%s","%");
  fprintf(fp,"confidence of fit is\t%e\n\n",confidence);
  for(i=0;i<npar;i++){ 
    if(singular_Hessian){
      fprintf(fp,"%s","%");
      fprintf(fp,"Hessian is singular\tno error estimates\n");
    }
    fprintf(fp,"%s","\n%");
    fprintf(fp,"FINAL PARAMETERS:\n\n");
    for(i = 0 ; i < npar ; i++){
      fprintf(fp,"%s","%");
      for(j=0;j<fixpar;j++)if(parlist[j]==i)break;
      if(j==fixpar)fprintf(fp,"par[%d]\t=\t%e\t+-\t%e\n",i,par[i],epar[i]);
      else       fprintf(fp,"par[%d]\t=\t%e\t+-\t(fixed)\n",i,par[i]);
    }
  }
  fprintf(fp,"\n\n");
  z = (xmax-xmin)/(2*ndata);
  xmin -= z; xmax += z;
  tx = (xmax-xmin) / 100. ;
  /* data and error bars */
  fprintf(fp,"# e c \"\\oc\"\n");
  for(i=0;i<ndata;i++){
    fprintf(fp,"%e %e %e\n",xdata[i],ydata[i],sigma[i]);
  }
  fprintf(fp,"# e0 c0\n");
  /* fitting function */
  for(i=0;i<101;i++){
    z=xmin + i*tx;
    fprintf(fp,"%e %e\n",z,F->f(z,par));
  }
  fclose(fp);
}

void BaseFitter::ReadParams(const char *file)
{


  FILE *file2 ;
  if((strcmp(file,"STNDIN")==0)||(strcmp(file,"stdin")==0))
    file2 =stdin ;
  else
    file2 = fopen(file,"r") ; 

  if(file2 == NULL) file2 = stdin ;

  //cerr<<"BaseFitter::ReadParams  "<<file2<<" "<<stdin<<endl ;  
  //cerr<<"BaseFitter::ReadParams  \n";
  //fprintf(stderr,"  %s\n",file) ;
  //cerr<<"BaseFitter::ReadParams  "<<file2<<" "<<stdin<<endl ;
  //fflush(stderr);

  F->f_init(file2); npar =F->Npar() ;
  
  bool interactive(file2==stdin) ;
  /*read in eps */
  if(interactive)fprintf(stderr,"Enter epsilon, maximum iterations\n");
  fscanf(file2,"%le",&eps);
  
  /*read in maxit */
  fscanf(file2,"%d",&maxit);
  
  /*read in display flag */
  if(interactive)fprintf(stderr,"Enter 1 for verbose, 0 for terse\n");
  fscanf(file2,"%d",&display);
  
  /* read in range of x to be fitted		*/
  if(interactive)fprintf(stderr,"Enter range of x to be fitted\n");
  fscanf(file2,"%le%le",&xlow,&xhigh);
  
  /*  read in number of parameters, allocate storage */
  /* f_init may have already set npar */
  if(npar < 0){
    if(interactive)fprintf(stderr,"Enter Number of parameters\n");
    fscanf(file2,"%d",&npar);
    F->setNpar(npar);
  }
  else {
    fprintf(stderr,"%d parameters\n",npar);
  }
  par= new double[npar];
  epar= new double[npar];
  parlist = new int[npar];
  
  /*  read in number of parameters to be held fixed */
  if(interactive)fprintf(stderr,"Enter number of parameters held fixed\n");
  fscanf(file2,"%d",&fixpar);
  printf("%d of %d  parameters are held fixed\n",fixpar,npar);
  for(int i=0;i<fixpar;i++){
    if(interactive)fprintf(stderr,
			   "Which parameter is fixed? 0 = first parameter\n");
    fscanf(file2,"%d",parlist+i);
  }
  
  /*  read in guess parameters */
  printf("\n\nInitial values:\n\n");
  for(int i = 0 ; i < npar ; i++){
    if(interactive)fprintf(stderr, "Enter initial value for parameter %d\n",i);
    if( 0 == fscanf(file2,"%le",par+i) ){
      printf("\nnot enough parameters\n");
      exit(0);
    }
    printf("par[%d]\t=\t%e\n",i,par[i]);
  }
  fclose(file2) ;
}

/*! compute inverse of dim by dim matrix at x, put result at y */
void BaseFitter::matinv(double *x,double *y, int dim){
double **a,*b,*c;
int     *perm,*perm2,i,j;
        perm=(int*)malloc(dim*sizeof(int));
        perm2=(int*)malloc(dim*sizeof(int));
        a=(double**)malloc(dim*sizeof(double*));
        b=(double*)malloc(dim*sizeof(double));
        c=(double*)malloc(dim*sizeof(double));
	a[0]=(double *)malloc(dim*dim*sizeof(double));
        for(i=0;i<dim;i++) a[i]=a[0]+(i*dim);
	for(i=0;i<dim*dim;i++)a[0][i]=x[i];
        if(0 ==factor(a,perm,perm2,dim)){ printf("singular matrix\n");
                exit(0);}
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++)b[j]=0.0;	b[i]=1.0;
        	subst(a,b,c,perm,perm2,dim);
		for(j=0;j<dim;j++)y[j*dim+i]=c[j];
	}
	free(a[0]); free(a); free(b); free(c); free(perm); free(perm2);
}


/*! solve the linear equation mat*ans = vec */
void BaseFitter::lineq(double *mat, double *vec,double *ans,int dim){
double **a;
int     *perm,*perm2,i;
        perm=(int*)malloc(dim*sizeof(int));
        perm2=(int*)malloc(dim*sizeof(int));
        a=(double**)malloc(dim*sizeof(double*));
	a[0]=(double *)malloc(dim*dim*sizeof(double));
        for(i=0;i<dim;i++) a[i]=a[0]+(i*dim);
	for(i=0;i<dim*dim;i++)a[0][i]=mat[i];
        if(0 ==factor(a,perm,perm2,dim)){ printf("singular matrix\n");
                exit(0);}
        subst(a,vec,ans,perm,perm2,dim);
	free(a[0]); free(a); free(perm); free(perm2);
}


/*!      Gaussian elimination with scaled TOTAL pivoting */
#define         max(a,b)        ( ((a)>(b))?(a) :(b) )
#define         abs(x)          ( ((x)>0)? (x) : -(x) )
int BaseFitter::factor(double **a, int  *p, int *q, int n)
{
int i,j,k,l,m,temp;
double  big,x,*d,epsilon;

        epsilon=n*1e-16;
                                /*      get space for row norms */
        if( NULL== ( d=(double*)malloc(n*sizeof(double)) )  ){ 
                printf("no space!!\n");
                return(0);
        }
        for(i=0;i<n;i++){
                p[i]=q[i]=i;    /*      initialize identity permutation */
                big=0;
                for(j=0;j<n;j++) big=max(big, abs(a[i][j]) );
                if(big==0){free(d);return(0);}
                d[i]=big;       /*      store row norms */
        }
        for(k=0;k<n-1;k++){
                j=k; m=k;
                big=0;
                for(i=k;i<n;i++)for(l=k;l<n;l++){
                        x=abs(a[p[i]][q[l]]/d[p[i]]);
                        if(x>big){
                                j=i;            /* find pivoting element*/
                                m=l;
                                big=x;
                        }
                }
                if(big==0){free(d);return(0);}
                        temp=p[k];
                        p[k]=p[j];      /* keep track of permutation    */
                        p[j]=temp;
                        temp=q[k];
                        q[k]=q[m];      /* keep track of permutation    */
                        q[m]=temp;
                for(i=k+1;i<n;i++){
                        /* store row multpliers in would-be zeroed location */
                        x=(a[p[i]][q[k]] /= a[p[k]][q[k]]);
                        for(j=k+1;j<n;j++) a[p[i]][q[j]] -= x*a[p[k]][q[j]];
                }
        }

        if(  abs(a[p[n-1]][q[n-1]]/d[p[n-1]])< epsilon)
	  {
	    free(d);
	    return(0);
	  }
	else
	  {
	    free(d);
	    return(1);
	  }
}
void BaseFitter::subst(double **a, double *b,double *x, int *p, int *q, int n)
{
int k,j;
double sum;
        for(k=0;k<n;k++){       /* foward elimination */
                sum=0;
                for(j=0;j<k;j++)sum += a[p[k]][q[j]] * x[q[j]];
                x[q[k]]=b[p[k]]-sum;
        }
        for(k=n-1;k>=0;k--){    /* back substitution */
                sum=0;
                for(j=k+1;j<n;j++)sum += a[p[k]][q[j]] * x[q[j]];
                x[q[k]]=(x[q[k]]-sum)/a[p[k]][q[k]];
        }
}

/*!	Minimization using variants of Newtons		*
 *	 and gradient descent methods			*
 *							*
 *	function to be minimized is 			*
 *							*
 *		phi( x )				*
 *							*
 *	whose gradient is given by			*
 *							*
 *		dphi( x, grad )				*
 *							*
 *	second partials are given by			*
 *							*
 *		ddphi( x , a )				*
 *							*
 *	minimization routine is				*
 *							*
 *		min(  )	                                *
 *							*/
double BaseFitter::min()
{
    double **a,*b,*y;
    int	*p,*q,i,count,flag;
    double *u,*grad,*xnew,step,nstep,grad_size,tol,phi_size,new_phi_size;
    double kno_tmp ;
    int kno_flag ;
    
    //Interface with Doug's code with out having to modify anything
    int n(npar);
    double *x=par ;
    //

    p=(int*)malloc(n*sizeof(int));
    q=(int*)malloc(n*sizeof(int));
    a=(double**)malloc(n*sizeof(double*));
    b=(double*)malloc(n*sizeof(double));
    y=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++) a[i]=(double*)malloc(n*sizeof(double));
    grad=(double*)malloc(n*sizeof(double));
    u=(double*)malloc(n*sizeof(double));
    xnew=(double*)malloc(n*sizeof(double));
    count=0;
    nstep=1;
    step=1;
    tol=1.00000001;
    phi_size=phi(x);
    dphi( x , grad );
    new_phi_size = phi_size ;
    while( (grad_size=norm(grad,n)) > eps && count <= maxit  ){
	if(display){
	    printf("iteration %d\n\t|phi| =\t%e\t|grad(phi)| =\t%e\n"
		,count,phi_size,grad_size);
	    for(i=0;i<n;i++)printf("\tx[%d] =\t%e\tgrad[%d] =\t%e\n"
		,i,x[i],i,grad[i]);
	}
	ddphi( x , a );
	if(0 !=(flag=factor(a,p,q,n))) {
	    subst(a,grad,y,p,q,n);
	    if( nstep*dot( y , grad ,n) < 0.0 ) nstep = - nstep;
	    /* Newton iteration	*/
	    for(i=0;i<n;i++){
		xnew[i]= x[i] - nstep* y[i];
	    }
/**
printf("min: y: (");
for(i=0;i<n;i++)
  printf("%g ", y[i]);
printf(")\n") ;
printf("min: xnew: (");
for(i=0;i<n;i++)
  printf("%g ", xnew[i]);
printf(")\n") ;
printf("min: phi(xnew): %g\n",phi(xnew));
**/
 
 kno_tmp = phi(xnew) ;
 kno_flag = 0 ;
 while( isnan(kno_tmp) || isinf(kno_tmp) ){ 
   /*take a smaller step*/
   /* we probably hit a very large argument for an exponential */
   printf("min: OOPS! We hit infinity or NaN in function evaluation. ") ;
   printf("Reducing nstep: %g\n",nstep );
    nstep /= 2.0;
    for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
    kno_tmp = phi(xnew) ;
    kno_flag ++;
    if(kno_flag>100)
      {
	printf("min: Still hitting  infinity or NaN\n I give up!!\n") ;
	exit(0);
      }
 }

	    if((new_phi_size=phi(xnew))*tol <phi_size){
		/* Try to take a bigger step */
		do {
		  /**
printf("min: bigger %g\n",nstep );
		  **/
		    phi_size=new_phi_size;
		    nstep *= 2;
if(fabs(nstep) > 10){printf("BREAKING-BIG, phi = %e\n",phi_size); break;}
		    for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
		} while( tol*phi_size > (new_phi_size=phi(xnew)) );
		new_phi_size = phi_size;
		nstep /= 2;
		for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
	    }
	    else{
		while( !( phi_size*tol > (new_phi_size=phi(xnew)))){
		    /* if phi increased take a smaller step */
		    /* comparison is !(cond) so will take correct action
		      if new_phi_size is NaN, for which comparisons fail */
		  /**
printf("min: smaller %g\n",nstep );
		  **/
		    nstep /= 2.0;
if(fabs(nstep) < 1e-4){printf("BREAKING-SMALL, phi = %e\n",phi_size); break;}
		    for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
		}
	    }

	    if( display ) printf("\tNewton descent\t\tnstep =\t%e\n",nstep);
	}
	else
	{	/* Gradient descent	(Last resort!!!)	*/
	  /**
printf("min: hit a singular matrix\n");
printf("min: grad: (");
for(i=0;i<n;i++){
  printf("%g ",grad[i]);
}  
printf(")\n");
	  **/
/**exit(0);**/
	    for(i=0;i<n;i++){
		u[i] = -grad[i] ;
		xnew[i]= x[i] + step * u[i];
	    }
	    if( (new_phi_size=phi(xnew))*tol <phi_size){
		do {	/* try to take a bigger byte! */
printf("min: taking bigger bite\n");
		    phi_size=new_phi_size;
		    step *= 2;
		    for(i=0;i<n;i++) xnew[i]= x[i] + step * u[i];
		}
		while(  phi_size*tol > (new_phi_size=phi(xnew)) );
		new_phi_size = phi_size;
		step /=2;
		for(i=0;i<n;i++) xnew[i]= x[i] + step * u[i];
	    }
	    else {
		while( tol*phi_size < (new_phi_size=phi(xnew)) ){
printf("min: taking smaller bite\n");
		    /* better take a smaller byte! */
		    step /= 2;
		    for(i=0;i<n;i++) xnew[i]= x[i] + step * u[i];
		}
	    }
	    if( display ) printf("\tgradient descent\t\tstep =\t%e\n",step);
	}
	for(i=0;i<n;i++) x[i]= xnew[i] ;

	/**
printf("min: x: (");
for(i=0;i<n;i++)
  printf("%g ", x[i]);
printf(")\n") ;
	**/
	phi_size=new_phi_size;
	dphi( x , grad );
	count++;
    }
    if(display){
	if(count >= maxit) {
	  printf("No ");
	  convergence = false ;
	}
	else{
	  convergence = true ;
	}
	printf("convergence after %d iterations\n|phi| = %e\t|grad(phi)| = %e\n"
	    ,count-1,phi_size,grad_size);
	printf("final parameters:\n");
	for(i=0;i<n;i++)printf("\tx[%d] =\t%e\tgrad[%d] =\t%e\n"
	    ,i,x[i],i,grad[i]);
	printf("\n");
    }

    free(p);
    free(q); 
    for(i=0;i<n;i++) free(a[i]) ; 
    free(a);
    free(b);
    free(y);
    free(grad);
    free(u);
    free(xnew);

    return(new_phi_size);
}


double BaseFitter::norm(double *x, int n)
{
    double sum;
    int i;
    sum=0;
    for(i=0;i<n;i++) sum+=x[i]*x[i];
    return( sqrt(sum) );
}
double BaseFitter::dot(double *x, double *y, int n)
{
    double sum;
    int i;
    sum=0;
    for(i=0;i<n;i++) sum+=x[i]*y[i];
    return( sum );
}


void BaseFitter::getParams(Array<Double>& p, Array<Double>& e){
  p.resize(npar);
  e.resize(npar);
  for(int i(0);i<npar;i++){
    p[i] = par[i] ;
    e[i] = epar[i] ;
  }
}

void BaseFitter::setParams(const Array<Double>& p ){
  if(p.size()!=npar){
    cerr<<"OOPS! "<<npar<<" != "<< p.size();
    exit(1001);
  }
  for(int i(0);i<npar;i++){
    par[i] = p[i] ;
  }
}

void BaseFitter::setFixedParams(const Array<int>& f ){
  fixpar = f.size() ;
  if(f.size()>npar){
    cerr<<"OOPS! "<<npar<<" < "<< f.size();
    exit(1002);
  }
  for(int i(0);i<fixpar;i++){
    parlist[i] = f[i] ;
  }
}

} // namespace CovFit
