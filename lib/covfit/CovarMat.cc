//$Id: CovarMat.cc,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: CovarMat.cc,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.3  2007/03/02 06:08:04  kostas
//fixed up the covariance matrix and the proplists so that they can handle the complex numbers hadspec is producing...
//
//Revision 1.2  2005/06/03 01:49:12  kostas
//jackknifed covariance matrix
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
//Revision 1.5  2004/05/19 22:11:29  kostas
//fixed another bug
//
//Revision 1.4  2004/05/19 21:53:12  kostas
//fixed bug
//
//Revision 1.3  2004/05/19 21:41:44  kostas
//added some Arrays in CovarMat
//
//Revision 1.2  2004/02/13 15:37:58  kostas
//Just comments
//

//
// This is just a c++ wrap around the c fitting codes provided by
// Doug Toussaint. A description of the algorithms used can be
// found in Doug's TASI lectures.
//
#include "covfit/covarmat.h"
#include <iostream>

namespace CovFit {

/*!
  read a list of propagators, compute average propagator and
   covariance matrix.
   Propagator is list of numbers
      0   value
      1   value
      2   value
      ...
      max value
      
   repeated N times.
   The code allows for a DeGrand header.
*/
void CovarMat::ReadDataList(char *file)
{ 
  //cerr<<"Entered CovarMat::ReadDataList(char *file)"<<endl ;
  FILE *fp ;
  if((strcmp(file,"STDIN")==0)||(strcmp(file,"stdin")==0))
    fp = stdin ;
  else
    fp = fopen(file,"r") ;
  //cerr<<"fp="<<fp<<"  stdin="<<stdin<<endl ;
  int ndata(0);
  double xx,yy;
  int f,fo,foo ;
  f=0 ;
  double tt ;
  bool check_adat(true);
  while(fscanf(fp,"%lg%le",&xx,&yy) == 2 ){
    if(check_adat){
      check_adat = false ;
      if(xx!=0) // DeGrand header found.... Ignore it!
	fscanf(fp,"%i%i%i",&f,&fo,&foo) ; // read the crap...
      //if f == 1 then complex numbers..
      else{ // we have to keep the data
	if(f==1) fscanf(fp,"%lg",&tt) ; //get rid of the imaginary part
	xdata[ndata]=xx ; // no checks needed.... it's the first time!
	ylist = (double *)realloc(ylist, (ndata+1)*sizeof(double));
	ylist[ndata]=yy;
	ndata++;
      }
    }
    else{
      if(f==1) fscanf(fp,"%lg",&tt) ; //get rid of the imaginary part
      //cerr<<"CovarMat::ReadDataList(char *file)  "<<xx<<" "<<yy<<endl ;
      if(ndata<npoints)
	xdata[ndata]=xx ;
      else
	if(fabs(xx-xdata[ndata%npoints])>1.0e-15){
	  /* unexpected entry - give no mercy */
	  fprintf(stderr,"Bad data in propagator %d\n",ndata/npoints);
	  exit(1);
	}
      ylist = (double *)realloc(ylist, (ndata+1)*sizeof(double));
      ylist[ndata]=yy;
      ndata++;
    }
  }

  fclose(fp);
  nprops = ndata/npoints ;
  //for(int i=0;i<ndata;i++)
  //  cerr<<ylist[i]<<endl;
}

void CovarMat::SetDataList(double *xd, double *ll,int length)
{
  if(length%npoints !=0 ){
    fprintf(stderr,"y-Lenght not good!\n");
    exit(1);
  }
  nprops = length/npoints ;

  for(int i(0);i<npoints;i++){
    xdata[i]=xd[i];
  }

  ylist = (double *)realloc(ylist, length*sizeof(double));
  for(int i(0);i<length;i++){
    ylist[i] = ll[i];
  }
}

void CovarMat::SetDataList(const Array<double>& xd, const PropList& ll)
{
  if(xd.size() != npoints){
    fprintf(stderr,"x-Lenght not good!\n");
    exit(1);
  }
  if(ll[0].size() != npoints){
    fprintf(stderr,"y-Lenght not good!\n");
    exit(1);
  }

  nprops = ll.size() ; 
  int length (nprops*npoints);

  for(int i(0);i<npoints;i++){
    xdata[i]=xd[i];
  }

  ylist = (double *)realloc(ylist, length*sizeof(double));
  for(int i(0);i<length;i++){
    ylist[i] = ll[i/npoints][i%npoints];
  }
}

void CovarMat::CalcCovarMat(){
  if(nprops<=0) return ; // Nothing to compute
  int i,j,dist ;
  
  double *current = (double*)malloc(npoints*sizeof(double));
  double *block   = (double*)malloc(npoints*sizeof(double));

  double *average = ydata ;
  
  /* initialize block, average and covariance */
  for(i=0;i<npoints;i++){
    block[i] = average[i]=0.0;
    for(j=0;j<npoints;j++)covar[i*npoints+j]=0.0;
  }

  int ndata = nprops*npoints ;
  /* read propagators, accumulate sums, check for proper format and order */
  int nprop(0);
  nblocks=0;
  for (int k(0);k<ndata;){
    for(dist=0;dist<npoints;dist++)
      current[dist] = ylist[k++] ;
    for(i=0;i<npoints;i++)
      block[i] += current[i];
    
    nprop++;
    //cerr<<"Doing prop "<<nprop<<endl ;
    if(nprop % blocksize == 0){	/* we have finished a block */
      for(i=0;i<npoints;i++)
	block[i] /= (double)blocksize;
      for(i=0;i<npoints;i++){
	average[i] += block[i];
	for(j=0;j<npoints;j++)
	  covar[i*npoints+j] += block[i]*block[j];
      }
      for(i=0;i<npoints;i++)block[i] = 0.0;
      nblocks++ ;
    }
  }
  
  free(ylist);
  ylist = (double *)malloc(sizeof(double)); // be ready for new data

  fprintf(stderr,"%d data read in, %d blocks\n",nprops,nblocks);
  
  /* normalize and subtract disconnected parts */
  /* covariance matrix is normalized for "error of mean" */
  /* no attempt is made to correct for difference between the
     sample average and the real average - subsequent fitting
     programs must take care of this. */
  double x(1.0/(double)nblocks);
  for(i=0;i<npoints;i++) average[i] *= x;
  for(i=0;i<npoints;i++)
    for(j=0;j<npoints;j++) 
      covar[i*npoints+j] =
	(x*covar[i*npoints+j] - average[i]*average[j])/((double)(nblocks));
  
  free(current);
  free(block);
}

void CovarMat::CalcJackCovarMat()
{
  CalcCovarMat() ;
  //scale the matrix with N^2 I think this is correct. Need to check
  double factor(nblocks*nblocks) ;
  //printf("CovarMat::CalcJackCovarMat() %g\n",factor);
  for(int i=0;i<npoints;i++)
    for(int j=0;j<npoints;j++) 
      {
//printf("CovarMat::CalcJackCovarMat(): %i %i\n",i,j);
	covar[i*npoints+j] *= factor ;
      }
}

/*!
  Output is
    PROPAGATOR
    i  average[i]
    COVARIANCES
    i j covar[i,j]
   Arguments are maximum distance of propagator and blocking factor
*/
void CovarMat::WriteCovarMat(char *file){
  if(npoints<=0) return ; // Nothing to write
  FILE *fp ;
  if((strcmp(file,"STDOUT")==0)||(strcmp(file,"stdout")==0))
    fp = stdout ;
  else
    fp = fopen(file,"w") ;

  /* print results */
  fprintf(fp,"DATA\n");
  for(int i=0;i<npoints;i++)
    fprintf(fp,"%f\t%.15e\n",xdata[i],ydata[i]);
  fprintf(fp,"COVARIANCE:  %d blocks\n",nblocks);
  for(int i=0;i<npoints;i++)for(int j=0;j<npoints;j++)
    fprintf(fp,"%d\t%d\t%.15e\n",i,j,covar[i*npoints+j]);
  fclose(fp);
}

/*!
Form of input:							*
 "DATA"	(or some other word)
 list of function points:
	x y
 "COVARIANCE"  blocks   "blocks" (or some other words), where blocks
     is the number of blocks used in computing the covariance matrix.
 list of covariance matrix elements:
	x1 x2 covar12
 Missing covariance matrix elements are taken to be zero,
  except diagonal ones must be specified.
 Extra covariance matrix elements are ignored.
*/
void CovarMat::ReadCovarMat(char *file){
  int i,j;

  FILE *fp = fopen(file,"r") ;
  if(fp == NULL){
    fprintf(stderr,"Cannot open CovarMat file %s\n",file);
    exit(1);
  }
  
  fscanf(fp,"%*s");	/* flush header word */
  npoints=0;
  double xx,yy,zz;
  while(fscanf(fp,"%le %le",&xx,&yy) == 2 ){
    xdata = (double *)realloc(xdata,(npoints+1)*sizeof(double));
    ydata = (double *)realloc(ydata,(npoints+1)*sizeof(double));
    xdata[npoints] = xx;   ydata[npoints] = yy;
    npoints++;
  }
  /* read covariances */
  fscanf(fp,"%*s%d%*s",&nblocks);	/* flush header words */
  fprintf(stderr,"\n\t read in %d data points\n",npoints);
  if( nblocks <= npoints ){
    fprintf(stderr,"Not enough blocks in covariance matrix\n");
    exit(0);
  }
  /* initialize matrix */
  covar = (double*)realloc(covar,npoints*npoints*sizeof(double));

  for(i=0;i<npoints;i++){
    for(j=0;j<npoints;j++)covar[i*npoints+j] = 0.0;
    covar[i*npoints+i] = -1.0;
    /*diagonal can not be < 0, use to make sure it is specified */
  }
  while(fscanf(fp,"%d%d%le",&i,&j,&zz) == 3){
    covar[i*npoints+j] = zz;
  }
  fclose(fp);
  /* check to make sure diagonal elements were specified */
  for(i=0,j=0;i<npoints;i++){
    if(covar[i*npoints+i] < 0.0){
      fprintf(stderr,
	      "Diagonal covariance missing or negative: x = %e\n",xdata[i]);
      j=1;
    }
  }
  if(j)exit(1);
}

} // namespace CovFit
