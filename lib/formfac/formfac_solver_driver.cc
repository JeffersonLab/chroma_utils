// $Id: formfac_solver_driver.cc,v 2.0 2008/12/05 04:43:37 edwards Exp $
/*! \file
 * \brief Driver for call the LLSq solver
 */

#include "formfac/formfac_solver_driver.h"
#include <iostream>

namespace FF
{
  //! Driver
  void llsq_driver(XMLWriter& xml_out, LLSqComponent<EnsemReal>& sys, const std::list<QsqVal>& qsq_val,
		   double cov_tol, double solve_tol)
  {
    // List of the unknowns. These are the column indices for the system solver.
    const std::vector<std::string> unknowns = sys.unknowns();

    //
    // Loop over each Qsq point
    //
    push(xml_out, "FormFactors");
    write(xml_out, "unknowns", unknowns);

    int n=0;
    for(std::list<QsqVal>::const_iterator pqsq = qsq_val.begin();
	pqsq != qsq_val.end();
	++pqsq,++n)
    {
      const QsqVal& qsq = *pqsq;
    
      push(xml_out, "elem");     // next array element
      write(xml_out, "n", n);
      write(xml_out, "Qsq", qsq.Qsq_c);
      write(xml_out, "Momenta", qsq);

      std::cout << "n=" << n << "  Qsq= " << qsq.Qsq_c << std::endl;
      sys.setQsq(qsq);

      std::cout << "llsq" << std::endl;
    
      push(xml_out, "Fit");

      const std::vector<TimeFitRange_t> fit_ranges = sys.getFitRanges();

      for(int tcnt=0; tcnt < fit_ranges.size(); ++tcnt)
      {
	push(xml_out, "elem");     // next array element

	const TimeFitRange_t& fit_range = fit_ranges[tcnt];
	sys.setTimeFit(fit_range);

        int t_i = fit_range.t_i;
        int t_f = fit_range.t_f;

	LLSqResults<EnsemReal> ret = llsq(sys, cov_tol, solve_tol);
	if (ret.count_singular_val > 0)
	  std::cerr << "Found " << ret.count_singular_val << " singular values at Qsq= " << qsq.Qsq_c 
                    << "  t_i= " << t_i << " t_f= " << t_f << std::endl;
      
	write(xml_out, "t_i", t_i);
	write(xml_out, "t_f", t_f);
	write(xml_out, "count_singular_val", ret.count_singular_val);
	write(xml_out, "chisq", ret.chisq);
	write(xml_out, "ndof", ret.ndof);
	write(xml_out, "Q", ret.Q);

	printf("t_i= %3d t_f= %3d   sing= %5d  chisq= %-13.5g  ndof= %4d  Q= %-10.2g\n",
	       t_i, t_f, ret.count_singular_val, toDouble(ret.chisq), ret.ndof, toDouble(ret.Q));

	for(int n=0; n < unknowns.size(); ++n)
	{
	  const EnsemReal& F = ret.x[n];

	  //  char ss[1024];
	  // sprintf(ss,"FF_n%d_ti%d_tf%d", n,t_i,t_f);
	  //  write(string(ss),F);

	  Real f = toFloat(mean(F));
	  Real e = toFloat(sqrt(variance(F)));

          printf("t_i= %2d t_f= %2d     F[ %2d ]= %- 13.5g +/- %- 13.5g\n",
                 t_i, t_f, n, toDouble(f), toDouble(e));

//	  std::cout << "t_i= " << t_i << " t_f= " << t_f 
//                    << "\tF[ " << n << " ]= " << f << " err= " << e 
//                    << "     chisq= " << c << std::endl;

	  push(xml_out, "elem");
	  write(xml_out, "name", unknowns[n]);
	  write(xml_out, "F", f);
	  write(xml_out, "F_err", e);
	  pop(xml_out);
	}
	//write out the correlation matrix for the fit params too
	Array2d<EnsemReal> covF = ret.x_cov;
	//only interested in the mean - normalise by the diags
	std::cout << "parameter correlations:" << std::endl;
	for(int n=0; n < unknowns.size(); ++n)
	  {
	    Real diag_n = sqrt( mean( covF[n][n] ) );
	    for(int m=0; m <= n; ++m)
	      {
		Real diag_m = sqrt( mean( covF[m][m] ) );
		double corr = toDouble( mean( covF[n][m] ) / ( diag_n * diag_m) );
		if( toDouble(diag_n*diag_m) == 0.0 ){ corr = 0.0;}
		std::cout << corr << "  ";
	      }
	    std::cout << std::endl;
	  }


	// Possibly save the file
	sys.save(ret.x, qsq.Qsq_c, t_i, t_f);

	pop(xml_out);  // elem
      }

      pop(xml_out);    // Fit
      pop(xml_out);    // elem
    } // end for(qsq)

    pop(xml_out);  // FormFactors
  }




}
