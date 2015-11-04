/*
 * Copyright 2014 Georgios Karagiannis
 *
 * Georgios Karagiannis
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 765 494-3405
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
*/

#include <math.h>
#include "nrutil.h"

void uniformdirectionrng(double *d, int n) ;
double uniformrng(void) ;
double normalrng(void) ;
double cost(double*,int) ;

void MH_HitAndRun( double *z, double *fz,
					int N_dimension, double temp,
					double scl, double *accpr_pop, double *z_new )
{

	int i ;
	double accpr ;
	double fy ;
	double r ;
	double un ;

	uniformdirectionrng( z_new , N_dimension ) ;
    for(i=1; i<=N_dimension; i++)
    	z_new[i] = z[i] + z_new[i]*scl*normalrng() ;

	fy = cost(z_new, N_dimension) ;

	r = (*fz-fy) / temp ;
	accpr = ( r>0.0 ? 1.0 : exp(r) ) ;
	un = uniformrng() ;
	if(accpr>un)
	{
		*fz = fy;
		for(i=1 ; i<=N_dimension ; i++) z[i] = z_new[i] ;
	}

	*accpr_pop = accpr ;

}
