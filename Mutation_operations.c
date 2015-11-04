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

double uniformrng( void ) ;

double normalrng( void ) ;

int integerrng( int, int ) ;

void uniformdirectionrng(double *, int) ;

double cost(double*,int) ;

void self_adj_index_search(int*,double,double*,int) ;

void Mutation_HitAndRun(double *z, double *fz,
				int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, double scl, double *accpr_pop,
				double *z_new){

	/* THIS IS THE HIT AND RUN METROPOLIS HASTINGS UPDATE */

	/* Smith, Robert L. "Efficient Monte Carlo procedures for generating points
	 * uniformly distributed over bounded regions." Operations Research 32.6
	 * (1984): 1296-1308. */

	/*Chen, Ming-Hui, and Bruce Schmeiser. "Performance of the Gibbs, hit-and-run,
	 * and Metropolis samplers." Journal of computational and graphical
	 * statistics 2.3 (1993): 251-272.*/

	int k_old ;
	int k_new ;
	int i;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */
	uniformdirectionrng( z_new , N_dimension ) ;
	un = normalrng()*scl ;
	for ( i = 1 ; i <= N_dimension ; i++ )
		z_new[i] = z[i] + z_new[i]*un ;

	fz_new = cost( z_new, N_dimension ) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}

/*K POINT OPERATION*/

void Mutation_Kpoint(double *z, double *fz,
				int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, double scl, double *accpr_pop,
				double *z_new){

	int k_old ;
	int k_new ;
	int i;
	int lab_1 ;
	int lab_2 ;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */

	for ( i=1 ; i<=N_dimension ; i++) z_new[i] = z[i] ;

	/*1st dimension*/
	lab_1 = integerrng(1,N_dimension) ;
	un = normalrng()*scl ;
	z_new[lab_1] += un ;

	/*2nd dimension*/
	if ( (uniformrng()<0.5) && (N_dimension>=2) ){
		lab_2 = integerrng(1,N_dimension-1) ;
		if ( lab_2 >= lab_1 ) lab_2++ ;
		un = normalrng()*scl ;
		z_new[lab_2] += un ;
	}

	fz_new = cost( z_new, N_dimension ) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}

/*METROPOLIS RANDOM WALK OPERATION*/

void Mutation_Metropolis(double *z, double *fz,
				int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, double scl, double *accpr_pop,
				double *z_new){

	int k_old ;
	int k_new ;
	int i;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */

	for ( i = 1 ; i <= N_dimension ; i++ )
		z_new[i] = z[i] +normalrng()*scl ;

	fz_new = cost( z_new, N_dimension ) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}









