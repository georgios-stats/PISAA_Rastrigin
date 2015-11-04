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

int integerrng(int, int) ;

double normalrng( void ) ;

double cost(double*,int) ;

void self_adj_index_search(int*,double,double*,int) ;

void CO_select_forward_0(double *prob, int *n1, int *n2, int N_population){

	/* draw n1 */

	*n1 = integerrng(1, N_population) ;

	/* draw n2|n1 */

	*n2 = integerrng(1, N_population-1) ;
	if ( *n2 >= *n1 ) (*n2)++ ; /*{1,...,N_population}-{n1}*/

	/* compute probability */

	*prob = -log((double) N_population-1) -log((double) N_population) ;

}

void CO_select_backward_0(double *prob, int N_population){
	*prob = -log((double) N_population-1) -log((double) N_population) ;
}

void CO_select_forward_1(double *prob, int *n1, int *n2,
				double *fx, int N_population, double temp){

/*	IN ORDER TO OVERCOME THE PDECISION ERRORS IN EXP(),
			WE APPLY THE MINUS-MAXIMUM TRICK,
			OR OTHERWISE THE PLUS-MINIMUM FOR THE NEGATIVE VALUES*/

	int i ;
	double fx_min ;
	double un ;
	double fx_sum ;
	double ss ;

	/* draw n1 */

	*n1 = integerrng(1, N_population) ;

	/* draw n2|n1 */

	fx_min = ( (*n1==1) ? fx[2] : fx[1] ) ;
	for ( i=1; i <= N_population; i++ )
		if ( (i != *n1) && fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i=1; i <= N_population; i++ )
		if ( i != *n1 ) fx_sum += exp( ( -fx[i] + fx_min ) / temp ) ;

	if (fx_sum == 0.0){
		*n2 = integerrng(1, N_population-1) ;
		if ( *n2 >= *n1 ) (*n2)++ ; /*{1,...,N_population}-{n1}*/
	}
	else{
		un = uniformrng() ;
		ss = 0.0 ; *n2 = 0 ;
		do {
			(*n2)++ ;
			if( *n2 != *n1 ) ss += exp((-fx[*n2]+fx_min)/temp -log(fx_sum));
		} while( ss < un && *n2 <= N_population ) ;
	}

	/* compute probability */

	if (fx_sum == 0.0)
		*prob = -log((double) N_population-1)
					-log((double) N_population) ;
	else
		*prob = (-fx[*n2]+fx_min)/temp -log(fx_sum)
					-log((double) N_population) ;

}

void CO_select_backward_1(double *prob, int n1, int n2,
				double *fx, int N_population, double temp){

	/*	IN ORDER TO OVERCOME THE PDECISION ERRORS IN EXP(),
				WE APPLY THE MINUS-MAXIMUM TRICK,
				OR OTHERWISE THE PLUS-MINIMUM FOR THE NEGATIVE VALUES*/

	int i ;
	double fx_min ;
	double un ;
	double fx_sum ;
	double ss ;

	fx_min = ( (n1==1) ? fx[2] : fx[1] ) ;
	for ( i=1; i <= N_population; i++ )
		if ( (i != n1) && (fx[i]<fx_min) ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1; i <= N_population; i++ )
		if ( i != n1 ) fx_sum += exp((-fx[i]+fx_min)/temp) ;

	if (fx_sum == 0.0)
		*prob = -log((double) N_population-1)
					-log((double) N_population) ;
	else
		*prob = (-fx[n2]+fx_min)/temp -log(fx_sum)
							-log((double) N_population) ;

}

void CO_select_forward_2(double *prob, int *n1, int *n2,
							double *fx, int N_population, double temp){

	/*	IN ORDER TO OVERCOME THE PDECISION ERRORS IN EXP(),
				WE APPLY THE MINUS-MAXIMUM TRICK,
				OR OTHERWISE THE PLUS-MINIMUM FOR THE NEGATIVE VALUES*/

	int i ;
	double fx_min ;
	double fx_min_22 ;
	double un ;
	double fx_sum ;
	double fx_sum_22 ;
	double ss ;

	double prob1 ;
	double prob2 ;
	double prob11 ;
	double prob12 ;
	double prob21 ;
	double prob22 ;

	/* get n1 ; comp Pr(n1) */

	fx_min = fx[1] ;
	for ( i=1 ; i <= N_population ; i++ )
		if ( fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1 ; i <= N_population ; i++ )
		fx_sum += exp( ( -fx[i] +fx_min ) / temp ) ;

	/*... store*/
	fx_min_22 = fx_min ;
	fx_sum_22 = fx_sum ;

	if ( fx_sum != 0.0 ){
		un = uniformrng() ;
		ss = 0.0 ; *n1 = 0 ;
		do {
			(*n1)++ ;
			ss += exp( (-fx[*n1]+fx_min)/temp -log(fx_sum) ) ;
		} while( ss < un && *n1 <= N_population ) ;

		prob11 = (-fx[*n1]+fx_min) / temp - log(fx_sum) ;
	}
	else{
		*n1 = integerrng(1, N_population) ;

		prob11 = -log( (double) N_population ) ;
	}

	/* get n2|n1 ; comp Pr(n2|n1) */

	fx_min = ( (*n1 == 1) ? fx[2] : fx[1] ) ;
	for ( i=1 ; i <= N_population ; i++ )
		if ( i != *n1 && fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1 ; i <= N_population ; i++ )
		if ( i != *n1) fx_sum += exp( ( -fx[i]+fx_min ) / temp ) ;

	if ( fx_sum != 0.0 ){
		un = uniformrng() ;
		ss = 0.0 ; *n2 = 0 ;
		do {
			(*n2)++ ;
			if( *n2 != *n1 ) ss += exp((-fx[*n2]+fx_min)/temp -log(fx_sum)) ;
		} while ( ss<un && *n2<=N_population ) ;

		prob12 = (-fx[*n2]+fx_min)/temp -log(fx_sum) ;
	}
	else{
		*n2 = (int) integerrng(1, N_population-1) ; /*{1,...,pop-1}*/
		if ( *n2 >= *n1 ) (*n2)++ ; /*{1,...,pop}-{k1}*/

		prob12 = -log( (double) (N_population-1) ) ;
	}

	/*compute Pr(n1)Pr(n2|n1)*/

	prob1 = prob11 + prob12 ;

	/* comp Pr(n2) */

	fx_min = fx_min_22 ;
	fx_sum = fx_sum_22 ;

	if ( fx_sum != 0.0 )
		prob22 = (-fx[*n2]+fx_min) / temp - log(fx_sum) ;
	else
		prob22 = -log( (double) N_population ) ;

	/* comp Pr(n1|n2) */

	fx_min = ( (*n2 == 1) ? fx[2] : fx[1] ) ;
	for ( i=1 ; i <= N_population ; i++ )
		if ( i != *n2 && fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1 ; i <= N_population ; i++ )
		if ( i != *n2) fx_sum += exp( (-fx[i]+fx_min ) / temp ) ;

	if ( fx_sum != 0.0 )
		prob21 = (-fx[*n1]+fx_min) / temp - log(fx_sum) ;
	else
		prob21 = - log( (double) (N_population-1) ) ;

	/* compute Pr(n2)Pr(n1|n2) */

	prob2 = prob22 + prob21 ;

	/* compute the whole probability */

	if ( prob1 > prob2 )
		*prob = prob1 + log( 1.0 + exp(prob2-prob1) ) ;
	else
		*prob = prob2 + log( 1.0 + exp(prob1-prob2) ) ;

}

void CO_select_backward_2(double *prob, int n1, int n2,
					double *fx, int N_population, double temp){

	/*	IN ORDER TO OVERCOME THE PDECISION ERRORS IN EXP(),
				WE APPLY THE MINUS-MAXIMUM TRICK,
				OR OTHERWISE THE PLUS-MINIMUM FOR THE NEGATIVE VALUES*/

	int i ;
	double fx_min ;
	double fx_sum ;

	double prob1 ;
	double prob2 ;
	double prob11 ;
	double prob12 ;
	double prob21 ;
	double prob22 ;

	/* with all */

	fx_min = fx[1] ;
	for ( i=1 ; i <= N_population ; i++ )
		if ( fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1 ; i <= N_population ; i++ )
		fx_sum += exp( ( -fx[i] +fx_min ) / temp ) ;

	if (fx_sum!=0.0)
		prob11 = (-fx[n1]+fx_min)/temp -log(fx_sum) ;
	else
		prob11 = -log( (double) N_population ) ;

	if (fx_sum!=0.0)
		prob22 = (-fx[n2]+fx_min)/temp -log(fx_sum) ;
	else
		prob22 = -log( (double) N_population ) ;

	/* n2|n1 */

	fx_min = ( (n1==1) ? fx[2] : fx[1] ) ;
	for ( i=1 ; i <= N_population ; i++ )
		if ( i != n1 && fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1 ; i <= N_population ; i++ )
		if ( i != n1) fx_sum += exp( ( -fx[i] +fx_min ) / temp ) ;

	if (fx_sum!=0.0)
		prob12 = (-fx[n2]+fx_min)/temp -log(fx_sum) ;
	else
		prob12 = -log( (double) (N_population-1) ) ;

	/* n1|n2 */

	fx_min = ( (n2==1) ? fx[2] : fx[1] ) ;
	for ( i=1 ; i <= N_population ; i++ )
		if ( i != n2 && fx[i] < fx_min ) fx_min = fx[i] ;

	fx_sum = 0.0 ;
	for ( i = 1 ; i <= N_population ; i++ )
		if ( i != n2) fx_sum += exp( (-fx[i] +fx_min ) / temp ) ;

	if (fx_sum!=0.0)
		prob21 = (-fx[n1]+fx_min)/temp -log(fx_sum) ;
	else
		prob21 = -log( (double) (N_population-1) ) ;

	prob1 = prob11 + prob12 ;
	prob2 = prob21 + prob22 ;

	if ( prob1 > prob2 )
		*prob = prob1 + log( 1.0 + exp(prob2-prob1) ) ;
	else
		*prob = prob2 + log( 1.0 + exp(prob1-prob2) ) ;

}

void CO_select_forward_3(double *prob, int *n1, int *n2,
			double *fx, int N_population, double temp){

/*	IN ORDER TO OVERCOME THE PDECISION ERRORS IN EXP(),
WE APPLY THE MINUS-MAXIMUM TRICK,
OR OTHERWISE THE PLUS-MINIMUM FOR THE NEGATIVE VALUES*/

	int i ;
	int j ;
	double pr_max ;
	double pr_sum ;
	double pr_cum ;
	double un ;

	/*get max*/
	pr_max = (-fx[2]-fx[1]) ;
	for (i=2; i<=N_population; i++)
		for (j=1; j<=i-1; j++)
			if ( (-fx[i]-fx[j]) > pr_max ) pr_max = (-fx[i]-fx[j]) ;

	/*get sum*/
	pr_sum = 0.0 ;
	for (i=2; i<=N_population; i++)
		for (j=1; j<=i-1; j++)
			pr_sum += exp( (-fx[i]-fx[j])/temp -pr_max/temp ) ;

	/*get the population labels*/
	if (pr_sum>0.0) {
		un = uniformrng() ;
		pr_cum = 0.0 ;
		for (i=2; i<=N_population; i++)
			for (j=1; j<=i-1; j++){
				pr_cum += exp( (-fx[i]-fx[j])/temp -pr_max/temp )/pr_sum ;
				if (pr_cum>=un) goto labels_found ;
			}
		labels_found :
		*n1 = i ;
		*n2 = j ;
	}
	else {
		*n1 = integerrng(1, N_population) ;
		*n2 = (int) integerrng(1, N_population-1) ; /*{1,...,pop-1}*/
	}

	/*get the probability*/
	if (pr_sum!=0)
		*prob = exp( (-fx[*n1]-fx[*n2])/temp -pr_max/temp )/pr_sum ;
	else
		*prob = log(2.0) -log(N_population) -log(N_population-1.0) ;
}

void CO_select_backward_3(double *prob, int n1, int n2,
			double *fx, int N_population, double temp){

/*	IN ORDER TO OVERCOME THE PDECISION ERRORS IN EXP(),
WE APPLY THE MINUS-MAXIMUM TRICK,
OR OTHERWISE THE PLUS-MINIMUM FOR THE NEGATIVE VALUES*/

	int i ;
	int j ;
	double pr_max ;
	double pr_sum ;
	double un ;

	/*get max*/
	pr_max = (-fx[2]-fx[1]) ;
	for (i=2; i<=N_population; i++)
		for (j=1; j<=i-1; j++)
			if ( (-fx[i]-fx[j]) > pr_max ) pr_max = (-fx[i]-fx[j]) ;

	/*get sum*/
	pr_sum = 0.0 ;
	for (i=2; i<=N_population; i++)
		for (j=1; j<=i-1; j++)
			pr_sum += exp( (-fx[i]-fx[j])/temp -pr_max/temp ) ;

	/*get the probability*/
	if (pr_sum!=0)
		*prob = exp( (-fx[n1]-fx[n2])/temp -pr_max/temp )/pr_sum ;
	else
		*prob = log(2.0) -log(N_population) -log(N_population-1.0) ;

}

void Crossover_snooker(double **x, double *fx,
		int N_dimension, int N_population,
		double *theta, double *grid_points, int grid_size,
		double temp, double scl, double *accpr_pop,
		double *z_new) {

	double temp_co ;

	int n1 ;
	int n2 ;
	int i ;
	int k_old ;
	int k_new ;
	double pr_fw ;
	double fz_old ;
	double fz_new ;
	double sumv ;

	double un ;

	double rat ;

	/*PARAMETERS*/

	temp_co = 10.1 ;

	/*CHOOSE INDIVIDUALS*/

	CO_select_forward_1(&pr_fw, &n1, &n2, fx, N_population, temp_co) ;

	/*FORWARD*/

	fz_old = fx[n1] ;

	self_adj_index_search(&k_old, fz_old, grid_points, grid_size) ;

	/* PROPOSE */
    sumv=0.0 ;
	for ( i=1 ; i<=N_dimension ; i++ ){
		z_new[i] = x[n2][i]-x[n1][i] ;
		sumv += z_new[i]*z_new[i] ;
	}
	for ( i=1 ; i<=N_dimension ; i++ ) z_new[i] /= sqrt(sumv) ;

	un = normalrng()*scl ;
	for ( i=1 ; i<=N_dimension ; i++ ) z_new[i] = x[n1][i] + un*z_new[i] ;

	/*BACKWARD*/

	fz_new = cost( z_new, N_dimension ) ;

    self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

    /*ACCEPT/REJECT*/

    rat = -theta[k_new] -fz_new/temp +theta[k_old] +fz_old/temp  ;

    *accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;
    un = uniformrng() ;

	if ( *accpr_pop>un ){
		fx[n1] = fz_new ;
		for( i=1 ; i<=N_dimension ; i++) x[n1][i] = z_new[i] ;
	}

}

void Crossover_linear(double **x, double *fx,
		int N_dimension, int N_population,
		double *theta, double *grid_points, int grid_size,
		double temp, double scl, double *accpr_pop,
		double *z_new){

	double temp_co ;

	int n1 ;
	int n2 ;
	int i ;
	int k_old ;
	int k_new ;
	double pr_fw ;
	double fz_old ;
	double fz_new ;
	double sumv ;

	double un ;

	double rat ;

	/*PARAMETERS*/

	temp_co = 10.1 ;

	/*CHOOSE INDIVIDUALS*/

	CO_select_forward_1(&pr_fw, &n1, &n2, fx, N_population, temp_co) ;

	/*FORWARD*/

	fz_old = fx[n1] ;

    self_adj_index_search(&k_old, fz_old, grid_points, grid_size) ;

	/* PROPOSE */

	un = (2*uniformrng()-1.0)*scl ;
	for ( i=1 ; i<=N_dimension ; i++ ) z_new[i] = x[n1][i] + un*x[n2][i] ;

	/*BACKWARD*/

	fz_new = cost( z_new, N_dimension ) ;

    self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

    /*ACCEPT/REJECT*/

    rat = -theta[k_new] -fz_new/temp +theta[k_old] +fz_old/temp  ;

    *accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;
    un = uniformrng() ;

	if ( *accpr_pop>un ){
		fx[n1] = fz_new ;
		for( i=1 ; i<=N_dimension ; i++) x[n1][i] = z_new[i] ;
	}

}

void Crossover_Kpoint(double **x, double *fx,
		int N_dimension, int N_population,
		double *theta, double *grid_points, int grid_size,
		double temp, double *accpr_pop,
		double *z_new_1, double *z_new_2 ) {

	double temp_co ;

	int n1 ;
	int n2 ;
	int i ;
	int K_CO ;
	int lab_1 ;
	int lab_2 ;
	int k_old_1 ;
	int k_old_2 ;
	int k_new_1 ;
	int k_new_2 ;
	double pr_fw ;
	double pr_bw ;
	double fz_old_1 ;
	double fz_old_2 ;
	double fz_new_1 ;
	double fz_new_2 ;
	double sumv ;

	double un ;

	double rat ;

	/*PARAMETERS*/

	temp_co = 10.1 ;

	/*CHOOSE INDIVIDUALS*/

	CO_select_forward_3(&pr_fw, &n1, &n2, fx, N_population, temp_co) ;

	/*FORWARD*/

	fz_old_1 = fx[n1] ;

	self_adj_index_search(&k_old_1, fz_old_1, grid_points, grid_size) ;

	fz_old_2 = fx[n2] ;

	self_adj_index_search(&k_old_2, fz_old_2, grid_points, grid_size) ;

	/* PROPOSE */

    K_CO = integerrng(1,2) ;
	if ( K_CO == 1 ){

		/* ONE POINT CO */

		lab_1 = integerrng(1,N_dimension-1) ;
		for ( i=1 ; i<=lab_1 ; i++ ){
			z_new_1[i] = x[n1][i] ;
			z_new_2[i] = x[n2][i] ;
		}
		for ( i=lab_1+1 ; i<=N_dimension; i++){
			z_new_1[i] = x[n2][i] ;
			z_new_2[i] = x[n1][i] ;
		}
	}
	else if ( K_CO == 2 ){

		/* TWO POINT CO */

		lab_1 = integerrng(1,N_dimension-1) ;
		lab_2 = integerrng(1,N_dimension-2) ;
		if (lab_2 >= lab_1) lab_2++ ;
		if (lab_1>lab_2) {
			i = lab_1; lab_1=lab_2 ; lab_2 =i ;
		}

		for ( i=1 ; i<=lab_1 ; i++ ){
			z_new_1[i] = x[n1][i] ;
			z_new_2[i] = x[n2][i] ;
		}

		for ( i=lab_1+1 ; i<=lab_2 ; i++ ){
			z_new_1[i] = x[n2][i] ;
			z_new_2[i] = x[n1][i] ;
		}

		for ( i=lab_2+1 ; i<=N_dimension ; i++ ){
			z_new_1[i] = x[n1][i] ;
			z_new_2[i] = x[n2][i] ;
		}
	}

	/*BACKWARD*/

	fz_new_1 = cost( z_new_1, N_dimension ) ;

	self_adj_index_search(&k_new_1, fz_new_1, grid_points, grid_size) ;

	fz_new_2 = cost( z_new_2, N_dimension ) ;

	self_adj_index_search(&k_new_2, fz_new_2, grid_points, grid_size) ;

	CO_select_backward_3(&pr_bw, n1, n2, fx, N_population, temp_co) ;

    /*ACCEPT/REJECT*/

	rat = -theta[k_new_1] -fz_new_1/temp
			-theta[k_new_2] -fz_new_2/temp
			+pr_bw
			+theta[k_old_1] +fz_old_1/temp
			+theta[k_old_2] +fz_old_2/temp
			-pr_fw ;

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;
    un = uniformrng() ;

	if ( *accpr_pop>un ){
		fx[n1] = fz_new_1 ;
		for( i=1 ; i<=N_dimension ; i++) x[n1][i] = z_new_1[i] ;
		fx[n2] = fz_new_2 ;
		for( i=1 ; i<=N_dimension ; i++) x[n2][i] = z_new_2[i] ;
	}

}












