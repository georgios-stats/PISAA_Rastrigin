/*
 * Copyrigtht 2014 Georgios Karagiannis
 *
 * This file is part of PISAA_Rastrigin.
 *
 * PISAA_Rastrigin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * PISAA_Rastrigin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISAA_Rastrigin.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Georgios Karagiannis 
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 (765) 496-1007
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
*/



#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#ifndef INFINITY
	#include <float.h>
	#define INFINITY DBL_MAX
#endif

#ifndef __ROTATE__
	#define __ROTATE__ 1
#endif

#if __ROTATE__
	void CreateRotationMatrix(double**,int) ;

	static double **data=NULL ;

#endif

/*

MULTIMODAL

Rastrigin(n)
minimize

subject to -2.56<=xi<=5.12 for i = 1; . . . ; n:
subject to -5.12<=xi<=5.12 for i = 1; . . . ; n:

Description:
Dimensions: d

The Rastrigin function has several local minima. It is highly multimodal, but locations of the minima are regularly distributed. It is shown in the plot above in its two-dimensional form.

Input Domain:
The function is usually evaluated on the hypercube xi ∈ [-5.12, 5.12], for all i = 1, …, d.

References:

Global Optimization Test Problems. Retrieved June 2013, from
http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.

Pohlheim, H. GEATbx Examples: Examples of Objective Functions (2005).
Retrieved June 2013, from http://www.geatbx.com/download/GEATbx_ObjFunExpl_v37.pdf.

*/


void get_data(int N_dimension){
#if __ROTATE__
	data = dmatrix(1,N_dimension,1,N_dimension) ;
	CreateRotationMatrix( data, N_dimension ) ;
#else
	;
#endif
}

void cost_bounds(double *z_min, double *z_max, int i){
	if (i > 0){
		*z_min = -5.12 ;
		*z_max = 5.12 ;
	}
}

double cost(double *y, int N_dimension){

	int i, j ;
	double ff ;
	double my_pi ;

	double *z ;

	double y_min ;
	double y_max ;

	my_pi = 3.14159265358979323846 ;

	for ( i = 1 ; i <= N_dimension ; i++){
		cost_bounds(&y_min, &y_max, i) ;
		if (y[i]>y_max || y[i]<y_min) return (double) INFINITY ;
	}

#if __ROTATE__
	z = dvector(1,N_dimension) ;
	for (i =1;i<= N_dimension;i++){
			z[i]=0.0;
			for (j =1;j<= N_dimension;j++)
				z[i] += data[i][j]*y[j] ;
	}

/*	for (i=1;i<=N_dimension;i++) {
		for (j=1;j<=N_dimension;j++)
			printf("%f ",data[i][j]) ;
			printf("\n") ;
	}
	printf("%d\n",dataQ) ;*/


#else
	z = y-1 ;
#endif

	ff = 10.0 * N_dimension ;
	for ( i = 1 ; i<= N_dimension ; i++ )
		ff += (z[i]*z[i] -10.0*cos(2.0*my_pi*z[i])) ;

#if __ROTATE__
	free_dvector(z,1,N_dimension) ;
#endif

	return ff ;

}

