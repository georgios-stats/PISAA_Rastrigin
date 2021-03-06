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


/*
 * References:
 *
 * P. Muller (1991),
 * A Generic Approach to Posterior Integration and Gibbs Sampling},
 * Purdue University, Department of Statistics
 * Indiana
}*/


#include <math.h>
#include "nrutil.h"

#include "RNG.h"

/* MATRIX MULTIPLICATION: C(N,M)=A(N,P)*B(P,M) */
static void MulMatrix( double **A, double **B, double **C, int n, int p, int m){

	int i ;
	int j ;
	int k ;

	for (i=1; i<=n; i++)
	   for (j=1; j<=m; j++){
			   C[i][j] = 0.0 ;
			   for (k=1; k<=p; k++) C[i][j] += A[i][k]*B[k][j] ;
	   }

}

static void MRot( double **MTurn, int i, int j, int N, double alpha){

	int k1 ;
	int k2 ;

	for (k1=1; k1<=N; k1++)
		for (k2=1; k2<=N; k2++)
			if (k1==i && k2==i)
				MTurn[k1][k2] = cos( alpha ) ;
			else if (k1==j && k2==j)
				MTurn[k1][k2] = cos( alpha ) ;
			else if (k1==i && k2==j)
				MTurn[k1][k2] = sin( alpha ) ;
			else if (k1==j && k2==i)
				MTurn[k1][k2] = -sin( alpha ) ;
			else
				MTurn[k1][k2] = ((k1==k2)?1.0:0.0) ;
}

void CreateRotationMatrix( double **MResult, int N )
{

	int i ;
	int k1 ;
	int k2 ;
	double **MTurn ;
	double **MWork ;
	double alpha ;
	double myPI ;

	myPI = 3.141592653589793;

	MTurn = dmatrix( 1, N, 1, N ) ;
	MWork = dmatrix( 1, N, 1, N ) ;

	/*Initialise MResult*/
	for (k1=1; k1<=N; k1++)
		for (k2=1; k2<=N; k2++)
			MResult[k1][k2] = ((k1==k2)?1.0:0.0) ;

	for (i=2; i<=N; i++){

		alpha = ( uniformrng()-0.5 ) * myPI * 0.5 ;
		MRot( MTurn, 1, i, N, alpha) ;

		for (k1=1; k1<=N; k1++)
			for (k2=1; k2<=N; k2++)
				MWork[k1][k2] = MResult[k1][k2] ;

		MulMatrix( MWork, MTurn, MResult, N, N, N ) ;
	}

	for ( i=2 ; i<=N-1 ; i++ ){

		alpha = ( uniformrng()-0.5 ) * myPI * 0.5 ;
		MRot( MTurn, i, N, N, alpha) ;

		for (k1=1; k1<=N; k1++)
			for (k2=1; k2<=N; k2++)
				MWork[k1][k2] = MResult[k1][k2] ;

		MulMatrix( MWork, MTurn, MResult, N, N, N ) ;

		}

	free_dmatrix( MTurn, 1, N, 1, N ) ;
	free_dmatrix( MWork, 1, N, 1, N ) ;

}


