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


void CO_select_forward_0(double *, int *, int *, int ) ;

void CO_select_backward_0(double *, int ) ;

void CO_select_forward_1(double *, int *, int *, double *, int , double ) ;

void CO_select_backward_1(double *, int , int , double *, int , double ) ;

void CO_select_forward_2(double *, int *, int *, double *, int , double ) ;

void CO_select_backward_2(double *, int , int , double *, int , double ) ;

void CO_select_forward_3(double *, int *, int *, double *, int , double ) ;

void CO_select_backward_3(double *, int , int , double *, int , double ) ;

void Crossover_snooker(double **, double *,
					int , int ,
					double *, double *, int ,
					double , double , double *,
					double *) ;

void Crossover_linear(double **, double *,
					int , int ,
					double *, double *, int ,
					double , double , double *,
					double *) ;

void Crossover_Kpoint(double **, double *,
					int , int ,
					double *, double *, int ,
					double , double *,
					double *, double * ) ;














