/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reaxc.h"
#include "reaxc_bonds.h"
#include "reaxc_bond_orders.h"
#include "reaxc_list.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"
#include <iostream>
#include <fstream>

void Bonds( reax_system *system, control_params * /*control*/,
            simulation_data *data, storage *workspace, reax_list **lists,
            output_controls * /*out_control*/ )
{
  int i, j, pj, natoms;
  int start_i, end_i;
  int type_i, type_j;
  double ebond, pow_BOs_be2, exp_be12, CEbo;
  double gp3, gp4, gp7, gp10, gp37;
  double exphu, exphua1, exphub1, exphuov, hulpov, estriph;
  double decobdbo, decobdboua, decobdboub;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij;
  reax_list *bonds;
  /* added code by ETHAN ---------------------------------------------------------------------*/
  double F1rest, F2rest, R12rest, Erest;
  int id1rest, id2rest;
  double Rijrest, dxrest, dyrest, dzrest;
  int nRowsrest;
  double** restMatrix; // pointer to pointer (dynamic memory allocation); is deallocated at end of file
  std::ifstream restfile ("rest-data.txt");
  if (restfile.is_open()) {
	  restfile >> nRowsrest;		//grabs number of atoms, listed by python
	  // allocate
	  restMatrix = new double*[nRowsrest];		//creates an array of arrays on the heap, so it does not go out of scope
	  for (i = 0; i < nRowsrest; i++) {
		  restMatrix[i] = new double[5];		//5? this is very odd.  I guess he only cares about the first five things in each row, which makes sense
	  }
	  // fill 2d array with infile parameters
	  for (i = 0; i < nRowsrest; i++) {
		  for (j = 0; j < 5; j++) {
			  restfile >> restMatrix[i][j];
		  }
	  }
	  restfile.close();
  }
  else { 
    std::cout << "Unable to open restraint parameter file" << std::endl;
	restMatrix = new double*[1]; // wouldn't be used in this case anyway?
  }
  /* END OF ADDED CODE BY ETHAN */

  bonds = (*lists) + BONDS;
  gp3 = system->reax_param.gp.l[3];
  gp4 = system->reax_param.gp.l[4];
  gp7 = system->reax_param.gp.l[7];
  gp10 = system->reax_param.gp.l[10];
  gp37 = (int) system->reax_param.gp.l[37];
  natoms = system->n;

  for( i = 0; i < natoms; ++i ) {
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);

    for( pj = start_i; pj < end_i; ++pj ) {
      j = bonds->select.bond_list[pj].nbr;

      if( system->my_atoms[i].orig_id > system->my_atoms[j].orig_id )
        continue;
      if( system->my_atoms[i].orig_id == system->my_atoms[j].orig_id ) {
        if (system->my_atoms[j].x[2] <  system->my_atoms[i].x[2]) continue;
        if (system->my_atoms[j].x[2] == system->my_atoms[i].x[2] &&
            system->my_atoms[j].x[1] <  system->my_atoms[i].x[1]) continue;
        if (system->my_atoms[j].x[2] == system->my_atoms[i].x[2] &&
            system->my_atoms[j].x[1] == system->my_atoms[i].x[1] &&
            system->my_atoms[j].x[0] <  system->my_atoms[i].x[0]) continue;
      }

      /* set the pointers */
      type_i = system->my_atoms[i].type;
      type_j = system->my_atoms[j].type;
      sbp_i = &( system->reax_param.sbp[type_i] );
      sbp_j = &( system->reax_param.sbp[type_j] );
      twbp = &( system->reax_param.tbp[type_i][type_j] );
      bo_ij = &( bonds->select.bond_list[pj].bo_data );

      /* calculate the constants */
      if (bo_ij->BO_s == 0.0) pow_BOs_be2 = 0.0;
      else pow_BOs_be2 = pow( bo_ij->BO_s, twbp->p_be2 );
      exp_be12 = exp( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
      CEbo = -twbp->De_s * exp_be12 *
        ( 1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2 );

      /* calculate the Bond Energy */
      data->my_en.e_bond += ebond =
        -twbp->De_s * bo_ij->BO_s * exp_be12
        -twbp->De_p * bo_ij->BO_pi
        -twbp->De_pp * bo_ij->BO_pi2; /*-*2--------PREVIOUSLY MODIFIED BY ETHAN--------*/

      /* ADD A RESTRAINT TERM TO THE BOND ENERGY */
      /* added code by ETHAN ---------------------------------------------------------------------*/
      // Comments: i = atom1; j = atom2 but loop only runs if j > i so a bond is never double counted
      //
      Erest = 0;
      for (int restk = 0; restk < nRowsrest; restk++) {
        id1rest = restMatrix[restk][0];
		id2rest = restMatrix[restk][1];
		if ( ((id1rest == system->my_atoms[i].orig_id) and (id2rest == system->my_atoms[j].orig_id)) or ((id1rest == system->my_atoms[j].orig_id) and (id2rest == system->my_atoms[i].orig_id)) ){
		  R12rest = restMatrix[restk][2];
		  F1rest = restMatrix[restk][3];
		  F2rest = restMatrix[restk][4];
          dxrest = system->my_atoms[i].x[0] - system->my_atoms[j].x[0];
          dyrest = system->my_atoms[i].x[1] - system->my_atoms[j].x[1];
          dzrest = system->my_atoms[i].x[2] - system->my_atoms[j].x[2];
          Rijrest = sqrt(dxrest*dxrest + dyrest*dyrest + dzrest*dzrest);
		  Erest = F1rest*(1.0 - exp(F2rest*(Rijrest - R12rest)*(Rijrest - R12rest) ));
		}		  
	  }
      ebond = ebond + Erest;
      data->my_en.e_bond = data->my_en.e_bond + Erest;
      /* end added code --------------------------------------------------------------------------*/

      /* tally into per-atom energy */
      if( system->pair_ptr->evflag)
        system->pair_ptr->ev_tally(i,j,natoms,1,ebond,0.0,0.0,0.0,0.0,0.0);

      /* calculate derivatives of Bond Orders */
      bo_ij->Cdbo += CEbo;
      bo_ij->Cdbopi -= (CEbo + twbp->De_p);
      bo_ij->Cdbopi2 -= (CEbo + twbp->De_pp);
      // ADDED CODE BY ETHAN -----------------
	  // "Cdbo, Cdbopi, Cdbopi2" are used to calculate forces in "reaxc_bond_orders.cpp" line 50
	  // need to include bond order derivative coefficient for restraint energy in order for forces to be updated
	  
	  
	  double restCdbo, restCdp, restCdpp; // twbp->r_s , twbp->r_pi , twbp->r_pi_pi
	  double term1, term2, term3, term4;
	  if (Erest != 0) {  // use boi or bei ?
		  term1 = pow(log(bo_ij->BO_s)/twbp->p_bo1, 1.0/twbp->p_bo2);
	      term2 = twbp->r_s*term1-R12rest;
	      term3 = exp(F2rest*pow(term2,2));
	      term4 = -2.0*F1rest*F2rest*twbp->r_s/twbp->p_bo2/bo_ij->BO_s/log(bo_ij->BO_s);
	      restCdbo = term1*term2*term3*term4; // sigmabond restraint coeff
		  term1 = pow(log(bo_ij->BO_pi)/twbp->p_bo3, 1.0/twbp->p_bo4);
	      term2 = twbp->r_p*term1-R12rest;
	      term3 = exp(F2rest*pow(term2,2));
	      term4 = -2.0*F1rest*F2rest*twbp->r_p/twbp->p_bo4/bo_ij->BO_pi/log(bo_ij->BO_pi);
	      restCdp = term1*term2*term3*term4; // pibond restraint coeff
		  // std::cout << bo_ij->BO_pi2 <<" " <<id1rest <<" " <<id2rest <<Rijrest << std::endl;  // I had to modify "reaxc_bond_orders.cpp" to remove e-10 = 0 condition
		  term1 = pow(log(bo_ij->BO_pi2)/twbp->p_bo5, 1.0/twbp->p_bo6);
	      term2 = twbp->r_pp*term1-R12rest;
	      term3 = exp(F2rest*pow(term2,2));
	      term4 = -2.0*F1rest*F2rest*twbp->r_pp/twbp->p_bo6/bo_ij->BO_pi2/log(bo_ij->BO_pi2);
	      restCdpp = term1*term2*term3*term4; //*term2*term3*term4; // pipibond restraint coeff
		  bo_ij->Cdbo += restCdbo;
	      bo_ij->Cdbopi += restCdp;
	      bo_ij->Cdbopi2 += restCdpp;
		  // std::cout << bo_ij->BO_s << " " << bo_ij->BO_pi << " " << bo_ij->BO_pi2 << std::endl;
	  }	
      	  
	  // END ADDED CODE BY ETHAN --------------

      /* Stabilisation terminal triple bond */
      if( bo_ij->BO >= 1.00 ) {
        if( gp37 == 2 ||
            (sbp_i->mass == 12.0000 && sbp_j->mass == 15.9990) ||
            (sbp_j->mass == 12.0000 && sbp_i->mass == 15.9990) ) {
          exphu = exp( -gp7 * SQR(bo_ij->BO - 2.50) );
          exphua1 = exp(-gp3 * (workspace->total_bond_order[i]-bo_ij->BO));
          exphub1 = exp(-gp3 * (workspace->total_bond_order[j]-bo_ij->BO));
          exphuov = exp(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
          hulpov = 1.0 / (1.0 + 25.0 * exphuov);

          estriph = gp10 * exphu * hulpov * (exphua1 + exphub1);
          data->my_en.e_bond += estriph;

          decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1) *
            ( gp3 - 2.0 * gp7 * (bo_ij->BO-2.50) );
          decobdboua = -gp10 * exphu * hulpov *
            (gp3*exphua1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));
          decobdboub = -gp10 * exphu * hulpov *
            (gp3*exphub1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));

          /* tally into per-atom energy */
          if( system->pair_ptr->evflag)
            system->pair_ptr->ev_tally(i,j,natoms,1,estriph,0.0,0.0,0.0,0.0,0.0);

          bo_ij->Cdbo += decobdbo;
          workspace->CdDelta[i] += decobdboua;
          workspace->CdDelta[j] += decobdboub;
        }
      }
    }
  }
  // delete the dynamically allocated 2d array 'restMatrix'
  for (i = 0; i < nRowsrest; i++){
	  delete[] restMatrix[i];
  }
  delete[] restMatrix;
}
