/*
    Copyright (C) 1998  Dennis Roddeman
    email: dennis.roddeman@feat.nl

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation 
    59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA
*/

#include "tochnog.h"

void parallel_new_dof_before( void )

{
  long int max_node=0, ipuknwn=0, iuknwn=0, idim=0, iloop=0, nloop=0, 
    inod=0, ldum=0, ithread=0, idum[1], *next_of_loop=NULL;
  double dtime=0., tmp=0., force_gravity[MDIM], node_damping[MDIM], 
    node_stiffness[MDIM], node_mass[MDIM], *node_dof=NULL, *node_dof_new=NULL, 
    *node_rhside=NULL, *node_lhside=NULL;

  if ( materi_velocity ) {
    array_set( node_damping, 0., ndim );
    array_set( node_stiffness, 0., ndim );
    array_set( node_mass, 0., ndim );
    array_set( force_gravity, 0., ndim );
    db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
    force_gravity_calculate( force_gravity );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    if ( max_node>=0 ) {
      next_of_loop = get_new_int(1+max_node);
      parallel_sys_next_of_loop( next_of_loop, max_node, nloop, ithread );
      for ( iloop=0; iloop<nloop; iloop++ ) {
        inod = next_of_loop[iloop];
        if ( inod>max_node )
          break;
        else if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
          node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
          node_dof_new = db_dbl( NODE_DOF, inod, VERSION_NEW );
          node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
          node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
          db( NODE_DAMPING, inod, idum, node_damping, ldum, 
            VERSION_NORMAL, GET_IF_EXISTS );
          db( NODE_STIFFNESS, inod, idum, node_stiffness, ldum, 
            VERSION_NORMAL, GET_IF_EXISTS );
          db( NODE_MASS, inod, idum, node_mass, ldum, 
            VERSION_NORMAL, GET_IF_EXISTS );
          for ( idim=0; idim<ndim; idim++ ) {
            ipuknwn = vel_indx/nder + idim;
            iuknwn = vel_indx + idim * nder;
            node_lhside[ipuknwn] += node_stiffness[idim]*dtime + node_damping[idim] + 
              node_mass[idim] / dtime;
            tmp = node_damping[idim] * node_dof_new[iuknwn] +
              node_mass[idim] * ( node_dof_new[iuknwn] - node_dof[iuknwn] ) / 
              dtime + node_mass[idim] * force_gravity[idim];
            if ( materi_displacement ) 
              tmp += node_stiffness[idim] * node_dof_new[dis_indx+idim*nder];
            else if ( materi_velocity_integrated ) 
              tmp += node_stiffness[idim] * node_dof_new[veli_indx+idim*nder];
            node_rhside[ipuknwn] -= tmp;
          }
        }
      }
      delete[] next_of_loop;
    }
  }

}

void parallel_new_dof_diagonal( void )

{
  long int inod=0, max_node=0, ipuknwn=0, iuknwn=0, iprinc=0, iloop=0, nloop=0, 
    ithread=0, icontrol=0, idim=0, ind=0, options_solver=-MATRIX_ITERATIVE_BICG, 
    ldum=0, idum[1], dof_principal[MUKNWN], *next_of_loop=NULL, *node_bounded=NULL;
  double tmp=0., dtime=0., ddum[1], options_relaxation[MPRINC], *node_dof=NULL, 
    *node_dof_new=NULL, *node_rhside=NULL, *node_lhside=NULL;

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET );
  db( OPTIONS_SOLVER, 0, &options_solver, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_SOLVER, icontrol, &options_solver, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_PRINCIPAL, 0, dof_principal, ddum, ldum, VERSION_NORMAL, GET );
  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  array_set( options_relaxation, 1., nprinc );
  db( OPTIONS_RELAXATION, 0, idum, options_relaxation, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_RELAXATION, icontrol, idum, options_relaxation, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS );

  if ( max_node>=0 ) {
    next_of_loop = get_new_int(1+max_node);
    parallel_sys_next_of_loop( next_of_loop, max_node, nloop, ithread );
    for ( iloop=0; iloop<nloop; iloop++ ) {
      inod = next_of_loop[iloop];
      if ( inod>max_node )
        break;
      else if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_NORMAL );
        node_dof_new = db_dbl( NODE_DOF, inod, VERSION_NEW );
        node_bounded = db_int( NODE_BOUNDED, inod, VERSION_NORMAL );
        node_lhside = db_dbl( NODE_LHSIDE, inod, VERSION_NORMAL );
        node_rhside = db_dbl( NODE_RHSIDE, inod, VERSION_NORMAL );
        iprinc = 0;
        for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
          iuknwn = ipuknwn*nder;
          if ( !node_bounded[ipuknwn] && node_lhside[ipuknwn]!=0. ) {
            if ( options_solver==-DIAGONAL || dof_principal[iuknwn]<0 ) {
              tmp = node_rhside[ipuknwn] / node_lhside[ipuknwn];
              if ( dof_principal[iuknwn]>=0 ) tmp *= options_relaxation[iprinc];
              node_dof_new[iuknwn] += tmp;
            }
          }
          if ( dof_principal[iuknwn]>=0 ) iprinc++;
        }
        if ( materi_damage ) {
          iuknwn = dam_indx;
          if ( node_dof_new[iuknwn]<0. ) node_dof_new[iuknwn] = 0.;
          if ( node_dof_new[iuknwn]>1. ) node_dof_new[iuknwn] = 1.;
        }
        if ( materi_displacement ) {
          for ( idim=0; idim<ndim; idim++ ) {
            iuknwn = dis_indx + idim * nder;
            ind = vel_indx+idim*nder;   
            node_dof_new[iuknwn] = node_dof[iuknwn] + node_dof_new[ind] * dtime;   
          }
        }
        if ( materi_plasti_kappa ) {
          iuknwn = kap_indx;
          if ( node_dof_new[iuknwn]<0. ) node_dof_new[iuknwn] = 0.;
        }
        if ( materi_velocity_integrated ) {
          for ( idim=0; idim<ndim; idim++ ) {
            iuknwn = veli_indx + idim * nder;
            ind = vel_indx+idim*nder;   
            node_dof_new[iuknwn] = node_dof[iuknwn] + node_dof_new[ind] * dtime;   
          }
        }
        if ( materi_void_fraction ) {
          iuknwn = void_indx;
          if ( node_dof_new[iuknwn]<0. ) node_dof_new[iuknwn] = 0.;
          if ( node_dof_new[iuknwn]>1. ) node_dof_new[iuknwn] = 1.;
        }
        if ( maxwell_e ) {
          for ( idim=0; idim<MDIM; idim++ ) {
            iuknwn = maxe_indx + idim*nder;
            node_dof_new[iuknwn] = node_dof[iuknwn] +
              node_dof_new[maxfe_indx+idim*nder] * dtime;   
          }
        }
        if ( wave_scalar ) {
          iuknwn = scal_indx;
          node_dof_new[iuknwn] = node_dof[iuknwn] +
            node_dof_new[fscal_indx] * dtime;   
        }
        if ( derivatives ) {
          for ( ipuknwn=0; ipuknwn<npuknwn; ipuknwn++ ) {
            if ( !node_bounded[ipuknwn] ) {
              iuknwn = ipuknwn*nder;
              node_dof_new[iuknwn+nder-1] = 
                ( node_dof_new[iuknwn] - node_dof[iuknwn] ) / dtime;
            }
          }
        }
        if ( materi_velocity ) {
          for ( idim=0; idim<ndim; idim++ ) {
            iuknwn = vel_indx + idim * nder;
          }
        }
      }
    }
    delete[] next_of_loop;
  }

}
