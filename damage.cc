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


void damage( long int gr, double new_epe[], double new_sig[], 
  double old_damage, double &new_damage )

{
  long int swit=0, ldum=0, idum[1];
  double damage_data[DATA_ITEM_SIZE];

  swit = set_swit(-1,-1,"damage");
  if ( swit ) pri( "In routine DAMAGE" );

  if ( db( GROUP_MATERI_DAMAGE_MAZARS, gr, idum, damage_data, 
    ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
    damage_mazars( damage_data, new_epe, new_sig, old_damage, new_damage );
  }
  else
    new_damage = old_damage;

  if ( swit ) pri( "Out routine DAMAGE" );

}

void damage_mazars( double materi_damage_mazars[], double new_epe[], double new_sig[], 
  double old_damage, double &new_damage )
 
{
  long int idim=0, swit=0, iwork[MDIM];
  double alpha=0., gamma=0., epseq=0., at=0., bt=0., ac=0., bc=0., beta=0.,
    dt=0., dc=0., tmp1=0., tmp2=0., eps0=0., epssiz=0., epsp[MDIM], 
    work1[MDIM*MDIM], work2[MDIM*MDIM];

  swit = set_swit(-1,-1,"damage_mazars");
  if ( swit ) pri( "In routine DAMAGE_MAZARS" );
 
  eps0 = materi_damage_mazars[0];
  at = materi_damage_mazars[1];
  bt = materi_damage_mazars[2];
  ac = materi_damage_mazars[3];
  bc = materi_damage_mazars[4];
  beta = materi_damage_mazars[5];
 
  tmp1 = array_size( new_epe, MDIM*MDIM );
  if ( tmp1<1.e-10 ) {
    alpha = 1.;
    epseq = 0.;
    dc = 0.;
    dt = 0.;
  }
  else {
      /* gamma */
    tmp1 = 0.;
    tmp2 = 0.;
    for ( idim=0; idim<MDIM; idim++ ) {
      if ( new_sig[idim]<0. ) {
        tmp1 += new_sig[idim]*new_sig[idim];
        tmp2 += new_sig[idim];
      }
    }
    if ( scalar_dabs(tmp2)<=1.e-6 )
      gamma = 1.;
    else
      gamma = -sqrt(scalar_dabs(tmp1))/tmp2;
        /* epqeq */
    array_move( new_epe, work1, MDIM*MDIM );
    matrix_jacobi( work1, MDIM, epsp, work2, iwork );
    epssiz = 0.;
    epseq = 0.;
    for ( idim=0; idim<MDIM; idim++ ) {
      epssiz += epsp[idim]*epsp[idim];
      if ( epsp[idim]>0. ) epseq += epseq + epsp[idim]*epsp[idim];
    }
    epssiz = sqrt( scalar_dabs(epssiz) );
    epseq = sqrt( scalar_dabs(epseq) );
    epseq = gamma * epseq;
      /* dt and dc */
    if ( epseq>eps0 ) {
      alpha = (epseq*epseq)/(epssiz*epssiz);
      if ( alpha<0. ) alpha = 0.;
      if ( alpha>1. ) alpha = 1.;
      dt = 1. - (eps0/epseq)*(1.-at) - at*exp(-bt*(epseq-eps0));
      dc = 1. - (eps0/epseq)*(1.-ac) - ac*exp(-bc*(epseq-eps0));
      if ( dt<0. ) dt = 0.;
      if ( dc<0. ) dc = 0.;
      if ( dt>1. ) dt = 1.;
      if ( dc>1. ) dc = 1.;
    }
    else {
      alpha = 0.;
      dt = 0.;
      dc = 0.;
    }
  }

    /* d */
  new_damage = dt*scalar_power(alpha,beta) + dc*scalar_power(1.-alpha,beta);
  if ( new_damage<old_damage ) new_damage = old_damage;

  if ( swit ) pri( "Out routine DAMAGE_MAZARS" );

}
