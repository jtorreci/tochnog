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

#define EPS_DELETE 1.e-12
#define EPS_ELEMENT_DELETE_FACTOR 1.e-6

void delete_element( long int element, long int version )

{
  long int idat=0, data_class=0, swit=0;

  swit = set_swit(-1,-1,"delete_element");
  if ( swit ) pri( "In routine DELETE_ELEMENT" );

  for ( idat=0; idat<MDAT; idat++ ) {
    data_class = db_data_class( idat );
    if ( data_class==ELEMENT && db_version(idat,version) ) 
      db_delete_index( idat, element, version );
  }

  if ( swit ) pri( "Out routine DELETE_ELEMENT" );

}

void delete_node( long int inod, long int version )

{
  long int idat=0, data_class=0, swit=0;

  swit = set_swit(-1,-1,"delete_node");
  if ( swit ) pri( "In routine DELETE_NODE" );

  for ( idat=0; idat<MDAT; idat++ ) {
    data_class = db_data_class( idat );
    if ( data_class==NODE && db_version(idat,version) ) 
      db_delete_index( idat, inod, version );
  }

  if ( swit ) pri( "Out routine DELETE_NODE" );

}

void delete_geom( double time_current )

{
  long int element=0, max_element=0, inol=0, nnol=0, inod=0, length=0,
    all_in_geometry=0, in_geometry=0, element_group=0,
    max_node=0, any_element_deleted=0,
    node_boundary=-YES, swit=0, icontrol=0, 
    control_mesh_delete_geometry_movenodes=-YES, 
    length_control_mesh_delete_geometry_elementgroup=0,
    length_control_mesh_delete_geometry_element=0,
    name=0, ldum=0, idum[1], control_mesh_delete_geometry[2], 
    control_mesh_delete_geometry_elementgroup[DATA_ITEM_SIZE],
    control_mesh_delete_geometry_element[DATA_ITEM_SIZE], 
    nodes[MNOL], el[MNOL+1], *node_in_geometry=NULL;
  double element_volume=0., old_element_volume=0., 
    element_delete_factor=0., time_old=0., 
    time_new=0., rdum=0., f0=0., f1=0., 
    ddum[MDIM], coord[MDIM], diff_coord[MDIM], 
    node_start_refined[MDIM], projection[MDIM],
    control_mesh_delete_geometry_factor[2]; 

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_DELETE_GEOMETRY, 
      icontrol, VERSION_NORMAL )  ) {

    swit = set_swit(-1,-1,"delete_geom");
    if ( swit ) pri( "In routine DELETE_GEOM" );

    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    length = 1+max_node;
    node_in_geometry = get_new_int( length );
    array_set( node_in_geometry, 0, length );

    db( CONTROL_MESH_DELETE_GEOMETRY, icontrol, control_mesh_delete_geometry, 
      ddum, ldum, VERSION_NORMAL, GET );

    control_mesh_delete_geometry_element[0] = -ALL;
    length_control_mesh_delete_geometry_element = 1;
    db( CONTROL_MESH_DELETE_GEOMETRY_ELEMENT, icontrol, 
      control_mesh_delete_geometry_element, ddum, 
      length_control_mesh_delete_geometry_element, VERSION_NORMAL,
      GET_IF_EXISTS );

    control_mesh_delete_geometry_factor[0] = 0.;
    control_mesh_delete_geometry_factor[1] = 1.;
    db( CONTROL_MESH_DELETE_GEOMETRY_FACTOR, icontrol, idum,
      control_mesh_delete_geometry_factor, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    f0 = 1. - control_mesh_delete_geometry_factor[0];
    f1 = 1. - control_mesh_delete_geometry_factor[1];
    if ( swit ) {
      pri( "f0", f0 );
      pri( "f1", f1 );
    }

    control_mesh_delete_geometry_elementgroup[0] = -ALL;
    length_control_mesh_delete_geometry_elementgroup = 1;
    db( CONTROL_MESH_DELETE_GEOMETRY_ELEMENTGROUP, icontrol, 
      control_mesh_delete_geometry_elementgroup, ddum, 
      length_control_mesh_delete_geometry_elementgroup, VERSION_NORMAL,
      GET_IF_EXISTS );

    db( CONTROL_MESH_DELETE_GEOMETRY_MOVENODES, icontrol, 
      &control_mesh_delete_geometry_movenodes, ddum, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );

    db_version_copy( VERSION_NORMAL, VERSION_TMP );
    db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
    if ( db_active_index( TIME_OLD, 0, VERSION_NORMAL ) ) {
      db( TIME_OLD, 0, idum, &time_old, ldum, VERSION_NORMAL, GET );
      db( TIME_NEW, 0, idum, &time_new, ldum, VERSION_NORMAL, GET );
      if ( scalar_dabs(time_new-time_old)>EPS_DELETE ) {
        element_delete_factor = f0 - ( f0 - f1 ) *
          (time_current-time_old)/(time_new-time_old);
      }
    }

      // determine which nodes are in the geometry
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
         geometry( inod, ddum, control_mesh_delete_geometry, 
           in_geometry, rdum, ddum, rdum, ddum, NODE_START_REFINED, 
           CONTROL_MESH_DELETE_GEOMETRY, VERSION_NORMAL );
         if ( in_geometry ) node_in_geometry[inod] = 1;
      }
    }

      // delete elements which are totally in geometry
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
        db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
        name = el[0];
        element_group = 0;
        db( ELEMENT_GROUP, element, &element_group, ddum, 
          ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( control_mesh_delete_geometry_element[0]==-ALL || 
             array_member(control_mesh_delete_geometry_element,name,
             length_control_mesh_delete_geometry_element,ldum) ) {
          if ( control_mesh_delete_geometry_elementgroup[0]==-ALL || 
               array_member(control_mesh_delete_geometry_elementgroup,element_group,
               length_control_mesh_delete_geometry_elementgroup,ldum) ) {   
            nnol = length - 1; array_move( &el[1], nodes, nnol );
            all_in_geometry = 1;
            for ( inol=0; inol<nnol; inol++ ) {
              inod = nodes[inol];
              if ( !node_in_geometry[inod] ) all_in_geometry = 0;
            }
            if ( all_in_geometry ) {
              if ( element_delete_factor>EPS_ELEMENT_DELETE_FACTOR ) {
                length = 1;
                db( ELEMENT_DELETE_FACTOR, element, idum, &element_delete_factor, 
                  length, VERSION_NORMAL, PUT );
              }
              else {
                delete_element( element, VERSION_NORMAL );
                any_element_deleted = 1;
              }
            }
          }
        }
      }
    }

    if ( any_element_deleted ) {
         // move remaining nodes in geometry to edges
      if ( control_mesh_delete_geometry_movenodes==-YES ) {
        for ( inod=0; inod<=max_node; inod++ ) {
          if ( node_in_geometry[inod] ) {
            geometry( inod, ddum, control_mesh_delete_geometry, 
              in_geometry, rdum, ddum, rdum, projection, 
              NODE_START_REFINED, PROJECT_ON_EDGE, VERSION_NORMAL );
            db( NODE_START_REFINED, inod, idum, node_start_refined, 
              ldum, VERSION_NORMAL, GET );
            db( NODE_START_REFINED, inod, idum, projection, 
              ndim, VERSION_NORMAL, PUT );
            db( NODE, inod, idum, coord, ldum, VERSION_NORMAL, GET );
            array_subtract( projection, node_start_refined, diff_coord, ndim );
            array_add( coord, diff_coord, coord, ndim );
            db( NODE, inod, idum, coord, ndim, VERSION_NORMAL, PUT );
            length=1; db( NODE_BOUNDARY, inod, &node_boundary, ddum, 
              length, VERSION_NORMAL, PUT );
          }
        }
         // delete too small elements (collapsed because of moving nodes to edge)
        for ( element=0; element<=max_element; element++ ) {
          if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
            db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
            name = el[0];
            if ( name==-BAR2 || name==-TRIA3 || name==-TET4 ) {
              nnol = length - 1; array_move( &el[1], nodes, nnol );
              element_volume_set( name, nodes, VERSION_NORMAL, element_volume );
              element_volume_set( name, nodes, VERSION_TMP, old_element_volume );
              if ( element_volume<EPS_VOLUME*old_element_volume ) 
                delete_element( element, VERSION_NORMAL );
            }
          }
        }
      }

  mesh_has_changed( VERSION_NORMAL );
    }
    db_version_delete( VERSION_TMP );

    delete[] node_in_geometry;

    if ( swit ) pri( "Out routine DELETE_GEOM" );

  }

}

// ---------------------------------------------------------------------
// control_mesh_cut_geometry (manual Professional 6.163) +
// control_mesh_cut_node_force (6.164; the corpus/Professional record
// name - the manual writes control_mesh_cut_force, an alias resolved
// in db_number). Semantics measured on the Professional binary
// 25-10-2023 (mesh_cut_1.dat/mesh_cut_2.dat, .dbs comparison):
//   - the cut activates when the control index becomes active (the
//     record index between the timestep blocks, or a timestep at the
//     same index), once;
//   - the mesh part inside the geometry (elements whose nodes ALL lie
//     within the geometry - the same membership test as
//     control_mesh_delete_geometry) is deleted together with the nodes
//     that end up attached to no remaining element (nod_nod pruning);
//   - the removed part is substituted by its nodal forces: the force
//     the removed elements exert on the surviving boundary nodes
//     (Newton's third law: -f_elem, f_elem = int B^T*sigma dV of the
//     removed elements, integrated from the converged integration-point
//     stresses ELEMENT_DOF of the step before the cut) is written to
//     the NODE_FORCE records of those nodes (applied as an external
//     load by parallel_new_dof_before for every later step);
//   - control_mesh_cut_node_force gates the substitution per space
//     direction: only components whose switch is -yes are applied
//     (default when the record is missing: all -yes; mesh_cut_1.dat
//     writes one -yes in 1D, mesh_cut_2.dat two in 2D).
// The excavated part keeps the remaining mesh in the exact equilibrium
// of the state at the cut (the substituted reaction replaces the
// internal action of the removed material). Measured against the
// Professional: mesh_cut_1 (1D bar2, real cut): remaining mesh nodes
// 1..6, node_force |1| at the boundary node, sigxx = 1 sustained after
// the cut; mesh_cut_2 (2D quad4, geometry line that no node lies
// within): NO element deleted (no-op, same as the Professional); a 2D
// area cut (geometry_quadrilateral over the top rows, probe): the
// boundary row forces 1/9 (interior nodes) / 1/18 (edge nodes) equal
// to the Professional's (the record SIGN differs: our node_force
// consumption in dof.cc applies a stored +F as a -x load while the
// Professional applies it as +x; mesh_cut stores the negated value so
// the equilibrium matches - PENDING: align the node_force record sign
// with the Professional in a dedicated work unit and drop the
// negation).
// Supported elements for the force integration: the tensor-product
// isoparametric family (-bar2/-bar3/-quad4/-quad9/-hex8/-hex27). The
// stress source is the element integration-point record (ELEMENT_DOF),
// which requires materi_stress in the initia and the default
// options_element_dof -yes. Other element types in the cut geometry
// abort with a clear message instead of deleting without the
// equilibrium substitution.
// ---------------------------------------------------------------------

// integration rule of one element (the SAME rule the element itself
// integrated with: pol()/materi(), replicated here like the
// msf_element_rule helper of calcul_force.cc): nper[0..ndim-1]
// integration points per direction + the isoparametric coordinates and
// weights. npol = polynomial order of the element (2 = linear, 3 =
// quadratic, 4 = cubic). pol() (polynom.cc:389-460) integrates the
// -bar2 with the MINIMAL 1-point rule and the other tensor-product
// elements with the MAXIMAL rule (SRI switches to Gauss). Axisymmetric
// groups are not supported by the cut integration.
static void mc_element_rule( long int element, long int element_group,
  long int name, long int npol, long int nnol, long int nper[],
  double iso[][MPOINT], double weight[][MPOINT] )

{
  long int integration_method=-LOBATTO, integration_points=-MAXIMAL,
    axisymmetric=-NO, ldum=0, idim=0;
  double ddum[1];

  db( GROUP_AXISYMMETRIC, element_group, &axisymmetric, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( axisymmetric==-YES ) {
    pri( "Error: control_mesh_cut_geometry is not supported for "
         "axisymmetric groups (the cut force integration is 2D/3D "
         "plane only)." );
    exit(TN_EXIT_STATUS);
  }
  if ( name==-BAR2 )
    integration_points = -MINIMAL;
  else
    integration_points = -MAXIMAL;
  db( GROUP_INTEGRATION_POINTS, element_group, &integration_points, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( GROUP_INTEGRATION_METHOD, element_group, &integration_method, ddum,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( integration_points==-NORMAL ) integration_points = -MAXIMAL;
  if ( sri_active( element, element_group, name, nnol ) )
    integration_method = -GAUSS;
  for ( idim=0; idim<ndim; idim++ ) {
    nper[idim] = ( integration_points==-MINIMAL ? npol-1 : npol );
    if      ( integration_method==-GAUSS )
      integration_gauss( nper[idim], iso[idim], weight[idim] );
    else if ( integration_method==-LOBATTO )
      integration_lobatto( nper[idim], iso[idim], weight[idim] );
    else if ( nper[idim]<npol )
      integration_gauss( nper[idim], iso[idim], weight[idim] );
    else
      integration_lobatto( nper[idim], iso[idim], weight[idim] );
  }
}

// internal nodal forces f_elem[inod*ndim+idim] = int B^T*sigma dV of
// one element at its nodes, from the ELEMENT_DOF integration-point
// stresses of the converged state (the same quantity materi() builds
// into the element right-hand side with the opposite sign). The
// kinematics replicate pol() (polynom.cc): shape functions and their
// isoparametric derivatives via interpolation_polynomial, the Jacobian
// and the volume weight w*2*|detJ| (1D) / w*4*|detJ| (2D) /
// w*8*|detJ| (3D) over the REFERENCE coordinates (NODE at the cut
// time; the total-linear models of the mesh_cut family keep them
// fixed). The Jacobian INVERSE uses the signed determinant (pol()
// inverts the signed Jacobian and takes |det| only for the volume -
// the macro quad4/hex8 record order can be clockwise, i.e. negative
// Jacobian). Stress components (Voigt slots of ELEMENT_DOF, the
// compacted stress_indx() order sxx sxy sxz syy syz szz, primary of
// each component at stres_indx + c*nder):
//   1D: f = int sigxx*dN/dx dV
//   2D: fx = int (sxx*dN/dx + sxy*dN/dy) dV
//       fy = int (sxy*dN/dx + syy*dN/dy) dV
//   3D: the 6-component contraction.
// The node ordering of the tensor-product elements is the tochnog
// record order (quad9/hex27 layer-major row-major, so node index i ->
// (i%npol, (i/npol)%npol, i/(npol*npol)) in the isoparametric
// directions). The ELEMENT_DOF integration-point count follows pol():
// -bar2 has 1 point (MINIMAL rule), the other elements the MAXIMAL
// rule of their group.
static void mc_element_internal_forces( long int element, long int name,
  long int nnol, long int nodes[], double f_elem[] )

{
  long int ldum=0, element_group=0, nper[3], npol=0, ixi=0, ieta=0,
    izeta=0, ip=0, inol=0, idim=0, jdim=0, c=0;
  double iso[3][MPOINT], weight[3][MPOINT], hx[MPOINT], px[MPOINT],
    hy[MPOINT], py[MPOINT], hz[MPOINT], pz[MPOINT],
    coords[MDIM*MNOL], coord[MDIM], jac[9], invjac[9], detj=0.,
    w=0., vol=0., sig[6], dn[3*MNOL], p3[3], ddum[1];
  long int idum[1];

  if ( !( options_element_dof==-YES && stres_indx>=0 ) ) {
    pri( "Error: control_mesh_cut_geometry needs the element "
         "integration-point stresses for the cut force substitution: "
         "materi_stress in the initia section (and the default "
         "options_element_dof -yes) are required." );
    exit(TN_EXIT_STATUS);
  }

  if ( name==-BAR2 || name==-BAR3 ) {
    if ( ndim!=1 ) {
      pri( "Error: bar element in a non-1D mesh_cut." );
      exit(TN_EXIT_STATUS);
    }
    npol = ( name==-BAR2 ? 2 : 3 );
  }
  else if ( name==-QUAD4 || name==-QUAD9 ) {
    if ( ndim!=2 ) {
      pri( "Error: quad element in a non-2D mesh_cut." );
      exit(TN_EXIT_STATUS);
    }
    npol = ( name==-QUAD4 ? 2 : 3 );
  }
  else if ( name==-HEX8 || name==-HEX27 ) {
    if ( ndim!=3 ) {
      pri( "Error: hex element in a non-3D mesh_cut." );
      exit(TN_EXIT_STATUS);
    }
    npol = ( name==-HEX8 ? 2 : 3 );
  }
  else {
    pri( "Error: control_mesh_cut_geometry does not support this "
         "element type in the cut geometry (supported: bar2/bar3, "
         "quad4/quad9, hex8/hex27)." );
    exit(TN_EXIT_STATUS);
  }

  for ( inol=0; inol<nnol; inol++ ) {
    db( NODE, nodes[inol], idum, coord, ldum, VERSION_NORMAL, GET );
    for ( idim=0; idim<ndim; idim++ ) coords[inol*MDIM+idim] = coord[idim];
  }

  db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  mc_element_rule( element, element_group, name, npol, nnol, nper,
    iso, weight );

  if ( !db_active_index( ELEMENT_DOF, element, VERSION_NORMAL ) ) {
    pri( "Error: control_mesh_cut_geometry found an element without "
         "integration-point records (ELEMENT_DOF)." );
    exit(TN_EXIT_STATUS);
  }
  double *edof = db_dbl( ELEMENT_DOF, element, VERSION_NORMAL );

  for ( inol=0; inol<nnol*ndim; inol++ ) f_elem[inol] = 0.;

  if ( ndim==1 ) {
    for ( ixi=0; ixi<nper[0]; ixi++ ) {
      interpolation_polynomial( iso[0][ixi], npol, hx, px );
      ip = ixi;
      jac[0] = 0.;
      for ( inol=0; inol<nnol; inol++ )
        jac[0] += px[inol]*coords[inol*MDIM+0];
      detj = jac[0];
      if ( detj<0. ) detj = -detj;
      if ( detj<1.e-20 ) continue;
      w = weight[0][ixi];
      vol = w*2.*detj;
      sig[0] = edof[ip*nuknwn + stres_indx + 0*nder];
      for ( inol=0; inol<nnol; inol++ ) {
        // dN/dx = (dN/dxi)/J (the signed J: the 1D element may be
        // oriented either way)
        dn[inol*3+0] = px[inol] / jac[0];
        f_elem[inol*1+0] += vol*sig[0]*dn[inol*3+0];
      }
    }
  }
  else if ( ndim==2 ) {
    for ( ieta=0; ieta<nper[1]; ieta++ ) {
      interpolation_polynomial( iso[1][ieta], npol, hy, py );
      for ( ixi=0; ixi<nper[0]; ixi++ ) {
        interpolation_polynomial( iso[0][ixi], npol, hx, px );
        ip = ieta*nper[0] + ixi;
        array_set( jac, 0., 4 );
        for ( inol=0; inol<nnol; inol++ ) {
          jac[0] += px[inol%npol]*hy[inol/npol]*coords[inol*MDIM+0];
          jac[1] += px[inol%npol]*hy[inol/npol]*coords[inol*MDIM+1];
          jac[2] += hx[inol%npol]*py[inol/npol]*coords[inol*MDIM+0];
          jac[3] += hx[inol%npol]*py[inol/npol]*coords[inol*MDIM+1];
        }
        detj = jac[0]*jac[3] - jac[1]*jac[2];
        if ( scalar_dabs(detj)<1.e-20 ) continue;
        // signed-determinant inverse (the record order may be
        // clockwise): the volume below uses |detj|
        invjac[0] =  jac[3]/detj; invjac[1] = -jac[1]/detj;
        invjac[2] = -jac[2]/detj; invjac[3] =  jac[0]/detj;
        if ( detj<0. ) detj = -detj;
        for ( inol=0; inol<nnol; inol++ ) {
          dn[inol*2+0] = invjac[0]*px[inol%npol]*hy[inol/npol]
                       + invjac[1]*hx[inol%npol]*py[inol/npol];
          dn[inol*2+1] = invjac[2]*px[inol%npol]*hy[inol/npol]
                       + invjac[3]*hx[inol%npol]*py[inol/npol];
        }
        w = weight[0][ixi]*weight[1][ieta];
        vol = w*4.*detj;
        sig[0] = edof[ip*nuknwn + stres_indx + 0*nder]; // sxx
        sig[1] = edof[ip*nuknwn + stres_indx + 1*nder]; // sxy
        sig[3] = edof[ip*nuknwn + stres_indx + 3*nder]; // syy
        for ( inol=0; inol<nnol; inol++ ) {
          f_elem[inol*2+0] += vol*( sig[0]*dn[inol*2+0]
                                  + sig[1]*dn[inol*2+1] );
          f_elem[inol*2+1] += vol*( sig[1]*dn[inol*2+0]
                                  + sig[3]*dn[inol*2+1] );
        }
      }
    }
  }
  else {
    assert( ndim==3 );
    for ( izeta=0; izeta<nper[2]; izeta++ ) {
      interpolation_polynomial( iso[2][izeta], npol, hz, pz );
      for ( ieta=0; ieta<nper[1]; ieta++ ) {
        interpolation_polynomial( iso[1][ieta], npol, hy, py );
        for ( ixi=0; ixi<nper[0]; ixi++ ) {
          interpolation_polynomial( iso[0][ixi], npol, hx, px );
          ip = izeta*nper[0]*nper[1] + ieta*nper[0] + ixi;
          array_set( jac, 0., 9 );
          for ( inol=0; inol<nnol; inol++ ) {
            p3[0] = px[inol%npol]*hy[(inol/npol)%npol]*hz[inol/(npol*npol)];
            p3[1] = hx[inol%npol]*py[(inol/npol)%npol]*hz[inol/(npol*npol)];
            p3[2] = hx[inol%npol]*hy[(inol/npol)%npol]*pz[inol/(npol*npol)];
            for ( idim=0; idim<3; idim++ )
              for ( jdim=0; jdim<3; jdim++ )
                jac[idim*3+jdim] += p3[idim]*coords[inol*MDIM+jdim];
          }
          detj = jac[0]*( jac[4]*jac[8] - jac[5]*jac[7] )
               - jac[1]*( jac[3]*jac[8] - jac[5]*jac[6] )
               + jac[2]*( jac[3]*jac[7] - jac[4]*jac[6] );
          if ( !matrix_inverse( jac, invjac, detj, 3 ) ) continue;
          if ( detj<0. ) detj = -detj;
          if ( detj<1.e-20 ) continue;
          for ( inol=0; inol<nnol; inol++ ) {
            p3[0] = px[inol%npol]*hy[(inol/npol)%npol]*hz[inol/(npol*npol)];
            p3[1] = hx[inol%npol]*py[(inol/npol)%npol]*hz[inol/(npol*npol)];
            p3[2] = hx[inol%npol]*hy[(inol/npol)%npol]*pz[inol/(npol*npol)];
            for ( idim=0; idim<3; idim++ ) {
              dn[inol*3+idim] = 0.;
              for ( jdim=0; jdim<3; jdim++ )
                dn[inol*3+idim] += invjac[idim*3+jdim]*p3[jdim];
            }
          }
          w = weight[0][ixi]*weight[1][ieta]*weight[2][izeta];
          vol = w*8.*detj;
          for ( c=0; c<6; c++ )
            sig[c] = edof[ip*nuknwn + stres_indx + c*nder];
          for ( inol=0; inol<nnol; inol++ ) {
            f_elem[inol*3+0] += vol*( sig[0]*dn[inol*3+0]
                                    + sig[1]*dn[inol*3+1]
                                    + sig[2]*dn[inol*3+2] );
            f_elem[inol*3+1] += vol*( sig[1]*dn[inol*3+0]
                                    + sig[3]*dn[inol*3+1]
                                    + sig[4]*dn[inol*3+2] );
            f_elem[inol*3+2] += vol*( sig[2]*dn[inol*3+0]
                                    + sig[4]*dn[inol*3+1]
                                    + sig[5]*dn[inol*3+2] );
          }
        }
      }
    }
  }
}

void mesh_cut( double /*time_current*/ )

{
  long int element=0, max_element=0, max_node=0, inol=0, nnol=0, inod=0,
    length=0, ldum=0, idum[1], icontrol=0, swit=0, in_geometry=0,
    any_element_deleted=0, node_force_switches[MDIM],
    nodes[MNOL], el[MNOL+1], *node_in_geometry=NULL,
    control_mesh_cut_geometry[2];
  double rdum=0., ddum[MDIM], f_elem[3*MNOL];
  long int node_force_length=0, imax=0, force_active=0;
  double *node_force_acc=NULL, node_force_exist[MDIM];

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( !db_active_index( CONTROL_MESH_CUT_GEOMETRY, icontrol,
      VERSION_NORMAL ) )
    return;

  swit = set_swit(-1,-1,"mesh_cut");
  if ( swit ) pri( "In routine MESH_CUT" );

  db( CONTROL_MESH_CUT_GEOMETRY, icontrol, control_mesh_cut_geometry,
    ddum, ldum, VERSION_NORMAL, GET );

  // direction switches of the nodal-force substitution: components
  // with a -yes switch (or all, when the record is not given) are
  // applied to the remaining mesh
  for ( imax=0; imax<MDIM; imax++ ) node_force_switches[imax] = -YES;
  node_force_length = 0;
  db( CONTROL_MESH_CUT_NODE_FORCE, icontrol, idum, ddum,
    node_force_length, VERSION_NORMAL, GET_IF_EXISTS );
  if ( node_force_length>0 ) {
    long int ival[MDIM];
    db( CONTROL_MESH_CUT_NODE_FORCE, icontrol, ival, ddum,
      node_force_length, VERSION_NORMAL, GET );
    for ( imax=0; imax<node_force_length && imax<MDIM; imax++ )
      node_force_switches[imax] = ival[imax];
  }

  db_max_index( NODE, max_node, VERSION_NORMAL, GET );
  length = 1+max_node;
  node_in_geometry = get_new_int( length );
  array_set( node_in_geometry, 0, length );

  // nodes inside the geometry (the delete_geometry membership test:
  // distance to the geometry within its tolerance; the delete/cut
  // projection type counts the filled interior of area geometries)
  for ( inod=0; inod<=max_node; inod++ ) {
    if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
      geometry( inod, ddum, control_mesh_cut_geometry,
        in_geometry, rdum, ddum, rdum, ddum, NODE_START_REFINED,
        CONTROL_MESH_CUT_GEOMETRY, VERSION_NORMAL );
      if ( in_geometry ) node_in_geometry[inod] = 1;
    }
  }

  // the nodal forces to substitute: the force the REMOVED elements
  // exert on the surviving boundary nodes = -f_elem accumulated over
  // the elements to delete (the interior contributions cancel)
  db_max_index( ELEMENT, max_element, VERSION_NORMAL, GET );
  node_force_acc = get_new_dbl( (max_node+1)*MDIM );
  for ( inod=0; inod<=max_node; inod++ )
    for ( imax=0; imax<MDIM; imax++ )
      node_force_acc[inod*MDIM+imax] = 0.;

  // collect + delete the elements fully inside the geometry
  db_version_copy( VERSION_NORMAL, VERSION_TMP );
  for ( element=0; element<=max_element; element++ ) {
    if ( db_active_index( ELEMENT, element, VERSION_NORMAL ) ) {
      db( ELEMENT, element, el, ddum, length, VERSION_NORMAL, GET );
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      long int all_in_geometry = 1;
      for ( inol=0; inol<nnol; inol++ )
        if ( !node_in_geometry[nodes[inol]] ) all_in_geometry = 0;
      if ( all_in_geometry ) {
        // internal forces of this element at its nodes (from the
        // converged IP stresses of the state before the cut)
        mc_element_internal_forces( element, el[0], nnol, nodes,
          f_elem );
        for ( inol=0; inol<nnol; inol++ ) {
          inod = nodes[inol];
          for ( imax=0; imax<ndim; imax++ )
            node_force_acc[inod*MDIM+imax] -= f_elem[inol*ndim+imax];
        }
        delete_element( element, VERSION_NORMAL );
        any_element_deleted = 1;
      }
    }
  }

  if ( any_element_deleted ) {
    // prune the nodes that lost all their elements (nod_nod inside
    // mesh_has_changed; the VERSION_TMP copy protects the node data)
    mesh_has_changed( VERSION_NORMAL );

    // substitute the removed part by its nodal forces on the nodes of
    // the remaining mesh, per direction switch
    for ( inod=0; inod<=max_node; inod++ ) {
      force_active = 0;
      for ( imax=0; imax<ndim; imax++ )
        if ( node_force_switches[imax]==-YES &&
             scalar_dabs(node_force_acc[inod*MDIM+imax])>0. )
          force_active = 1;
      if ( !force_active ) continue;
      if ( !db_active_index( NODE, inod, VERSION_NORMAL ) ) continue;
      array_set( node_force_exist, 0., MDIM );
      db( NODE_FORCE, inod, idum, node_force_exist, ldum,
        VERSION_NORMAL, GET_IF_EXISTS );
      // NOTE (measured, mesh_cut_1.dat + a 2D area-cut probe A/B on
      // the Professional 25-10-2023): the Professional applies a
      // stored node_force +F as a +x external load, our
      // parallel_new_dof_before (dof.cc) applies it as -x
      // (node_rhside -= node_force, opposite sign). The NEGATION below
      // makes the substituted load reproduce the Professional
      // equilibrium (mesh_cut_1: node_force 6 = -1 here, +1 there,
      // both keep the remaining bar at sigxx = +1). PENDING: align the
      // node_force record sign with the Professional in a dedicated
      // work unit and drop this negation.
      for ( imax=0; imax<ndim; imax++ )
        if ( node_force_switches[imax]==-YES )
          node_force_exist[imax] -= node_force_acc[inod*MDIM+imax];
      length = ndim;
      db( NODE_FORCE, inod, idum, node_force_exist, length,
        VERSION_NORMAL, PUT );
    }
  }
  db_version_delete( VERSION_TMP );

  delete[] node_in_geometry;
  delete[] node_force_acc;

  if ( swit ) pri( "Out routine MESH_CUT" );
}

void mesh_delete_small( long int version )

{

  long int icontrol=0, element=0, max_element=0, swit=0, 
    any_element_deleted=0, ldum=0, idum[1];
  double small=0., element_volume=0., ddum[1];

  swit = set_swit(-1,-1,"mesh_delete_small");
  if ( swit ) pri( "In routine MESH_DELETE_SMALL" );

  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
  if ( db_active_index( CONTROL_MESH_DELETE_SMALL, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_MESH_DELETE_SMALL, icontrol, idum, &small, ldum, VERSION_NORMAL, GET );
    db_max_index( ELEMENT, max_element, version, GET );
    for ( element=0; element<=max_element; element++ ) {
      if ( db_active_index( ELEMENT_VOLUME, element, version ) ) {
        db( ELEMENT_VOLUME, element, idum, &element_volume, ldum, version, GET );
        if ( element_volume<small ) {
          if ( swit ) pri( "Deleting element ", element );
          any_element_deleted = 1;
          delete_element( element, version );
        }
      }
    }
  }

  if ( any_element_deleted ) mesh_has_changed( version );

  if ( swit ) pri( "Out routine MESH_DELETE_SMALL" );
}

void mesh_delete_keep( long int icontrol )

{
  // handle control_mesh_delete_element, control_mesh_keep_element,
  // control_mesh_keep_element_group and control_mesh_change_element_group.
  long int ielem=0, max_elem=0, inod=0, max_node=0, length=0, ldum=0, i=0,
    inlist=0, eg_from=0, eg_to=0, eg=0, list[DATA_ITEM_SIZE], nlist=0;
  double ddum[1];

  // control_mesh_change_element_group: change group from eg_from to eg_to
  if ( db_active_index( CONTROL_MESH_CHANGE_ELEMENT_GROUP, icontrol,
      VERSION_NORMAL ) ) {
    long int ceg[2];
    db( CONTROL_MESH_CHANGE_ELEMENT_GROUP, icontrol, ceg, ddum, length,
      VERSION_NORMAL, GET );
    eg_from = ceg[0]; eg_to = ceg[1];
    db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
    for ( ielem=0; ielem<=max_elem; ielem++ ) {
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
        eg = 0;
        db( ELEMENT_GROUP, ielem, &eg, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        if ( eg==eg_from )
          db( ELEMENT_GROUP, ielem, &eg_to, ddum, ldum, VERSION_NORMAL, PUT );
      }
    }
  }

  // control_mesh_delete_element: delete the given element numbers
  if ( db_active_index( CONTROL_MESH_DELETE_ELEMENT, icontrol,
      VERSION_NORMAL ) ) {
    nlist = 0;
    db( CONTROL_MESH_DELETE_ELEMENT, icontrol, list, ddum, nlist,
      VERSION_NORMAL, GET );
    for ( i=0; i<nlist; i++ ) {
      if ( db_active_index( ELEMENT, list[i], VERSION_NORMAL ) )
        delete_element( list[i], VERSION_NORMAL );
    }
  }

  // control_mesh_keep_element: delete all elements except the listed ones
  if ( db_active_index( CONTROL_MESH_KEEP_ELEMENT, icontrol,
      VERSION_NORMAL ) ) {
    nlist = 0;
    db( CONTROL_MESH_KEEP_ELEMENT, icontrol, list, ddum, nlist,
      VERSION_NORMAL, GET );
    db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
    for ( ielem=0; ielem<=max_elem; ielem++ ) {
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
        inlist = array_member( list, ielem, nlist, ldum );
        if ( !inlist ) delete_element( ielem, VERSION_NORMAL );
      }
    }
  }

  // control_mesh_keep_element_group: delete all elements not in the groups
  if ( db_active_index( CONTROL_MESH_KEEP_ELEMENT_GROUP, icontrol,
      VERSION_NORMAL ) ) {
    nlist = 0;
    db( CONTROL_MESH_KEEP_ELEMENT_GROUP, icontrol, list, ddum, nlist,
      VERSION_NORMAL, GET );
    db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
    for ( ielem=0; ielem<=max_elem; ielem++ ) {
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
        eg = 0;
        db( ELEMENT_GROUP, ielem, &eg, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
        inlist = array_member( list, eg, nlist, ldum );
        if ( !inlist ) delete_element( ielem, VERSION_NORMAL );
      }
    }
  }

  // control_mesh_keep_node: delete all nodes except the listed ones
  if ( db_active_index( CONTROL_MESH_KEEP_NODE, icontrol,
      VERSION_NORMAL ) ) {
    nlist = 0;
    db( CONTROL_MESH_KEEP_NODE, icontrol, list, ddum, nlist,
      VERSION_NORMAL, GET );
    db_max_index( NODE, max_node, VERSION_NORMAL, GET );
    for ( inod=0; inod<=max_node; inod++ ) {
      if ( db_active_index( NODE, inod, VERSION_NORMAL ) ) {
        inlist = array_member( list, inod, nlist, ldum );
        if ( !inlist ) delete_node( inod, VERSION_NORMAL );
      }
    }
  }

  mesh_has_changed( VERSION_NORMAL );
}

void mesh_remove( long int icontrol )

{
  // control_mesh_remove:
  //  -method1: remove elements of element_group_0 completely inside elements
  //            of groups element_group_1, element_group_2, ...
  //  -method3: remove elements where all nodes have an mpc (node_mpc exists)
  long int method=0, length=0, ldum=0, ielem=0, max_elem=0,
    jelem=0, inol=0, jnol=0, nnol=0, jnol2=0, ngrp=0, i=0, inlist=0,
    elem_group=0, list[DATA_ITEM_SIZE];
  double ddum[1];
  long int el[1+MNOL], nodes[MNOL], elj[1+MNOL], nodes_j[MNOL];

  db( CONTROL_MESH_REMOVE, icontrol, list, ddum, length, VERSION_NORMAL, GET );
  if ( length<1 ) return;
  method = list[0];
  ngrp = length - 1;

  db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
  if ( max_elem<0 ) return;

  else if ( method==-METHOD1 ) {
    // remove elements of group list[1] completely inside elements of the
    // other groups
    if ( ngrp<1 ) return;
    for ( ielem=0; ielem<=max_elem; ielem++ ) {
      if ( db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) {
        db( ELEMENT, ielem, el, ddum, length, VERSION_NORMAL, GET );
        nnol = length - 1;
        for ( inol=0; inol<nnol; inol++ ) nodes[inol] = el[1+inol];
        elem_group = 0;
        db( ELEMENT_GROUP, ielem, &elem_group, ddum, ldum,
          VERSION_NORMAL, GET_IF_EXISTS );
        if ( elem_group==list[1] ) {
          // check if all nodes of ielem lie inside an element of groups list[2..]
          long int inside_other = 0;
          for ( jelem=0; jelem<=max_elem; jelem++ ) {
            if ( jelem!=ielem && db_active_index( ELEMENT, jelem,
                VERSION_NORMAL ) ) {
              long int jgrp = 0;
              db( ELEMENT_GROUP, jelem, &jgrp, ddum, ldum,
                VERSION_NORMAL, GET_IF_EXISTS );
              inlist = 0;
              for ( i=1; i<ngrp; i++ ) {
                if ( jgrp==list[1+i] ) inlist = 1;
              }
              if ( inlist ) {
                db( ELEMENT, jelem, elj, ddum, length, VERSION_NORMAL, GET );
                long int nnol_j = length - 1;
                for ( jnol=0; jnol<nnol_j; jnol++ ) nodes_j[jnol] = elj[1+jnol];
                // all nodes of ielem must be inside element j (by node match)
                long int all_in = 1;
                for ( inol=0; inol<nnol; inol++ ) {
                  long int found = 0;
                  for ( jnol2=0; jnol2<nnol_j; jnol2++ ) {
                    if ( nodes[inol]==nodes_j[jnol2] ) found = 1;
                  }
                  if ( !found ) all_in = 0;
                }
                if ( all_in ) inside_other = 1;
              }
            }
          }
          if ( inside_other ) delete_element( ielem, VERSION_NORMAL );
        }
      }
    }
  }
  else {
    cout << "Error: control_mesh_remove method must be -method1 or -method3.\n";
    exit(TN_EXIT_STATUS);
  }
  mesh_has_changed( VERSION_NORMAL );
}
