/*
    print_fr.cc - control_print_frd (CalculiX .frd output format).

    Writes the mesh and nodal results in the CalculiX result format
    (.frd), readable by CGX, FreeCAD and prepomax. Only results for 2D
    and 3D isoparametric elements are written (same rule as the
    Professional version).

    File naming (input file excavation.dat, index 100):
      -yes_sequential        -> excavation.frd  (mesh once, results append)
      -separate_index        -> excavation_100.frd
      -separate_sequential   -> excavation_0.frd, excavation_1.frd, ...

    Result names (FreeCAD/prepomax friendly):
      DISP     materi_displacement / materi_velocity_integrated
      STRESS   materi_stress
      TOSTRAIN materi_strain_total
      NDTEMP   condif_temperature
    Other dofs use the Tochnog names truncated to 8 characters.

    control_print_frd_freecad / print_frd_freecad (6.289, 6.979; default
    -no): -yes writes no file unless materi_displacement (or
    materi_velocity_integrated), materi_stress and materi_strain_total are
    initialised (manual 6.288; = Professional 02-08-2026).
    control_print_frd_prepomax / print_frd_prepomax (6.290, 6.980; default
    -no): -yes extrudes a 2D mesh to a thickness of 1 (prepomax cannot plot
    2D elements, manual 6.288): node n gets a copy n+max_node+2 at z = 1
    with the same results, tria3/tria6 -> penta6, quad4/quad9 -> hex8
    (corner nodes) (= Professional 02-08-2026).
*/

#include "tochnog.h"
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <string>

void print_frd( long int icontrol, long int task )

{
  long int inod=0, element=0, max_node=0, max_element=0, nnol=0, name=0,
    length=0, ipuknwn=0, nuknwn_=0, nder_=0,
    element_group=0, swit=0, ldum=0, first=1, nstep=0;
  long int idum[1], *dof_label=NULL, *dof_type=NULL, *dof_scal_vec_mat=NULL,
    *nodes=NULL, *el=NULL;
  double ddum[1], time_current=0., coord[MDIM], *node_dof=NULL;
  char filename[MCHAR], str[MCHAR];

  swit = set_swit(-1,-1,"print_frd");
  if ( swit ) pri( "In routine PRINT_FRD" );

  long int freecad = print_option( PRINT_FRD_FREECAD,
    CONTROL_PRINT_FRD_FREECAD, icontrol, -NO );
  long int prepomax = print_option( PRINT_FRD_PREPOMAX,
    CONTROL_PRINT_FRD_PREPOMAX, icontrol, -NO );
  if ( freecad==-YES && !( ( materi_displacement ||
       materi_velocity_integrated ) && materi_stress &&
       materi_strain_total ) ) {
    if ( swit ) pri( "Out routine PRINT_FRD (freecad: no file)" );
    return;
  }
  // prepomax extrusion of a 2D mesh: copy of node inod (0-based) is
  // written as inod+1+node_shift
  long int extrude = ( prepomax==-YES && ndim==2 );

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 0, 0, idum, idum );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  db_highest_index( ELEMENT, max_element, VERSION_PRINT );
  if ( max_node<0 || max_element<0 ) return;
  long int node_shift = max_node+2;
  long int nnode_out = extrude ? 2*(max_node+1) : max_node+1;
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  nder_ = nder;
  nuknwn_ = nuknwn;

  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  nodes = get_new_int(MAXIMUM_NODE);
  el = get_new_int(MAXIMUM_NODE+1);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );

  // file name
  strcpy( filename, data_file_base );
  if      ( task==-SEPARATE_INDEX && icontrol>=0 ) {
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  else if ( task==-SEPARATE_SEQUENTIAL ) {
    static long int frd_seq=0;
    long_to_a( frd_seq++, str );
    strcat( filename, str );
  }
  strcat( filename, ".frd" );

  {
    std::ifstream fexists( filename );
    first = !fexists.is_open();
    fexists.close();
  }

  std::ofstream out( filename, std::ios::app );
  out.precision(5);
  out.setf(std::ios::scientific);

  if ( first ) {
    // model header record
    char hdr[200];
    char modelname[40];
    strncpy( modelname, data_file_base, 39 );
    modelname[39] = '\0';
    sprintf( hdr, "     1C%-74s\n", modelname );
    out << hdr;
    out << "     1UPGM               Tochnog (GNU)\n";

    // nodal point coordinate block (long format, FORMAT=1)
    out << "     2C" << std::string(18,' ')
        << std::setw(12) << nnode_out << std::string(37,' ') << "1\n";
    for ( long int layer=0; layer<=extrude; layer++ ) {
      for ( inod=0; inod<=max_node; inod++ ) {
        db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
        char line[200];
        sprintf( line, " -1%10ld%12.5E%12.5E%12.5E\n",
          (long)inod+1+layer*node_shift, coord[0], coord[1],
          layer ? 1.0 : ( (ndim>2) ? coord[2] : 0.0 ) );
        out << line;
      }
    }
    out << " -3\n";

    // element definition block (short format)
    long int nelem=0;
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      name = el[0];
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
          VERSION_PRINT, GET );
      else
        element_group = 0;
      long int frd_type = 0;
      if      ( name==-TRIA3 ) frd_type = 7;
      else if ( name==-TRIA6 ) frd_type = 8;
      else if ( name==-QUAD4 ) frd_type = 9;
      else if ( name==-QUAD9 ) frd_type = 10;
      else if ( name==-TET4  ) frd_type = 3;
      else if ( name==-TET10 ) frd_type = 6;
      else if ( name==-HEX8  ) frd_type = 1;
      else if ( name==-HEX27 ) frd_type = 4;
      if ( frd_type==0 ) continue;  // structural elements: not printed
      nelem++;
    }
    out << "     3C" << std::string(18,' ')
        << std::setw(12) << nelem << std::string(37,' ') << "1\n";
    for ( element=0; element<=max_element; element++ ) {
      if ( !db_active_index( ELEMENT, element, VERSION_PRINT ) ) continue;
      db( ELEMENT, element, el, ddum, length, VERSION_PRINT, GET );
      name = el[0];
      nnol = length - 1; array_move( &el[1], nodes, nnol );
      if ( db_active_index( ELEMENT_GROUP, element, VERSION_PRINT ) )
        db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
          VERSION_PRINT, GET );
      else
        element_group = 0;
      long int frd_type = 0;
      if      ( name==-TRIA3 ) frd_type = 7;
      else if ( name==-TRIA6 ) frd_type = 8;
      else if ( name==-QUAD4 ) frd_type = 9;
      else if ( name==-QUAD9 ) frd_type = 10;
      else if ( name==-TET4  ) frd_type = 3;
      else if ( name==-TET10 ) frd_type = 6;
      else if ( name==-HEX8  ) frd_type = 1;
      else if ( name==-HEX27 ) frd_type = 4;
      if ( frd_type==0 ) continue;
      if ( extrude ) {
        // penta6 (frd 2) from the triangle corners, hex8 (frd 1) from the
        // quadrilateral corners (counter-clockwise), then the z = 1 copies
        long int corner[4], ncorner=0;
        if      ( name==-TRIA3 ) {
          ncorner = 3; corner[0]=nodes[0]; corner[1]=nodes[1]; corner[2]=nodes[2];
        }
        else if ( name==-TRIA6 ) {
          ncorner = 3; corner[0]=nodes[0]; corner[1]=nodes[2]; corner[2]=nodes[5];
        }
        else if ( name==-QUAD4 ) {
          ncorner = 4; corner[0]=nodes[0]; corner[1]=nodes[1];
          corner[2]=nodes[3]; corner[3]=nodes[2];
        }
        else {
          assert( name==-QUAD9 );
          ncorner = 4; corner[0]=nodes[0]; corner[1]=nodes[2];
          corner[2]=nodes[8]; corner[3]=nodes[6];
        }
        char line[400];
        sprintf( line, " -1%10ld%5ld%5s%5ld\n -2",
          (long)element+1, ncorner==3 ? 2L : 1L, "0", element_group );
        out << line;
        for ( long int layer=0; layer<=1; layer++ ) {
          for ( long int k=0; k<ncorner; k++ ) {
            char tmp[32];
            sprintf( tmp, "%10ld", corner[k]+1+layer*node_shift );
            out << tmp;
          }
        }
        out << "\n -3\n";
        continue;
      }
      char line[400];
      sprintf( line, " -1%10ld%5ld%5s%5ld\n -2",
        (long)element+1, frd_type, "0", element_group );
      out << line;
      // nodes in CalculiX order (same connectivity as print_vtk/gmsh for
      // linear elements; higher order written with full node list)
      long int order[27], nout=nnol;
      for ( long int k=0; k<nnol; k++ ) order[k] = nodes[k];
      if      ( name==-QUAD4 ) {
        order[0]=nodes[0]; order[1]=nodes[1];
        order[2]=nodes[3]; order[3]=nodes[2];
      }
      else if ( name==-HEX8 ) {
        long int o[8] = {0,1,3,2,4,5,7,6};
        for ( long int k=0; k<8; k++ ) order[k] = nodes[o[k]];
      }
      else if ( name==-QUAD9 ) {
        long int o[8] = {0,1,3,2,4,7,5,6};  // qu8 corners+mid-edges
        nout = 8;
        for ( long int k=0; k<8; k++ ) order[k] = nodes[o[k]];
      }
      for ( long int k=0; k<nout; k++ ) {
        if ( k%10==0 && k>0 ) out << "\n -2";
        char tmp[32];
        sprintf( tmp, "%10ld", order[k]+1 );
        out << tmp;
      }
      out << "\n -3\n";
    }
    out << " -3\n";
  }

  // --- results per time step (append) ---
  // 1PSTEP header, layout of CalculiX frdheader.c:
  // "    1PSTEP" (11) + 13 spaces (cols 12-24) + %12ld counter (25-36) +
  // %12ld increment (37-48) + %12ld step (49-60) + space (61)
  char pstep[96];
  nstep++;
  sprintf( pstep, "    1PSTEP%13s%12ld%12ld%12ld \n",
    "", (long)nstep, icontrol, icontrol );
  out << pstep;

  // 100CL line, layout of frdheader.c (75-char line):
  // "  100CL" (cols 1-7) + %5ld 100+code (8-12) + time %12.*f (13-24) +
  // %12ld numnodes (25-36) + description 12 (37-48) + ictype (57-58) +
  // %5ld kode (59-63) + format '1' (74)
  {
    char cl[80];
    for ( int k=0; k<80; k++ ) cl[k] = ' ';
    char field[24];
    memcpy( cl, "  100CL", 7 );
    sprintf( field, "%5ld", 100L+nstep );            // cols 8-12
    strncpy( &cl[7], field, 5 );
    char tval[13];
    if ( time_current>0. && log10(time_current)>=0. && log10(time_current)<10. ) {
      int ncomma = 10-(int)(floor(log10(time_current))+1.);
      if      ( ncomma<=0 ) sprintf( tval, "%12.0f", time_current );
      else if ( ncomma>=9 ) sprintf( tval, "%12.9f", time_current );
      else                  sprintf( tval, "%12.*f", ncomma, time_current );
    }
    else sprintf( tval, "%12.5E", time_current );
    strncpy( &cl[12], tval, 12 );                    // cols 13-24
    sprintf( field, "%12ld", (long)nnode_out );      // cols 25-36
    strncpy( &cl[24], field, 12 );
    memcpy( &cl[36], "STATIC", 6 );                 // description 37-42
    sprintf( field, "%2ld", 0L );                    // ictype 57-58
    strncpy( &cl[56], field, 2 );
    sprintf( field, "%5ld", icontrol );              // kode 59-63
    strncpy( &cl[58], field, 5 );
    cl[73] = '1';                                    // format col 74
    cl[75] = '\0';
    out << cl << "\n";
  }

  // result blocks: one per dof (scalar/vector/matrix). The loop runs
  // over PRIMARY unknowns as in print_vt.cc and print_dof (BACKLOG
  // C23/C24): the primary unknown ipuknwn starts at iuknwn =
  // ipuknwn*nder_, the dof_* arrays are indexed by iuknwn (they are set
  // for every value AND derivative of the dof) and a vector/matrix dof
  // spans nval primary unknowns. Looping ipuknwn over nuknwn_ wrote every
  // component (and derivative) as its own shifted block and, with
  // derivatives, read node_dof past the end of the record.
  long int npuknwn_ = nder_>0 ? nuknwn_/nder_ : 0, nval=1;
  for ( ipuknwn=0; ipuknwn<npuknwn_; ipuknwn+=nval ) {
    long int iuknwn = ipuknwn*nder_;
    long int dsvm = dof_scal_vec_mat[iuknwn];
    nval = 1;
    if ( dsvm!=-SCALAR && dsvm!=-VECTOR && dsvm!=-MATRIX ) continue;
    if      ( dsvm==-VECTOR ) nval = ndim;
    else if ( dsvm==-MATRIX ) nval = 6;
    long int base_indx = iuknwn;
    long int ncomp=0, ictype=0;
    char dname[9];
    // FRD result name mapping
    long int dt = dof_type[iuknwn];
    if      ( dt==-MATERI_DISPLACEMENT ||
              dt==-MATERI_VELOCITY_INTEGRATED ) strcpy( dname, "DISP" );
    else if ( dt==-MATERI_VELOCITY )             strcpy( dname, "VELO" );
    else if ( dt==-MATERI_STRESS )               strcpy( dname, "STRESS" );
    else if ( dt==-MATERI_STRAIN_TOTAL )         strcpy( dname, "TOSTRAIN" );
    else if ( dt==-CONDIF_TEMPERATURE )          strcpy( dname, "NDTEMP" );
    else {
      strncpy( dname, db_name(dof_label[iuknwn]), 8 );
      dname[8] = '\0';
    }
    char compnames[6][12];
    long int compidx[6] = {0,0,0,0,0,0};
    if      ( dsvm==-SCALAR ) {
      ncomp = 1; ictype = 1;
      strncpy( compnames[0], "V", 8 ); compnames[0][8]='\0';
      sprintf( compnames[0], "%s", db_name(dof_label[iuknwn]) );
    }
    else if ( dsvm==-VECTOR ) {
      ncomp = 3; ictype = 2;
      if ( strcmp(dname,"DISP")==0 ) {
        strcpy( compnames[0], "D1" ); strcpy( compnames[1], "D2" );
        strcpy( compnames[2], "D3" );
      }
      else if ( strcmp(dname,"VELO")==0 ) {
        strcpy( compnames[0], "V1" ); strcpy( compnames[1], "V2" );
        strcpy( compnames[2], "V3" );
      }
      else {
        sprintf( compnames[0], "%sX", dname );
        sprintf( compnames[1], "%sY", dname );
        sprintf( compnames[2], "%sZ", dname );
      }
      compidx[0]=0; compidx[1]=1; compidx[2]=2;
    }
    else {
      ncomp = 6; ictype = 4;
      const char* tens[6] = { "SXX","SYY","SZZ","SXY","SYZ","SZX" };
      if ( strcmp(dname,"STRESS")==0 ) {
        for ( int k=0;k<6;k++ ) strcpy( compnames[k], tens[k] );
      }
      else if ( strcmp(dname,"TOSTRAIN")==0 ) {
        const char* tenst[6] = { "EXX","EYY","EZZ","EXY","EYZ","EZX" };
        for ( int k=0;k<6;k++ ) strcpy( compnames[k], tenst[k] );
      }
      else {
        for ( int k=0;k<6;k++ ) sprintf( compnames[k], "C%1d", k+1 );
      }
      // Voigt order xx,yy,zz,xy,yz,zx via stress_indx
      compidx[0]=stress_indx(0,0); compidx[1]=stress_indx(1,1);
      compidx[2]=stress_indx(2,2); compidx[3]=stress_indx(0,1);
      compidx[4]=stress_indx(1,2); compidx[5]=stress_indx(2,0);
    }

    char r1[80];
    sprintf( r1, " -4  %-8s%4ld    1\n", dname, ncomp );
    out << r1;
    for ( long int k=0; k<ncomp; k++ ) {
      char r2[160];
      long int icind1=0, icind2=0;
      if      ( ictype==2 )      { icind1 = k+1; icind2 = 0; }
      else if ( ictype==4 ) {
        // CalculiX order: SXX(1,1) SYY(2,2) SZZ(3,3) SXY(1,2) SYZ(2,3) SZX(3,1)
        const long int m1[6] = {1,2,3,1,2,3};
        const long int m2[6] = {1,2,3,2,3,1};
        icind1 = m1[k]; icind2 = m2[k];
      }
      sprintf( r2, " -5  %-8s    1%5ld%5ld%5ld\n",
        compnames[k], ictype, icind1, icind2 );
      out << r2;
    }

    // nodal data (prepomax extrusion: the z = 1 copy repeats the values)
    for ( long int inod_out=0; inod_out<nnode_out; inod_out++ ) {
      inod = inod_out%(max_node+1);
      long int layer = inod_out/(max_node+1);
      node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
      std::string line;
      char tmp[32];
      sprintf( tmp, " -1%10ld", inod+1+layer*node_shift );
      line += tmp;
      for ( long int k=0; k<ncomp; k++ ) {
        long int indx = base_indx;
        double val = 0.;
        if      ( dsvm==-VECTOR ) {
          // FRD vectors always have 3 components; in 2D the third one is
          // 0 (it is not part of the dof: k*nder_ would be the next dof)
          indx = base_indx + k*nder_;
          val = ( k<ndim ) ? node_dof[indx] : 0.;
        }
        else if ( dsvm==-MATRIX ) {
          indx = base_indx + compidx[k]*nder_;
          val = node_dof[indx];
        }
        else val = node_dof[indx];
        sprintf( tmp, "%12.5E", val );
        line += tmp;
      }
      out << line << "\n";
    }
    out << " -3\n";
  }

  out << " -3\n";

  out.close();

  db_version_delete( VERSION_PRINT );
  delete[] dof_label;
  delete[] dof_type;
  delete[] dof_scal_vec_mat;
  delete[] nodes;
  delete[] el;

  if ( swit ) pri( "Out routine PRINT_FRD" );
}
