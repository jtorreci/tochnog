/*
    print_tb.cc - control_print_tabular
    Tabular export of nodal results for programmatic post-processing.

    P5-T1: CSV of the current (last) state + optional SQLite primary table.
    Exports the dofs present in the model (displacement, stress, strain,
    pressure, ...), detected via dof_label/dof_scal_vec_mat (same approach
    as print_vtk.cc).

    CSV always available. SQLite only if compiled with SQLITE_USE=1; else a
    warning is printed and the CSV is still produced.
*/

#include "tochnog.h"
#include "sqlite.h"
#include <fstream>
#include <string>

void print_tabular( long int icontrol )
{
  long int i=0, inod=0, idim=0, swit=0, max_node=0, nder_=0;
  long int ldum=0, idum[1], length=0, nuknwn_=0;
  long int *dof_label=NULL, *dof_type=NULL, *dof_scal_vec_mat=NULL;
  double ddum[1], coord[MDIM], *node_dof=NULL, time_current=0.;
  char filename[MCHAR];
  std::string csv_file, sqlite_file;

  swit = set_swit(-1,-1,"print_tabular");
  if ( swit ) pri( "In routine PRINT_TABULAR" );

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  renumbering( VERSION_PRINT, NO, 0, 0, idum, idum );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  nder_ = nder;
  nuknwn_ = nuknwn;

  // Detect present dofs via dof_label / dof_scal_vec_mat
  dof_label = get_new_int(MUKNWN);
  dof_type = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_TYPE, 0, dof_type, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // Count active scalar/vector/matrix dofs
  for ( i=0; i<nuknwn_; i++ ) {
  }

  // CSV output: base<icontrol>.csv
  strcpy( filename, data_file_base );
  if ( icontrol!=-1 ) {
    char str[MCHAR];
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  strcat( filename, ".csv" );
  csv_file = filename;
  for (i=0;i<nuknwn_;i++) fprintf(stderr, "TB dof[%ld] label=%ld svm=%ld\n", i, dof_label[i], dof_scal_vec_mat[i]);
  {
    std::ofstream out(csv_file);
    out.precision(TN_PRECISION);
    out << "node,t";
    for ( i=0; i<nuknwn_; i++ ) {
      if ( dof_scal_vec_mat[i]!=-NO ) {
        if      ( dof_scal_vec_mat[i]==-SCALAR )
          out << "," << db_name(dof_label[i]);
        else if ( dof_scal_vec_mat[i]==-VECTOR ) {
          for ( idim=0; idim<ndim; idim++ )
            out << "," << db_name(dof_label[i]) << "_" << idim;
        }
        else if ( dof_scal_vec_mat[i]==-MATRIX ) {
          long int kdim, ldim;
          for ( kdim=0; kdim<MDIM; kdim++ )
            for ( ldim=0; ldim<MDIM; ldim++ ) {
              if ( stress_indx(kdim,ldim)>=0 && dof_label[stres_indx+stress_indx(kdim,ldim)*nder_]>=0 )
                out << "," << db_name(dof_label[i]) << "_" << kdim << ldim;
            }
        }
      }
    }
    out << "\n";

    for ( inod=1; inod<=max_node; inod++ ) {
      db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
      node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
      out << inod << "," << time_current;
      for ( i=0; i<nuknwn_; i++ ) {
        if ( dof_scal_vec_mat[i]!=-NO ) {
          if      ( dof_scal_vec_mat[i]==-SCALAR )
            out << "," << node_dof[i];
          else if ( dof_scal_vec_mat[i]==-VECTOR ) {
            for ( idim=0; idim<ndim; idim++ )
              out << "," << node_dof[i+idim*nder_];
          }
          else if ( dof_scal_vec_mat[i]==-MATRIX ) {
            long int kdim, ldim;
            for ( kdim=0; kdim<MDIM; kdim++ )
              for ( ldim=0; ldim<MDIM; ldim++ ) {
                long int indx = stress_indx(kdim,ldim);
                if ( indx>=0 && dof_label[stres_indx+indx*nder_]>=0 )
                  out << "," << node_dof[stres_indx+indx*nder_];
              }
          }
        }
      }
      out << "\n";
    }
  }

  // SQLite primary table (optional)
  sqlite_file = csv_file.substr(0, csv_file.size()-4) + ".sqlite";
  {
    SqliteDB* db = sqlite_db_open( sqlite_file.c_str() );
    if ( !db ) {
      pri( "Warning: SQLite not available in this build; CSV written only." );
    } else {
      char sql[2048];
      for ( inod=0; inod<=max_node; inod++ ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        sprintf(sql, "INSERT OR REPLACE INTO primary_data (node,t) VALUES (%ld,%g);",
          inod, time_current);
        db->exec(sql);
      }
      sqlite_db_close( db );
    }
  }
  if ( swit ) pri( "Out routine PRINT_TABULAR" );
  (void)dof_type; (void)length; (void)max_node;
}
