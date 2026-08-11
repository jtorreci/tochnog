/*
    print_tb.cc - control_print_tabular
    Tabular export of nodal results for programmatic post-processing.

    P5-T2: time series. The CSV is written in append mode (header only on
    the first step), each call adds a row per node with the current time t.
    The SQLite primary table uses INSERT OR REPLACE on the (node,t) key, so
    each time step accumulates rows.

    Exports the dofs present in the model (displacement, stress, strain,
    pressure, ...), detected via dof_scal_vec_mat (same approach as
    print_vtk.cc: dof type is -SCALAR / -VECTOR / -MATRIX).

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
  long int ldum=0, idum[1], nuknwn_=0;
  long int *dof_label=NULL, *dof_scal_vec_mat=NULL;
  double ddum[1], coord[MDIM], *node_dof=NULL, time_current=0.;
  char filename[MCHAR];
  std::string csv_file, sqlite_file;
  std::ifstream fexists;

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
  dof_scal_vec_mat = get_new_int(MUKNWN);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // CSV output: base<icontrol>.csv (append per time step)
  strcpy( filename, data_file_base );
  if ( icontrol!=-1 ) {
    char str[MCHAR];
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  strcat( filename, ".csv" );
  csv_file = filename;

  {
    // write header only if the file does not exist yet (first call)
    fexists.open(csv_file.c_str());
    bool first = !fexists.is_open();
    fexists.close();

    std::ofstream out(csv_file.c_str(), std::ios::app);
    out.precision(TN_PRECISION);
    if ( first ) {
      out << "node,t";
      for ( i=0; i<nuknwn_; i++ ) {
        if ( dof_scal_vec_mat[i]!=-NO ) {
          if      ( dof_scal_vec_mat[i]==-SCALAR )
            out << "," << db_name(dof_label[i]);
          else if ( dof_scal_vec_mat[i]==-VECTOR ) {
            for ( idim=0; idim<ndim; idim++ ) {
              char str[MCHAR];
              sprintf(str, "%s_%ld", db_name(dof_label[i]), idim);
              out << "," << str;
            }
          }
          else if ( dof_scal_vec_mat[i]==-MATRIX ) {
            long int kdim, ldim;
            for ( kdim=0; kdim<MDIM; kdim++ )
              for ( ldim=0; ldim<MDIM; ldim++ ) {
                char str[MCHAR];
                sprintf(str, "%s_%ld%ld", db_name(dof_label[i]), kdim, ldim);
                out << "," << str;
              }
          }
        }
      }
      out << "\n";
    }

    for ( inod=1; inod<=max_node; inod++ ) {
      db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
      node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
      out << inod << "," << time_current;
      for ( i=0; i<nuknwn_; i++ ) {
        if ( dof_scal_vec_mat[i]==-SCALAR ) {
          out << "," << node_dof[i];
        }
        else if ( dof_scal_vec_mat[i]==-VECTOR ) {
          for ( idim=0; idim<ndim; idim++ )
            out << "," << node_dof[i+idim*nder_];
        }
        else if ( dof_scal_vec_mat[i]==-MATRIX ) {
          long int kdim, ldim;
          for ( kdim=0; kdim<MDIM; kdim++ )
            for ( ldim=0; ldim<MDIM; ldim++ ) {
              long int indx = stress_indx(kdim,ldim);
              out << "," << node_dof[i+indx*nder_];
            }
        }
      }
      out << "\n";
    }
  }

  // SQLite primary table (optional): long format (node,dof,t,value)
  sqlite_file = csv_file.substr(0, csv_file.size()-4) + ".sqlite";
  {
    SqliteDB* db = sqlite_db_open( sqlite_file.c_str() );
    if ( !db ) {
      pri( "Warning: SQLite not available in this build; CSV written only." );
    } else {
      char sql[512];
      for ( inod=1; inod<=max_node; inod++ ) {
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        for ( i=0; i<nuknwn_; i++ ) {
          if ( dof_scal_vec_mat[i]==-SCALAR ) {
            sprintf(sql,
              "INSERT OR REPLACE INTO primary_data (node,dof,t,value) VALUES (%ld,'%s',%g,%g);",
              inod, db_name(dof_label[i]), time_current, node_dof[i]);
            db->exec(sql);
          }
          else if ( dof_scal_vec_mat[i]==-VECTOR ) {
            for ( idim=0; idim<ndim; idim++ ) {
              sprintf(sql,
                "INSERT OR REPLACE INTO primary_data (node,dof,t,value) VALUES (%ld,'%s_%ld',%g,%g);",
                inod, db_name(dof_label[i]), idim, time_current, node_dof[i+idim*nder_]);
              db->exec(sql);
            }
          }
          else if ( dof_scal_vec_mat[i]==-MATRIX ) {
            long int kdim, ldim;
            for ( kdim=0; kdim<MDIM; kdim++ )
              for ( ldim=0; ldim<MDIM; ldim++ ) {
                long int indx = stress_indx(kdim,ldim);
                sprintf(sql,
                  "INSERT OR REPLACE INTO primary_data (node,dof,t,value) VALUES (%ld,'%s_%ld%ld',%g,%g);",
                  inod, db_name(dof_label[i]), kdim, ldim, time_current,
                  node_dof[i+indx*nder_]);
                db->exec(sql);
              }
          }
        }
      }
      sqlite_db_close( db );
    }
  }

  if ( swit ) pri( "Out routine PRINT_TABULAR" );
}
