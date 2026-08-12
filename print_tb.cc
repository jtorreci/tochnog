/*
    print_tb.cc - control_print_tabular
    Tabular export of nodal results for programmatic post-processing.

    P5-T3: time series + derived magnitudes (von Mises, Tresca, principal
    stresses) computed in C++ when a stress tensor dof is present.

    The CSV is written in append mode (header only on the first step), each
    call adds a row per node with the current time t. The SQLite uses long
    format primary(node,dof,t,value) and a derived(node,t,...) table.

    Exports the dofs present in the model, detected via dof_scal_vec_mat
    (same approach as print_vtk.cc: dof type is -SCALAR/-VECTOR/-MATRIX).

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
  long int ldum=0, idum[1], nuknwn_=0, sig_indx=-1;
  long int *dof_label=NULL, *dof_scal_vec_mat=NULL;
  long int *old_node_numbers=NULL;
  double ddum[1], coord[MDIM], *node_dof=NULL, time_current=0.;
  char filename[MCHAR];
  std::string csv_file, sqlite_file;

  swit = set_swit(-1,-1,"print_tabular");
  if ( swit ) pri( "In routine PRINT_TABULAR" );

  db_version_copy( VERSION_NORMAL, VERSION_PRINT );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  old_node_numbers = get_new_int(max_node+1);
  renumbering( VERSION_PRINT, YES, 0, 0, old_node_numbers, idum );
  db_highest_index( NODE, max_node, VERSION_PRINT );
  db( TIME_CURRENT, 0, idum, &time_current, ldum, VERSION_NORMAL, GET );
  nder_ = nder;
  nuknwn_ = nuknwn;

  dof_label = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  db( DOF_LABEL, 0, dof_label, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  // stress tensor dof (if present) for derived magnitudes
  for ( i=0; i<nuknwn_; i++ )
    if ( dof_scal_vec_mat[i]==-MATRIX ) { sig_indx = i; break; }

  strcpy( filename, data_file_base );
  if ( icontrol!=-1 ) {
    char str[MCHAR];
    long_to_a( icontrol, str );
    strcat( filename, str );
  }
  strcat( filename, ".csv" );
  csv_file = filename;

  // --- CSV (append per time step) ---
  {
    std::ifstream fexists(csv_file.c_str());
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
      if ( sig_indx>=0 ) out << ",vmises,tresca,sig1,sig2,sig3";
      out << "\n";
    }

    for ( inod=0; inod<=max_node; inod++ ) {
      long int node_nr = old_node_numbers[inod];
      db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
      node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
      out << node_nr << "," << time_current;
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
      if ( sig_indx>=0 ) {
        double sig[6], dout[5];
        sig[0]=node_dof[sig_indx+stress_indx(0,0)*nder_];
        sig[1]=node_dof[sig_indx+stress_indx(1,1)*nder_];
        sig[2]=node_dof[sig_indx+stress_indx(2,2)*nder_];
        sig[3]=node_dof[sig_indx+stress_indx(0,1)*nder_];
        sig[4]=node_dof[sig_indx+stress_indx(0,2)*nder_];
        sig[5]=node_dof[sig_indx+stress_indx(1,2)*nder_];
        calc_derived( sig, dout );
        for ( i=0; i<5; i++ ) out << "," << dout[i];
      }
      out << "\n";
    }
  }

  // --- SQLite (optional) ---
  sqlite_file = csv_file.substr(0, csv_file.size()-4) + ".sqlite";
  {
    SqliteDB* sdb = sqlite_db_open( sqlite_file.c_str() );
    if ( !sdb ) {
      pri( "Warning: SQLite not available in this build; CSV written only." );
    } else {
      char sql[512];
      // metadata (idempotent): mesh, convention, units
      sprintf(sql,
        "INSERT OR REPLACE INTO meta (key,value) VALUES ('ndim','%ld');", ndim);
      sdb->exec(sql);
      sprintf(sql,
        "INSERT OR REPLACE INTO meta (key,value) VALUES ('convention','compression-negative');");
      sdb->exec(sql);
      sprintf(sql,
        "INSERT OR REPLACE INTO meta (key,value) VALUES ('file_base','%s');", data_file_base);
      sdb->exec(sql);
      for ( inod=0; inod<=max_node; inod++ ) {
        long int node_nr = old_node_numbers[inod];
        db( NODE, inod, idum, coord, ldum, VERSION_PRINT, GET );
        double cz = ndim>2 ? coord[2] : 0.0;
        sprintf(sql,
          "INSERT OR REPLACE INTO coords (node,x,y,z) VALUES (%ld,%g,%g,%g);",
          node_nr, coord[0], coord[1], cz);
        sdb->exec(sql);
        node_dof = db_dbl( NODE_DOF, inod, VERSION_PRINT );
        for ( i=0; i<nuknwn_; i++ ) {
          if ( dof_scal_vec_mat[i]==-SCALAR ) {
            sprintf(sql,
              "INSERT OR REPLACE INTO primary_data (node,dof,t,value) VALUES (%ld,'%s',%g,%g);",
              node_nr, db_name(dof_label[i]), time_current, node_dof[i]);
            sdb->exec(sql);
          }
          else if ( dof_scal_vec_mat[i]==-VECTOR ) {
            for ( idim=0; idim<ndim; idim++ ) {
              sprintf(sql,
                "INSERT OR REPLACE INTO primary_data (node,dof,t,value) VALUES (%ld,'%s_%ld',%g,%g);",
                node_nr, db_name(dof_label[i]), idim, time_current, node_dof[i+idim*nder_]);
              sdb->exec(sql);
            }
          }
          else if ( dof_scal_vec_mat[i]==-MATRIX ) {
            long int kdim, ldim;
            for ( kdim=0; kdim<MDIM; kdim++ )
              for ( ldim=0; ldim<MDIM; ldim++ ) {
                long int indx = stress_indx(kdim,ldim);
                sprintf(sql,
                  "INSERT OR REPLACE INTO primary_data (node,dof,t,value) VALUES (%ld,'%s_%ld%ld',%g,%g);",
                  node_nr, db_name(dof_label[i]), kdim, ldim, time_current,
                  node_dof[i+indx*nder_]);
                sdb->exec(sql);
              }
          }
        }
        if ( sig_indx>=0 ) {
          double sig[6], dout[5];
          sig[0]=node_dof[sig_indx+stress_indx(0,0)*nder_];
          sig[1]=node_dof[sig_indx+stress_indx(1,1)*nder_];
          sig[2]=node_dof[sig_indx+stress_indx(2,2)*nder_];
          sig[3]=node_dof[sig_indx+stress_indx(0,1)*nder_];
          sig[4]=node_dof[sig_indx+stress_indx(0,2)*nder_];
          sig[5]=node_dof[sig_indx+stress_indx(1,2)*nder_];
          calc_derived( sig, dout );
          sprintf(sql,
            "INSERT OR REPLACE INTO derived (node,t,vmises,tresca,sig1,sig2,sig3)"
            " VALUES (%ld,%g,%g,%g,%g,%g,%g);",
            node_nr, time_current, dout[0], dout[1], dout[2], dout[3], dout[4]);
          sdb->exec(sql);
        }
      }
      sqlite_db_close( sdb );
    }
  }

  if ( swit ) pri( "Out routine PRINT_TABULAR" );

  delete[] old_node_numbers;
  delete[] dof_label;
  delete[] dof_scal_vec_mat;
}
