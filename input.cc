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

#define MCHAR_WORDS 60
#define MSTRING 800
#define MDEFINE 200
#define MARITHMETIC 500
#define MLENGTH 10000

int counter_a=0, counter_b=0, counter_c=0, counter_d=0;
long int reading_define=0, using_define=0, idefine=0, ndefine=0, istring=0, define_nstring[MDEFINE];
char *define_words[MDEFINE], *define_strings[MDEFINE][MSTRING];
long int reading_arithmetic=0, using_arithmetic=0, iarithmetic=0, narithmetic=0, using_if=0;
double arithmetic_values[MARITHMETIC];
char *arithmetic_words[MARITHMETIC];

void input( )

{
  long int i=0, j=0, n=0, iv=0, length=0, max=0, idat=0, itmp=0, idim=0,
    ready=0, d_is_set=0, index=0, range_length=0, last_data_value=0,
    unknown_indx=0, ninitia=0, istr=0, nstr=0,
    print_define=-NO, print_arithmetic=-NO,
    ldum=0, *integer_range=NULL, *range=NULL, 
    *ival=NULL, *dof_label=NULL, *dof_type=NULL, 
    *dof_principal=NULL, *dof_scal_vec_mat=NULL, 
    *dof_amount=NULL, *initialization_values=NULL;
  double d=0., ddum[1], *dval=NULL;
  char str_total[MCHAR], str[MCHAR], str_tmp[MCHAR];

  integer_range = get_new_int(MRANGE);
  range = get_new_int(MRANGE);
  ival = get_new_int(MLENGTH);
  dof_label = get_new_int(DATA_ITEM_SIZE);
  dof_type = get_new_int(MUKNWN);
  dof_principal = get_new_int(MUKNWN);
  dof_scal_vec_mat = get_new_int(MUKNWN);
  dof_amount = get_new_int(MUKNWN);
  initialization_values = get_new_int(MUKNWN);
  dval = get_new_dbl(MLENGTH);

  for ( idefine=0; idefine<MDEFINE; idefine++ ) {
    define_words[idefine] = new char[MCHAR_WORDS];
    for ( istring=0; istring<MSTRING; istring++ ) {
      define_strings[idefine][istring] = new char[MCHAR_WORDS];
    }
  }
  for ( iarithmetic=0; iarithmetic<MARITHMETIC; iarithmetic++ ) {
    arithmetic_words[iarithmetic] = new char[MCHAR_WORDS];
  }
  
  array_set( initialization_values, -EMPTY, DATA_ITEM_SIZE );
  array_set( dof_principal, -NO, MUKNWN );
  array_set( define_nstring, 0, MDEFINE );
  array_set( dof_amount, 0, MUKNWN );

  ofstream out( "tn.log", ios::app );
  out << "\n\nCalculation with data file " << data_file << " starts.";
  out.close();

  if ( freopen( data_file, "r", stdin ) == NULL ) {
    pri( "\nError: cannot open ", data_file );
    exit(TN_EXIT_STATUS);
  }

    /* read initialization part */
  strcpy(str,"");
  while( strcmp(str,"end_initia") ) {
    npuknwn += n; 
    nuknwn = npuknwn*nder;
    if ( npuknwn>MPUKNWN ) {
      pri( "\nError: MPUKNWN too small. Increase it in tochnog.h and recompile.\n" );
      exit(TN_EXIT_STATUS);
    }
    if ( nuknwn>MUKNWN ) {
      pri( "\nError: MUKNWN too small. Increase it in tochnog.h and recompile.\n" );
      exit(TN_EXIT_STATUS);
    }
    if ( nprinc>MPRINC ) {
      pri( "\nError: MPRINC too small. Increase it in tochnog.h and recompile.\n" );
      exit(TN_EXIT_STATUS);
    }
    unknown_indx += n * nder;
    if ( !(cin >> str) ) {
      pri( "\nError in initialization part." );
      exit(TN_EXIT_STATUS);
    }
    input_skip_comment( str ); if ( echo ) cout << str << " ";
    ninitia++; 
    if ( ninitia>DATA_ITEM_SIZE ) {
      pri( "\nError in initialization part. DATA_ITEM_SIZE too small." );
      exit(TN_EXIT_STATUS);
    }
    strcpy(initialization_names[ninitia-1],str);
    n = 0;
    if      ( !strcmp(str,"echo") ) {
      if ( ninitia!=1 ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      if ( !(cin >> str_tmp) ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      input_skip_comment( str_tmp );
      if ( !strcmp(str_tmp,"-yes") ) {
        echo = 1;
        initialization_values[ninitia-1] = -YES;
        cout << "echo -yes";
      }
      else {
        if ( strcmp(str_tmp,"-no") ) {
          pri( "\nError in initialization part." );
          exit(TN_EXIT_STATUS);
        }
        initialization_values[ninitia-1] = -NO;
      }
    }
    else if ( !strcmp(str,"number_of_space_dimensions") || !strcmp(str,"ndim")) {
      if ( ninitia!=2 ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      if ( !(cin >> ndim) ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      if ( echo ) cout << ndim;
      initialization_values[ninitia-1] = ndim;
      if ( ndim<=0 || ndim>MDIM ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
    }
    else if ( !strcmp(str,"derivatives") ) {
      if ( ninitia!=3 ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      derivatives = 1;
      nder = 1 + ndim + 1;
    }
    else if ( !strcmp(str,"beam_rotation") ) {
      beam_rotation = 1;
      rot_indx = unknown_indx;
      if ( ndim==2 ) {
        n = 1;
        array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      }
      else {
        n = 3;
        array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      }
      array_set( &dof_type[rot_indx], -BEAM_ROTATION, n*nder );
      for ( idim=0; idim<n; idim++ ) {
        array_set( &dof_principal[rot_indx+idim*nder], nprinc, nder );
        nprinc += 1;
      }
    }
    else if ( !strcmp(str,"condif_temperature") ) {
      condif_temperature = 1;
      temp_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[temp_indx], -CONDIF_TEMPERATURE, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      array_set( &dof_principal[temp_indx], nprinc, n*nder );
      nprinc += n;
    }
    else if ( !strcmp(str,"groundflow_pressure") ) {
      groundflow_pressure = 1;
      pres_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[pres_indx], -GROUNDFLOW_PRESSURE, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      array_set( &dof_principal[pres_indx], nprinc, n*nder );
      nprinc += n;
    }
    else if ( !strcmp(str,"groundflow_velocity") ) {
      groundflow_velocity = 1;
      gvel_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[gvel_indx], -GROUNDFLOW_VELOCITY, n*nder );
      array_set( &dof_scal_vec_mat[gvel_indx], -VECTOR, n*nder );
    }
    else if ( !strcmp(str,"materi_damage") ) {
      materi_damage = 1;
      dam_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[dam_indx], -MATERI_DAMAGE, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_density") ) {
      materi_density = 1;
      dens_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[dens_indx], -MATERI_DENSITY, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_diffusion") ) {
      materi_diffusion = 1;
      diff_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[diff_indx], -MATERI_DIFFUSION, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }             
    else if ( !strcmp(str,"materi_displacement") ) {
      materi_displacement = 1;
      dis_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[dis_indx], -MATERI_DISPLACEMENT, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
    }
    else if ( !strcmp(str,"materi_history_variables") ) {
      if ( !(cin >> materi_history_variables) ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      assert( materi_history_variables>0 );
      if ( echo ) cout << materi_history_variables;
      initialization_values[ninitia-1] = materi_history_variables;
      hisv_indx = unknown_indx;
      n = materi_history_variables;
      array_set( &dof_type[hisv_indx], -MATERI_HISTORY_VARIABLES, 
        materi_history_variables*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_maxwell_stress") ) {
      if ( !(cin >> materi_maxwell_stress) ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      assert( materi_maxwell_stress>0 );
      assert( materi_maxwell_stress<MMAXWELL );
      if ( echo ) cout << materi_maxwell_stress;
      initialization_values[ninitia-1] = materi_maxwell_stress;
      mstres_indx = unknown_indx;
      n = 6 * materi_maxwell_stress;
      array_set( &dof_type[mstres_indx], -MATERI_MAXWELL_STRESS, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_f") ) {
      materi_plasti_f = 1;
      f_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[f_indx], -MATERI_PLASTI_F, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_f_nonlocal") ) {
      materi_plasti_f_nonlocal = 1;
      fn_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[fn_indx], -MATERI_PLASTI_F_NONLOCAL, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_intergranular") ) {
      materi_strain_intergranular = 1;
      epi_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[epi_indx], 
        -MATERI_STRAIN_INTERGRANULAR, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_incremental_substeps") ) {
      materi_plasti_incremental_substeps = 1;
      substeps_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[substeps_indx], -MATERI_PLASTI_INCREMENTAL_SUBSTEPS, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_kappa") ) {
      materi_plasti_kappa = 1;
      kap_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[kap_indx], -MATERI_PLASTI_KAPPA, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_rho") ) {
      materi_plasti_rho = 1;
      rho_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[rho_indx], -MATERI_PLASTI_RHO, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_softvar_local") ) {
      materi_plasti_softvar_local = 1;
      svloc_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[svloc_indx], -MATERI_PLASTI_SOFTVAR_LOCAL, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_softvar_nonlocal") ) {
      materi_plasti_softvar_nonlocal = 1;
      svnonloc_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[svnonloc_indx], -MATERI_PLASTI_SOFTVAR_NONLOCAL, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_strainenergy") ) {
      materi_strainenergy = 1;
      ener_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[ener_indx], -MATERI_STRAINENERGY, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_elasti") ) {
      materi_strain_elasti = 1;
      epe_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[epe_indx], -MATERI_STRAIN_ELASTI, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti") ) {
      materi_strain_plasti = 1;
      epp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[epp_indx], -MATERI_STRAIN_PLASTI, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_total") ) {
      materi_strain_total = 1;
      ept_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[ept_indx], -MATERI_STRAIN_TOTAL, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_stress") ) {
      materi_stress = 1;
      stres_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[stres_indx], -MATERI_STRESS, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_velocity") ) {
      materi_velocity = 1;
      vel_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[vel_indx], -MATERI_VELOCITY, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
      for ( idim=0; idim<ndim; idim++ ) {
        array_set( &dof_principal[vel_indx+idim*nder], nprinc, nder );
        nprinc += 1;
      }
    }
    else if ( !strcmp(str,"materi_velocity_integrated") ) {
      materi_velocity_integrated = 1;
      veli_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[veli_indx], -MATERI_VELOCITY_INTEGRATED, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
    }
    else if ( !strcmp(str,"materi_void_fraction") ) {
      materi_void_fraction = 1;
      void_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[void_indx], -MATERI_VOID_FRACTION, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_work") ) {
      materi_work = 1;
      work_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[work_indx], -MATERI_WORK, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"maxwell_e") ) {
      maxwell_e = 1;
      maxe_indx = unknown_indx;
      n = MDIM;
      array_set( &dof_type[maxe_indx], -MAXWELL_E, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
    }                      
    else if ( !strcmp(str,"maxwell_fe") ) {
      maxwell_fe = 1;
      maxfe_indx = unknown_indx;
      n = MDIM;
      array_set( &dof_type[maxfe_indx], -MAXWELL_FE, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
      array_set( &dof_principal[unknown_indx], nprinc, n*nder );
      nprinc += n;
    }                      
    else if ( !strcmp(str,"maxwell_er") ) {
      maxwell_er = 1;
      maxer_indx = unknown_indx;
      n = MDIM;
      array_set( &dof_type[maxer_indx], -MAXWELL_ER, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
      array_set( &dof_principal[unknown_indx], nprinc, n*nder );
      nprinc += n;
    }                      
    else if ( !strcmp(str,"maxwell_ei") ) {
      maxwell_ei = 1;
      maxei_indx = unknown_indx;
      n = MDIM;
      array_set( &dof_type[maxei_indx], -MAXWELL_EI, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
      array_set( &dof_principal[unknown_indx], nprinc, n*nder );
      nprinc += n;
    }                      
    else if ( !strcmp(str,"residue") ) {
      residue = 1;
      res_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[res_indx], -RESIDUE, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"wave_scalar") ) {
      wave_scalar = 1;
      scal_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[scal_indx], -WAVE_SCALAR, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"wave_fscalar") ) {
      wave_fscalar = 1;
      fscal_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[fscal_indx], -WAVE_FSCALAR, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      array_set( &dof_principal[fscal_indx], nprinc, n*nder );
      nprinc += n;
    }
    else if ( !strcmp(str,"number_of_integration_points") ) {
      if ( !(cin >> npointmax) ) {
        pri( "\nError in initialization part." );
        exit(TN_EXIT_STATUS);
      }
      if ( echo ) cout << ndim;
    }
    else if ( strcmp(str,"end_initia") ) {
      pri( "\nError in initialization part." );
      exit(TN_EXIT_STATUS);
    }
    array_set( &dof_amount[unknown_indx], n, n*nder );
    if ( echo ) cout << "\n";
  }
  nuknwn = npuknwn * nder;
  if ( ndim<0 || ndim>MDIM ) {
    pri( "\nError: number_of_space_dimensions is not correct." );
    exit(TN_EXIT_STATUS);
  }

  db_initialize( dof_type, dof_label );
  db( INITIALIZATION_VALUES, 0, initialization_values, ddum, 
    ninitia, VERSION_NORMAL, PUT );

    /* check initialization part */
  if ( residue ) 
    check_unknown( "derivatives", YES, CHECK_USAGE_AND_ERROR );
  if ( materi_damage )  {
    check_unknown( "materi_strain_elasti", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
  }
  if ( materi_displacement ) 
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
  if ( materi_work ) {
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
  }
  if ( materi_plasti_f ) 
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
  if ( materi_plasti_f_nonlocal ) 
    check_unknown( "materi_plasti_f", YES, CHECK_USAGE_AND_ERROR );
  if ( maxwell_er )
    check_unknown( "maxwell_ei", YES, CHECK_USAGE_AND_ERROR );
  if ( maxwell_ei )
    check_unknown( "maxwell_er", YES, CHECK_USAGE_AND_ERROR );
  if ( maxwell_e )
    check_unknown( "maxwell_fe", YES, CHECK_USAGE_AND_ERROR );
  if ( maxwell_fe )
    check_unknown( "maxwell_e", YES, CHECK_USAGE_AND_ERROR );   
  if ( wave_fscalar ) 
    check_unknown( "wave_scalar", YES, CHECK_USAGE_AND_ERROR );
  if ( wave_scalar ) 
    check_unknown( "wave_fscalar", YES, CHECK_USAGE_AND_ERROR );
  if ( materi_diffusion ) 
    check_unknown( "materi_density", NO, CHECK_USAGE_AND_ERROR );
  if ( materi_density ) 
    check_unknown( "materi_diffusion", NO, CHECK_USAGE_AND_ERROR );
  if ( materi_strainenergy ) {
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_strain_elasti", YES, CHECK_USAGE_AND_ERROR );
  }


    /* read data part */
  input_read_string( echo, str, d, d_is_set );
  input_skip_comment( str );
  while ( strcmp(str,"end_data") ) {

      // data item name
    input_skip_comment( str );
    if ( echo ) cout << str << " ";
    idat = db_number( str ); 
    if ( idat<0 ) {
      pri( "\nError in data part." );
      pri( "I do not know ", str );
      exit(TN_EXIT_STATUS);
    }
    check( idat, CHECK_USAGE_AND_ERROR );

      // index
    if ( db_no_index(idat ) )
      range[0] = 0;
    else {
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment(str);
      if ( echo ) cout << " " << str << " ";
      if ( !strcmp(str,"-ra") ) 
         range_scan( echo, d, d_is_set, range, range_length );
      else {
        if ( !string_isinteger(str) ) {
          pri( "\nError in data part." );
          pri( "Illegal index ", str );
          exit(TN_EXIT_STATUS);
        }
        range[0] = atoi(str);
      }
    }


      // data item values
    length = last_data_value = istr = nstr = 0;
    for ( iv=0; !last_data_value; iv++ ) {
      if ( istr>=nstr ) {
        input_read_string( echo, str_total, d, d_is_set );
        istr = 0;
        nstr = strlen(str_total);
      }
      ready = 0;
      while ( !ready ) {
        for ( j=0, i=istr; !ready; i++ ) {
          if ( str_total[i]==',' )
            ready = 1;
          else {
            str[j] = str_total[i];
            str[j+1] = '\0';
          }
          istr = i+1;
          j++;
          if ( istr==nstr ) ready = 1;
        }
      }
      if ( !d_is_set ) input_skip_comment( str );
      if      ( db_fixed_length(idat) ) {
        last_data_value = ( iv==(db_data_length(idat)) );
      }
      else if ( !strcmp(str,"end_data") )
        last_data_value = 1;
      else if ( !strcmp(str,"start_if") )
        last_data_value = 1;
      else if ( !d_is_set ) 
       last_data_value = ( db_number(str)>=0 );
      if ( !last_data_value ) {
        if ( echo ) {
          if ( d_is_set ) cout << " " << d << " ";
          else cout << " " << str << " ";
        }
        if ( db_type(idat)==INTEGER ) {
          if ( str[0]=='-' ) {
            itmp = db_number( &str[1] );
            if ( itmp<0 ) {
              pri( "\nError in data part." );
              pri( "I do not know ", str );
              exit(TN_EXIT_STATUS);
            }
            check( itmp, CHECK_USAGE_AND_ERROR );
            ival[iv] = -itmp;
          }
          else {
            if ( d_is_set )
              ival[iv] = (long int) d;
            else {
              if ( !string_isinteger(str) && !string_isdouble(str) ) {
                pri( "\n\nError in data part." );
                pri( "Problem reading ", db_name(idat) );
                pri( "I don't know what to do with ", str );
                if ( db_fixed_length(idat) )
                  pri( "Number of data values expected ", db_data_length(idat) );
                exit(TN_EXIT_STATUS);
              }
              ival[iv] = atoi(str);
            }
          }
        }
        else if ( db_type(idat)==DOUBLE_PRECISION ) {
          if ( d_is_set )
            dval[iv] = d;
          else {
            if ( !string_isdouble(str) ) {
              pri( "\n\nError in data part." );
              pri( "Problem reading ", db_name(idat) );
              pri( "I don't know what to do with ", str );
              if ( db_fixed_length(idat) )
                pri( "Number of data values expected ", db_data_length(idat) );
              exit(TN_EXIT_STATUS);
            }
            dval[iv] = atof(str);
          }
        }
        length = iv + 1;
        if ( length>MLENGTH-1 ) {
          pri( "\nError in data part. Record length exceeds MLENGTH." );
          pri( "\nIncrease MLENGTH in input.cc." );
          exit(TN_EXIT_STATUS);
        }
      }
    }
    if ( echo ) cout << "\n";

      // store 
    if ( range[0]==-RA ) range_expand( range, integer_range, ldum, range_length );
    db_max_index( idat, max, VERSION_NORMAL, GET );
    for ( i=0, ready=0; !ready; i++) {
      if ( range[0]==-RA ) 
        index = integer_range[i];
      else 
        index = range[0];
      if ( db_active_index(idat,index,VERSION_NORMAL) ) db_error( idat, index );
      db( idat, index, ival, dval, length, VERSION_NORMAL, PUT );
      if ( range[0]==-RA ) 
        ready = ( (i+1)==range_length );
      else 
        ready = 1;
    }
  }
  if ( echo ) cout << "end_data" << "\n\n";

  if ( using_define ) {
    pri( "\n\nError: start_define not closed." );
    exit(TN_EXIT_STATUS);
  }
  if ( using_arithmetic ) {
    pri( "\n\nError: start_arithmetic not closed." );
    exit(TN_EXIT_STATUS);
  }
  if ( using_if ) {
    pri( "\n\nError: start_if not closed." );
    exit(TN_EXIT_STATUS);
  }

  if ( nuknwn>0 ) {
    db( DOF_AMOUNT, 0, dof_amount, ddum, nuknwn, VERSION_NORMAL, PUT );
    db( DOF_LABEL, 0, dof_label, ddum, nuknwn, VERSION_NORMAL, PUT );
    db( DOF_TYPE, 0, dof_type, ddum, nuknwn, VERSION_NORMAL, PUT );
    db( DOF_SCAL_VEC_MAT, 0, dof_scal_vec_mat, ddum, nuknwn, VERSION_NORMAL, PUT );
    db( DOF_PRINCIPAL, 0, dof_principal, ddum, nuknwn, VERSION_NORMAL, PUT );
  }

  db( PRINT_ARITHMETIC, 0, &print_arithmetic, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( print_arithmetic==-YES && narithmetic>0 ) {
    cout << "\nArithmetic:\n";
    for ( iarithmetic=0; iarithmetic<narithmetic; iarithmetic++ ) {
      cout << arithmetic_words[iarithmetic] << " ";
      cout << arithmetic_values[iarithmetic] << "\n";
    }
    cout << "\n";
  }

  db( PRINT_DEFINE, 0, &print_define, ddum, ldum, 
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( print_define==-YES && ndefine>0 ) {
    cout << "\nDefine:\n";
    for ( idefine=0; idefine<ndefine; idefine++ ) {
      cout << define_words[idefine] << " ";
      for ( istring=0; istring<define_nstring[idefine]; istring++ ) {
        cout << define_strings[idefine][istring] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  }

  input_check_required();

    // de-allocate

  for ( idefine=0; idefine<MDEFINE; idefine++ ) {
    delete[] define_words[idefine];
    for ( istring=0; istring<MSTRING; istring++ ) {
      delete[] define_strings[idefine][istring];
    }
  }

  if ( echo ) date();

  delete[] integer_range;
  delete[] range;
  delete[] ival;
  delete[] dof_label;
  delete[] dof_type;
  delete[] dof_principal;
  delete[] dof_scal_vec_mat;
  delete[] dof_amount;
  delete[] initialization_values;
  delete[] dval;

}

void input_read_string( long int echo, char str[], double &d, long int &d_is_set )

{

  long int i=0, length=0, ready=0, found=0, minus=0, plus=0, multiply=0, divide=0,
    apply_if=0;
  double d_new=0.;
  char *s, str_tmp[MCHAR];

  start_of_input_read_string:
  d_is_set = 0;

  if ( !using_define ) {
    if ( !(cin >> str) ) {
      pri( "\nError in data part." );
      pri( "Unexpected end of input detected (internal location a)." );
      pri( "Last word read ", str );
      exit(TN_EXIT_STATUS);
    }
    string_convert_to_lower_case( str );
  }

  if ( !strcmp(str,"start_define") ) {
    if ( echo ) cout << "\nstart_define\n";
    if( ndefine==MDEFINE ) {
      cout << "\nError: maximal number of defines is " << MDEFINE << "\n";
      pri( "Increase MDEFINE in input.cc" );
      exit(TN_EXIT_STATUS);
    }
    if ( reading_define ) {
      pri( "\nError in data part." );
      pri( "Error detected for define: ", define_words[ndefine] );
      pri( "You cannot start a define inside a define block." );
      exit(TN_EXIT_STATUS);
    }
    if ( reading_arithmetic ) {
      pri( "\nError in data part." );
      pri( "Error detected for define: ", define_words[ndefine] );
      pri( "You cannot start a define inside an arithmetic block." );
      exit(TN_EXIT_STATUS);
    }
    reading_define = 1;
      // read word
    if ( !(cin >> str) ) {
      pri( "\nError in data part." );
      pri( "Unexpected end of input (internal location b)." );
      pri( "Last word read ", str );
      exit(TN_EXIT_STATUS);
    }
    string_convert_to_lower_case( str );
    input_skip_comment( str );
    if ( echo ) cout << str << " ";
    if ( strlen(str)>MCHAR_WORDS-1 ) {
      pri( "Word to long ", str );
      length = MCHAR_WORDS-1; pri( "maximum allowed length is ", length );
      exit(TN_EXIT_STATUS);
    }
    for ( idefine=0; idefine<ndefine; idefine++ ) {
      if ( !strcmp(str,define_words[idefine]) ) {
        pri( "\nDefined word already exists ", str );
        exit(TN_EXIT_STATUS);
      }
    }
    strcpy( define_words[ndefine], str );
      // read strings
    ready = 0;
    do {
      if ( !(cin >> str) ) {
        pri( "\nError in data part." );
        pri( "Unexpected end of input (internal location c)." );
        pri( "Last word read ", str );
        exit(TN_EXIT_STATUS);
      }
      string_convert_to_lower_case( str );
      input_skip_comment( str );
      if      ( !strcmp("end_arithmetic",str) ) {
        pri( "\nError in data part." );
        pri( "Error detected for define: ", define_words[ndefine] );
        pri( "You cannot close a start_define with an end_arithmetic." );
        exit(TN_EXIT_STATUS);
      }
      else if ( !strcmp("start_define",str) ) {
        pri( "\nError in data part." );
        pri( "Error detected for define: ", define_words[ndefine] );
        pri( "You cannot start a define inside a define block." );
        exit(TN_EXIT_STATUS);
      }
      else if ( !strcmp("start_arithmetic",str) ) {
        pri( "\nError in data part." );
        pri( "Error detected for define: ", define_words[ndefine] );
        pri( "You cannot start an arithmetic inside a define block." );
        exit(TN_EXIT_STATUS);
      }
      else if ( strcmp("end_define",str) ) {
        if ( echo ) cout << " " << str << " ";
        if ( !strcmp(str,"plus") || !strcmp(str,"minus") ||
           !strcmp(str,"multiply") || !strcmp(str,"divide") ) {
          pri( "\nError in data part. Do you want start_arithmetic ..  end_arithmetic?" );
          exit(TN_EXIT_STATUS);
        }
        istring = define_nstring[ndefine];
        if ( istring==MSTRING-1 ) {
          cout << "\nError: maximal number of strings is " << MSTRING << "\n";
          pri( "Error detected for define: ", define_words[ndefine] );
          pri( "Maybe the define is not correct." );
          pri( "Else, increase MSTRING in input.cc and recompile." );
          exit(TN_EXIT_STATUS);
        }
        strcpy( define_strings[ndefine][istring], str );
        istring++;
        define_nstring[ndefine] = istring;
      }
      else {
        if ( define_nstring[ndefine]==0 ) {
          pri( "\n\nError: no strings defined for define block.\n\n" );
          exit(TN_EXIT_STATUS);
        }
        if ( echo ) cout << "\nend_define\n";
        reading_define = 0;
        ready = 1;
      }
    }
    while ( !ready );
    ndefine++;
    goto start_of_input_read_string;
  }
  else if ( !strcmp(str,"start_arithmetic") ) {
    if ( echo ) cout << str << "\n";
    if ( reading_define ) {
      pri( "\nError in data part." );
      pri( "Error detected for arithmetic: ", arithmetic_words[narithmetic] );
      pri( "You cannot start an arithmetic inside a define block." );
      exit(TN_EXIT_STATUS);
    }
    if ( reading_arithmetic ) {
      pri( "\nError in data part." );
      pri( "Error detected for arithmetic: ", arithmetic_words[narithmetic] );
      pri( "You cannot start an arithmetic inside an arithmetic block." );
      exit(TN_EXIT_STATUS);
    }
    reading_arithmetic = 1;
    if ( !(cin >> str) ) {
      pri( "\nError in data part." );
      pri( "Unexpected end of input (internal location d)." );
      pri( "Last word read ", str );
      exit(TN_EXIT_STATUS);
    }
    string_convert_to_lower_case( str );
    input_skip_comment( str );         
    if ( echo ) cout << " " << str << " ";
    if ( strlen(str)>MCHAR_WORDS-1 ) {
      pri( "Word to long ", str );
      length = MCHAR_WORDS-1; pri( "maximum allowed length is ", length );
      exit(TN_EXIT_STATUS);
    }
    for ( iarithmetic=0; iarithmetic<narithmetic; iarithmetic++ ) {
      if ( !strcmp(str,arithmetic_words[iarithmetic]) ) {
        pri( "\nArithmetic word already exists ", str );
        exit(TN_EXIT_STATUS);
      }
    }
    strcpy( arithmetic_words[narithmetic], str );
    if ( !(cin >> str) ) {
      pri( "\nError in data part." );
      pri( "Unexpected end of input (internal location e)." );
      pri( "Last word read ", str );
      exit(TN_EXIT_STATUS);
    }            
    string_convert_to_lower_case( str );
    input_skip_comment( str );         
    if ( echo ) cout << " " << str << " ";
    found = 0;
    for ( i=0; i<narithmetic; i++ ) {
      if ( !strcmp( arithmetic_words[i], str ) ) {
        d = arithmetic_values[i];
        found = 1;
      }
    }             
    if ( !found ) {
      if ( !string_isdouble(str) ) {
        pri( "\nError in data part." );
        pri( "Unexpected end of input (internal location f)." );
        pri( "Last word read ", str );
        exit(TN_EXIT_STATUS);
      }
      d = atof(str);
    }
    ready = 0;
    do {
      if ( !(cin >> str) ) {
        pri( "\nError in data part." );
        pri( "Unexpected end of input (internal location g)." );
        pri( "Last word read ", str );
        exit(TN_EXIT_STATUS);
      }
      string_convert_to_lower_case( str );
      input_skip_comment( str );         
      if      ( !strcmp("end_define",str) ) {
        pri( "\nError in data part." );
        pri( "Error detected for arithmetic: ", arithmetic_words[narithmetic] );
        pri( "You cannot close a start_arithmetic with an end_define." );
        exit(TN_EXIT_STATUS);
      }
      else if ( !strcmp("start_define",str) ) {
        pri( "\nError in data part." );
        pri( "Error detected for arithmetic: ", arithmetic_words[narithmetic] );
        pri( "You cannot start a define inside a arithmetic block." );
        exit(TN_EXIT_STATUS);
      }
      else if ( !strcmp("start_arithmetic",str) ) {
        pri( "\nError in data part." );
        pri( "Error detected for arithmetic: ", arithmetic_words[narithmetic] );
        pri( "You cannot start a arithmetic inside a define block." );
        exit(TN_EXIT_STATUS);
      }
      else if ( !strcmp("end_arithmetic",str) ) {
        if ( echo ) cout << "\n" << str << "\n";
        arithmetic_values[narithmetic] = d;
        reading_arithmetic = 0;
        ready = 1;
      }     
      else {                               
        if ( echo ) cout << " " << str << " ";
        if      ( !strcmp(str,"plus") ) 
          plus = 1;
        else if ( !strcmp(str,"minus") ) 
          minus = 1;
        else if ( !strcmp(str,"multiply") ) 
          multiply = 1;
        else if ( !strcmp(str,"divide") ) 
          divide = 1;
        else {
          found = 0;
          for ( i=0; i<narithmetic; i++ ) {
            if ( !strcmp( arithmetic_words[i], str ) ) {
              d_new = arithmetic_values[i];
              found = 1;
            }
          }
          if ( !found ) {
            if ( !string_isdouble(str) ) {
              pri( "\nError in data part." );
              pri( "I don't know what to do with ", str );
              exit(TN_EXIT_STATUS);
            }
            d_new = atof(str);
          }                   
          if      ( plus )
            d += d_new;
          else if ( minus )
            d -= d_new;
          else if ( multiply )
            d *= d_new;
          else {
            if ( !divide ) {
              pri( "Error in arithmetic." );
              exit(TN_EXIT_STATUS );
            }
            d /= d_new;
          }
          plus = minus = multiply = divide = 0;
        }
      }
    }
    while ( !ready );
    narithmetic++;                                  
    if ( narithmetic>MARITHMETIC-1 ) {
      pri( "Error: MARITHMETIC in input.cc too small. Increase it." );
      exit(TN_EXIT_STATUS );
    }
    goto start_of_input_read_string;
  }
  else if ( !strcmp(str,"start_if") ) {
    if ( using_if ) {
      pri( "Error, start_if cannot be nested." );
      exit(TN_EXIT_STATUS);
    }
    using_if = 1;
    if ( !(cin >> str) ) {
      pri( "\nError in data part." );
      pri( "Unexpected end of input (internal location g)." );
      pri( "Last word read ", str );
      exit(TN_EXIT_STATUS);
    }
    string_convert_to_lower_case( str );
    s = str;
    if ( *s=='(' ) {
      pri( "Error, comment not allowed inside start_if ... end_if" );
      exit(TN_EXIT_STATUS);
    }
    found = apply_if = 0;
    for ( idefine=0; idefine<ndefine; idefine++ ) {
      if ( !strcmp( define_words[idefine], str ) ) {
        if      ( !strcmp(define_strings[idefine][0],"true" ) ) {
          found = 1;
          apply_if = 1;
        }
        else if ( !strcmp(define_strings[idefine][0],"false" ) ) {
          found = 1;
          apply_if = 0;
        }
        else {
          pri( "\n\nError, illegal value for ", str );
          pri( "Use either true or false." );
          exit(TN_EXIT_STATUS);
        }
      }
    }
    if ( !found ) {
      pri( "\n\nError, start_if cannot find", str );
      exit(TN_EXIT_STATUS);
    }
    if ( apply_if ) {
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
    }
    else {
      loop_if:
      if ( !(cin >> str) ) {
        pri( "\nError in data part." );
        pri( "Unexpected end of input (internal location g)." );
        pri( "Last word read ", str );
        exit(TN_EXIT_STATUS);
      }
      string_convert_to_lower_case( str );
      s = str;
      if ( *s=='(' ) {
        pri( "Error, comment not allowed inside start_if ... end_if" );
        exit(TN_EXIT_STATUS);
      }
      if ( strcmp(str,"end_if") ) {
        goto loop_if;
      }
      else {
        using_if = 0;
        input_read_string( echo, str, d, d_is_set );
        input_skip_comment( str );
      }
    }
  }
  else if ( !strcmp(str,"end_if") && using_if ) {
    using_if = 0;
    input_read_string( echo, str, d, d_is_set );
    input_skip_comment( str );
  }
  else {
    strcpy( str_tmp, &str[0] );          
    if ( !using_define ) {
      minus = 0;
      if      ( str[0]=='-' ) {
        minus = 1;
        strcpy( str_tmp, &str[1] );
      }
      for ( i=0; i<ndefine; i++ ) {
        if ( !strcmp( define_words[i], str_tmp ) ) {
          using_define = 1;
          idefine = i;
          istring = 0;
        }
      }
    }
    if ( using_define ) {
      if      ( minus )
        strcpy( str, "-" );
      else
        strcpy( str, "" );
      strcat( str, define_strings[idefine][istring] );           
      istring++;
      if ( istring==define_nstring[idefine] ) {
        using_define = 0;
        istring = 0;
      }
    }
    else {
      for ( i=0; i<narithmetic; i++ ) {
        if ( !strcmp( arithmetic_words[i], str ) ) {
          d = arithmetic_values[i];
          d_is_set = 1;
        }
      }                 
    }
  }
  if ( !strcmp(str,"counter_a") ) {
    itoa(counter_a,str);
    counter_a++;
  }
  if ( !strcmp(str,"counter_a_apply") ) {
    itoa(counter_a,str);
  }
  if ( !strcmp(str,"counter_b") ) {
    itoa(counter_b,str);
    counter_b++;
  }
  if ( !strcmp(str,"counter_b_apply") ) {
    itoa(counter_b,str);
  }
  if ( !strcmp(str,"counter_c") ) {
    itoa(counter_c,str);
    counter_c++;
  }
  if ( !strcmp(str,"counter_c_apply") ) {
    itoa(counter_c,str);
    counter_c++;
  }
  if ( !strcmp(str,"counter_d") ) {
    itoa(counter_d,str);
    counter_d++;
  }
  if ( !strcmp(str,"counter_d_apply") ) {
    itoa(counter_d,str);
  }

}

void input_skip_comment( char str[] )

{
  long int i=0, l=0, nleft=0, nright=0, d_is_set=0;
  char *s, tmp_str[MCHAR];
  double d=0.;

  s = str;

  if ( *s=='(' ) {
    loop_comment:
    l = strlen(str);
    if ( cin.eof() ) {
      pri( "\nError in data file. Comment not closed." );
      exit(TN_EXIT_STATUS);
    }
    for ( i=0; i<=l-1; i++ ) {
      if ( *(s+i)=='(' ) nleft++;
      else if ( *(s+i)==')' ) nright++;
      if ( nright==nleft ) {
        if ( i==l-1 ) {
          if ( cin.eof() ) {
            pri( "\nError in data file. Comment not closed." );
            exit(TN_EXIT_STATUS);
          }
          else {
            input_read_string( 0, str, d, d_is_set );
              /* new part also starts with comment? */
            if ( str[0]=='(' ) goto loop_comment;
            string_convert_to_lower_case( str );
            return;
          }
        }
        else {
          strcpy( tmp_str, str );
          strcpy( str, &tmp_str[i+1] );
            /* new part also starts with comment? */
          if ( str[0]=='(' ) goto loop_comment;
          string_convert_to_lower_case( str );
          return;
        }
      }
    }
    if ( !(cin >> str) ) {
      pri( "\nError in data part." );
      pri( "Unexpected end of input (internal location h)." );
      pri( "Last word read ", str );
      exit(TN_EXIT_STATUS);
    }
    s = str;
    goto loop_comment;
  }

  string_convert_to_lower_case( str );
  if ( !strcmp(str,"counter_a") ) {
    itoa(counter_a,str);
    counter_a++;
  }
  if ( !strcmp(str,"counter_a_apply") ) {
    itoa(counter_a,str);
  }
  if ( !strcmp(str,"counter_b") ) {
    itoa(counter_b,str);
    counter_b++;
  }
  if ( !strcmp(str,"counter_b_apply") ) {
    itoa(counter_b,str);
  }
  if ( !strcmp(str,"counter_c") ) {
    itoa(counter_c,str);
    counter_c++;
  }
  if ( !strcmp(str,"counter_c_apply") ) {
    itoa(counter_c,str);
  }
  if ( !strcmp(str,"counter_d") ) {
    itoa(counter_d,str);
    counter_d++;
  }
  if ( !strcmp(str,"counter_d_apply") ) {
    itoa(counter_d,str);
  }


  return;

}

void input_runtime( void ) 

  /* 
    Read the runtime file at the start of each step.
    This allows for changing of data on the fly.
    E.g. use a control_print_gid record if you want 
    to make a plot in some long taking calculation.
  */

{
  long int iv=0, idat=0, index=0, length=0, itmp=0,
    last_data_value=0, ival[DATA_ITEM_SIZE];
  double dval[DATA_ITEM_SIZE];
  char filename[MCHAR], str[MCHAR];
  ifstream in;

  set_swit(-1,-1,"input_runtime");
  strcpy( filename, data_file_base );
  strcat( filename, ".run" );
  in.open( filename );
  if ( !in ) return;

  in >> str;
  string_convert_to_lower_case( str );
  while ( strcmp(str,"end_data" ) ) {
    idat = db_number( str ); 
    if ( idat<0 ) {
      goto close_runtime_file;
    }
    if ( db_no_index(idat ) ) 
       index = 0;
    else
       in >> index;
    length = last_data_value = 0;
    for ( iv=0; !last_data_value; iv++ ) {
      if ( !(in >> str) ) {
        goto close_runtime_file;
      }
      string_convert_to_lower_case( str );
      if      ( db_fixed_length(idat) ) {
        last_data_value = ( iv==(db_data_length(idat)) );
      }
      else if ( !strcmp(str,"end_data") )
        last_data_value = 1;
      else
        last_data_value = ( db_number(str)>=0 );
      if ( !last_data_value ) {
        if ( db_type(idat)==INTEGER ) {
          if ( str[0]=='-' ) {
            itmp = db_number( &str[1] );
            if ( itmp<0 ) {
              goto close_runtime_file;
            }
            check( itmp, CHECK_USAGE_AND_ERROR );
            ival[iv] = -itmp;
          }
          else {
            if ( !string_isinteger(str) ) {
              goto close_runtime_file;
            }
            ival[iv] = atoi(str);
          }
        }
        else if ( db_type(idat)==DOUBLE_PRECISION ) {
          if ( !string_isdouble(str) ) {
            goto close_runtime_file;
          }
          dval[iv] = atof(str);
        }
        length = iv + 1;
        if ( length>DATA_ITEM_SIZE ) {
          goto close_runtime_file;
        }
      }
    }
    db( idat, index, ival, dval, length, VERSION_NORMAL, PUT );
    any_runtime = 1;
  }

  close_runtime_file:
  in.close();
}

void input_check_required( void ) 

{
  long int idat=0, jdat=0, index=0, max_index=0, check_combination=-YES, ldum=0;
  double ddum[1];

  db( CHECK_COMBINATION, 0, &check_combination, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );

  if (  check_combination==-YES ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      jdat = db_data_required(idat);
      if ( jdat>=0 ) {
        db_max_index( idat, max_index, VERSION_NORMAL, GET );
        for ( index=0; index<=max_index; index++ ) {
          if ( db_active_index( idat, index, VERSION_NORMAL ) ) {
            if ( !db_active_index( jdat, index, VERSION_NORMAL ) ) {
              pri( "Error in data part" );
              pri( "Data item ", -idat );
              pri( "can only be used in combination with ", -jdat );
              pri( "To suppress this checking use CHECK_COMBINATION -NO in the data part" );
              exit(TN_EXIT_STATUS);
            }
          }
        }
      }
    }
  }
}
