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
#include <fstream>

#define MCHAR_WORDS 60
#define MSTRING 800
#define MDEFINE 200
#define MARITHMETIC 500
#define MLENGTH 10000

int counter_a=0, counter_b=0, counter_c=0, counter_d=0;
long int reading_define=0, using_define=0, idefine=0, ndefine=0, istring=0, define_nstring[MDEFINE];
std::ifstream include_file_stream;
long int include_reading=0;
long int input_abaqus_switch_global=0;
char *define_words[MDEFINE], *define_strings[MDEFINE][MSTRING];
long int reading_arithmetic=0, using_arithmetic=0, iarithmetic=0, narithmetic=0, using_if=0, using_if_not=0;
double arithmetic_values[MARITHMETIC];
char *arithmetic_words[MARITHMETIC];
// one-token pushback for the Professional no-index spelling
// "axisymmetric -yes"/-no (manual 6.17): the switch token is read where
// the parser expects the record index and is re-injected as the first
// data value (group 0), see the data part parser in input().
static char saved_first_value[MCHAR];
static long int has_saved_first_value = 0;
// mrange / mstring (manual Professional 4.55 / 4.56, initialization part):
// run-time limits below the compiled MRANGE (tochnog.h) and MSTRING (this
// file) buffer sizes. mrange caps the words of an index range -ra ... -ra
// and its expanded list; mstring caps the strings of a define block.
static long int mrange_limit = MRANGE, mrange_given = 0;
static long int mstring_limit = MSTRING-1, mstring_given = 0;

// Reads one whitespace-delimited word of at most MCHAR-1 characters into
// str (BACKLOG C56: a plain `in >> str` wrote past the MCHAR buffer for a
// longer word). Returns 0 at the end of the input, as `in >> str` did.
// A longer word is split: the rest stays in the stream and is read as the
// next word. That is harmless inside a comment (the parenthesis counting
// of input_skip_comment() sees every character), so a word that starts a
// comment (first character '(') or is read inside one (in_comment) may
// be long; any other over-long word is an input error.
static long int input_word( std::istream &in, char str[], long int in_comment )

{
  long int length=0;
  int c=0;

  if ( !(in >> std::setw(MCHAR) >> str) ) return 0;
  if ( in_comment || str[0]=='(' || in.eof() ) return 1;
  c = in.rdbuf()->sgetc();
  if ( c!=EOF && !isspace(c) ) {
    pri( "\nError in input file." );
    pri( "Word too long, it starts with ", str );
    length = MCHAR-1; pri( "maximum allowed length is ", length );
    exit(TN_EXIT_STATUS);
  }
  return 1;
}

// Reads the positive integer value of the initialization keyword name
// (mrange, mstring); compiled_max is the largest value the GNU buffers
// hold and constant names the constant to raise above it.
static long int input_initia_limit( const char *name, long int compiled_max,
  const char *constant )
{
  long int value=0;
  char str[MCHAR];
  if ( !input_word( cin, str, 0 ) || !string_isinteger( str ) ) {
    pri( "\nError in initialization part." );
    pri( "An integer value is expected after ", name );
    exit(TN_EXIT_STATUS);
  }
  value = atoi( str );
  if ( value<1 ) {
    pri( "\nError in initialization part." );
    pri( "The value should be larger than 0 for ", name );
    exit(TN_EXIT_STATUS);
  }
  if ( value>compiled_max ) {
    cout << "\nError in initialization part: " << name << " " << value
         << " exceeds the compiled maximum " << compiled_max << ".\n";
    cout << "Else, increase " << constant << " and recompile.\n";
    exit(TN_EXIT_STATUS);
  }
  return value;
}

// Index range length against mrange (Professional 02-08-2026: "Error:
// mrange exceeded."); only checked when mrange is specified, the
// compiled MRANGE limit is checked in range.cc.
static void input_mrange_check( long int length )
{
  if ( mrange_given && length>mrange_limit ) {
    pri( "\nError: mrange exceeded." );
    pri( "Increase it in the initialisation part." );
    exit(TN_EXIT_STATUS);
  }
}

void input( )

{
  long int in_value_range=0, value_range_words=0,
    i=0, j=0, n=0, iv=0, length=0, max=0, idat=0, itmp=0, idim=0,
    ready=0, d_is_set=0, index=0, range_length=0, last_data_value=0,
    unknown_indx=0, ninitia=0, istr=0, nstr=0,
    print_define=-NO, print_arithmetic=-NO,
    ldum=0, *integer_range=NULL, *range=NULL, 
    *ival=NULL, *dof_label=NULL, *dof_type=NULL, 
    *dof_principal=NULL, *dof_scal_vec_mat=NULL, 
    *dof_amount=NULL, *initialization_values=NULL;
  double d=0., ddum[1], *dval=NULL;
  char str_total[MCHAR], str[MCHAR], str_tmp[MCHAR];
  long int include_depth = 0;

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
    if ( !input_word( cin, str, 0 ) ) {
      pri( "\nError in initialization part." );
      exit(TN_EXIT_STATUS);
    }
    input_skip_comment( str ); if ( echo ) cout << str << " ";
    // Professional 2026 manual names of initialization records (4.14,
    // 4.28): normalised to the GNU name before it is stored, so every
    // check on the GNU name keeps working
    if ( !strcmp(str,"materi_history_variable") )
      strcpy( str, "materi_history_variables" );
    else if ( !strcmp(str,"materi_strain_energy") )
      strcpy( str, "materi_strainenergy" );
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
      if ( !input_word( cin, str_tmp, 0 ) ) {
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
    else if ( !strcmp(str,"groundflow_saturation") ) {
      groundflow_saturation = 1;
      gsat_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[gsat_indx], -GROUNDFLOW_SATURATION, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"groundflow_pressure_gradient") ) {
      // manual Professional 4.7: the gradient of the hydraulic pressure
      // head dh/dx dh/dy dh/dz is added to the node_dof records (names
      // pres_gradx..). The groundflow element (groundfl.cc) computes the
      // pressure gradient at the integration points; the dof is filled
      // there by recovery when declared.
      groundflow_pressure_gradient = 1;
      pres_grad_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[pres_grad_indx], -GROUNDFLOW_PRESSURE_GRADIENT,
        n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
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
    else if ( !strcmp(str,"materi_displacement_relative") ) {
      materi_displacement_relative = 1;
      dis_rel_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[dis_rel_indx], -MATERI_DISPLACEMENT_RELATIVE, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
    }
    else if ( !strcmp(str,"materi_history_variables") ||
              !strcmp(str,"materi_plasti_diprisco_history") ||
              !strcmp(str,"materi_plasti_hypo_history") ||
              !strcmp(str,"materi_plasti_camclay_history") ) {
      // materi_plasti_diprisco_history (manual Professional 4.18) is
      // the per-model name of materi_history_variables: same mechanism,
      // same shared hisv dof (basenames hisv0..hisv(n-1)). The manual
      // prescribes n = 11 for group_materi_plasti_diprisco and n = 12
      // for group_materi_plasti_diprisco_density.
      //
      // materi_plasti_hypo_history (manual Professional 4.23) declares
      // EIGHT hypoplasticity history variables named hyhis0..hyhis7
      // (void ratio, substep size, mobilized friction angle, stiffness
      // measure, structure s, OCR, density index, intergranular rho).
      // It uses the same shared hisv dof mechanism but with the hyhis
      // basename so target/print names match the Professional.
      if ( !strcmp(str,"materi_plasti_diprisco_history") )
        materi_plasti_diprisco_history = 1;
      if ( !strcmp(str,"materi_plasti_hypo_history") )
        materi_plasti_hypo_history = 1;
      if ( !strcmp(str,"materi_plasti_camclay_history") )
        materi_plasti_camclay_history = 1;
      if ( materi_plasti_camclay_history ) {
        // materi_plasti_camclay_history (manual Professional 4.16): the
        // history variables e0 (void ratio) and p0 (preconsolidation
        // pressure) of the camclay model are added to the node_dof
        // records with the basenames cchis0/cchis1. Two fixed history
        // variables; no number follows the keyword. The camclay law
        // (plasti.cc/stress.cc) reads them through the shared hisv
        // mechanism (hisv_indx = this dof, like the hypo hyhis dofs).
        materi_plasti_camclay_history = 1;
        materi_history_variables = 2;
        if ( echo ) cout << materi_history_variables;
        initialization_values[ninitia-1] = materi_history_variables;
        hisv_indx = unknown_indx;
        n = materi_history_variables;
        array_set( &dof_type[hisv_indx], -MATERI_PLASTI_CAMCLAY_HISTORY,
          materi_history_variables*nder );
        array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      }
      else if ( materi_plasti_hypo_history ) {
        // fixed 8 history variables, no number follows the keyword
        materi_history_variables = 8;
        if ( echo ) cout << materi_history_variables;
        initialization_values[ninitia-1] = materi_history_variables;
        hisv_indx = unknown_indx;
        n = materi_history_variables;
        array_set( &dof_type[hisv_indx], -MATERI_PLASTI_HYPO_HISTORY,
          materi_history_variables*nder );
        array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
      }
      else {
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
    }
    else if ( !strcmp(str,"materi_plasti_kappa_shear") ) {
      // manual Professional 4.25: the size of the shear part of the
      // plastic strain kappa_shear is added to the node_dof records.
      materi_plasti_kappa_shear = 1;
      kapsh_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[kapsh_indx], -MATERI_PLASTI_KAPPA_SHEAR, n*nder );
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
    else if ( !strcmp(str,"materi_strain_isa_c") ) {
      // Professional syntax layout: six symmetric back-strain components.
      // Registration is parse/storage only; no material law reads this DOF.
      n = 6;
      array_set( &dof_type[unknown_indx], -MATERI_STRAIN_ISA_C, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_isa_eacc") ) {
      // Professional syntax layout: one accumulated-strain scalar.
      // Registration is parse/storage only; no material law reads this DOF.
      n = 1;
      array_set( &dof_type[unknown_indx], -MATERI_STRAIN_ISA_EACC, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
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
    else if ( !strcmp(str,"materi_plasti_cap1_history") ) {
      materi_plasti_cap1_history = 1;
      cap1_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[cap1_indx], -MATERI_PLASTI_CAP1_HISTORY, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_plasti_hardsoil_history") ) {
      // manual Professional 4.22: the history variable abs(p) (maximum
      // pressure history) for the hardsoil model. Same concept as
      // materi_stress_pressure_history (4.50): the dof sph is SHARED
      // (same running max |p| update in dof.cc, same basename). The
      // loading/unloading switch of group_materi_elasti_hardsoil
      // (E50 vs Eur) reads this dof (old_unknowns[sph_indx]).
      materi_plasti_hardsoil_history = 1;
      sph_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[sph_indx], -MATERI_PLASTI_HARDSOIL_HISTORY, n*nder );
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
    else if ( !strcmp(str,"materi_plasti_phimob") ) {
      // F5 T3 C3: mobilized friction angle (manual Professional 4.26)
      materi_plasti_phimob = 1;
      phimob_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[phimob_indx], -MATERI_PLASTI_PHIMOB, n*nder );
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
    else if ( !strcmp(str,"materi_strain_plasti_hardsoil") ) {
      // manual Professional 4.40: the plastic strain specifically for
      // the hardsoil model, added to the node_dof records (same
      // 6-component layout as materi_strain_plasti, dedicated dof).
      // Filled like materi_strain_plasti: the dof integrates the
      // plastic strain increment (RHS in materi.cc).
      materi_strain_plasti_hardsoil = 1;
      hsepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[hsepp_indx], -MATERI_STRAIN_PLASTI_HARDSOIL, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_cap") ) {
      // manual Professional 4.35: the plastic strain specifically for
      // cap models, added to the node_dof records. Same mechanism as
      // materi_strain_plasti_hardsoil: dedicated dof, filled with the
      // plastic strain increment (RHS in materi.cc, consolidated block).
      materi_strain_plasti_cap = 1;
      capepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[capepp_indx], -MATERI_STRAIN_PLASTI_CAP, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_compression") ) {
      // manual Professional 4.36: idem for the compression model.
      materi_strain_plasti_compression = 1;
      cepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[cepp_indx], -MATERI_STRAIN_PLASTI_COMPRESSION, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_diprisco") ) {
      // manual Professional 4.37: idem for the di Prisco model.
      materi_strain_plasti_diprisco = 1;
      depp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[depp_indx], -MATERI_STRAIN_PLASTI_DIPRISCO, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_druckprag") ) {
      // manual Professional 4.39: idem for the drucker-prager model.
      materi_strain_plasti_druckprag = 1;
      dpepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[dpepp_indx], -MATERI_STRAIN_PLASTI_DRUCKPRAG, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_camclay") ) {
      // manual Professional 4.34: idem for the camclay model (F5 T2 V7).
      materi_strain_plasti_camclay = 1;
      camepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[camepp_indx], -MATERI_STRAIN_PLASTI_CAMCLAY, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_mohr_coul") ) {
      // manual Professional 4.41: idem for the mohr-coulomb model (F5 T2 V7).
      materi_strain_plasti_mohr_coul = 1;
      mcepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[mcepp_indx], -MATERI_STRAIN_PLASTI_MOHR_COUL, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_tension") ) {
      // manual Professional 4.42: idem for the tension model (F5 T2 V7).
      materi_strain_plasti_tension = 1;
      tnepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[tnepp_indx], -MATERI_STRAIN_PLASTI_TENSION, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_plasti_vonmises") ) {
      // manual Professional 4.43: idem for the von Mises model (F5 T2 V7).
      materi_strain_plasti_vonmises = 1;
      vmepp_indx = unknown_indx;
      n = 6;
      array_set( &dof_type[vmepp_indx], -MATERI_STRAIN_PLASTI_VONMISES, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -MATRIX, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_total_kappa") ) {
      // manual Professional 4.45: history maximum of a total strain
      // measure, added to the node_dof records (F5 T2 V7, materi.cc).
      materi_strain_total_kappa = 1;
      eptk_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[eptk_indx], -MATERI_STRAIN_TOTAL_KAPPA, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_total_compression_kappa") ) {
      // manual Professional 4.46: history maximum of a total strain
      // measure, added to the node_dof records (F5 T2 V7, materi.cc).
      materi_strain_total_compression_kappa = 1;
      eptkc_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[eptkc_indx], -MATERI_STRAIN_TOTAL_COMPRESSION_KAPPA, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_total_shear_kappa") ) {
      // manual Professional 4.47: history maximum of a total strain
      // measure, added to the node_dof records (F5 T2 V7, materi.cc).
      materi_strain_total_shear_kappa = 1;
      eptks_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[eptks_indx], -MATERI_STRAIN_TOTAL_SHEAR_KAPPA, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_strain_total_tension_kappa") ) {
      // manual Professional 4.48: history maximum of a total strain
      // measure, added to the node_dof records (F5 T2 V7, materi.cc).
      materi_strain_total_tension_kappa = 1;
      eptkt_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[eptkt_indx], -MATERI_STRAIN_TOTAL_TENSION_KAPPA, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
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
    else if ( !strcmp(str,"materi_stress_pressure_history") ) {
      // manual Professional 4.50: the maximum of the absolute value of
      // the pressure over time is stored in the node_dof records (see
      // group_materi_elasti_stress_pressure_history_factor, manual 6.655).
      materi_stress_pressure_history = 1;
      sph_indx = unknown_indx;
      n = 1;
      array_set( &dof_type[sph_indx], -MATERI_STRESS_PRESSURE_HISTORY, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
    }
    else if ( !strcmp(str,"materi_acceleration") ) {
      // manual Professional 4.11: the accelerations are added to the
      // node_dof records. The GNU dynamics scheme solves velocities;
      // the acceleration dofs are derived records updated at the end of
      // each step as a = (v_new - v_old)/dt (parallel_new_dof_diagonal
      // in dof.cc). Prescribing -accx in bounda_dof imposes the
      // velocity bound v_new = v_old + a*dt (bounda.cc).
      materi_acceleration = 1;
      acc_indx = unknown_indx;
      n = ndim;
      array_set( &dof_type[acc_indx], -MATERI_ACCELERATION, n*nder );
      array_set( &dof_scal_vec_mat[unknown_indx], -VECTOR, n*nder );
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
    else if ( !strcmp(str,"mrange") ) {
      mrange_limit = input_initia_limit( "mrange", MRANGE,
        "MRANGE in tochnog.h" );
      mrange_given = 1;
      if ( echo ) cout << mrange_limit;
      initialization_values[ninitia-1] = mrange_limit;
    }
    else if ( !strcmp(str,"mstring") ) {
      mstring_limit = input_initia_limit( "mstring", MSTRING-1,
        "MSTRING in input.cc" );
      mstring_given = 1;
      if ( echo ) cout << mstring_limit;
      initialization_values[ninitia-1] = mstring_limit;
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
  // F5 T2 V7: the total strain kappas are history maxima of the total
  // strain dof (materi.cc), which therefore needs to exist
  if ( materi_strain_total_kappa || materi_strain_total_compression_kappa ||
       materi_strain_total_shear_kappa || materi_strain_total_tension_kappa )
    check_unknown( "materi_strain_total", YES, CHECK_USAGE_AND_ERROR );
  if ( materi_displacement_relative ) {
    check_unknown( "materi_displacement", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_velocity", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_velocity_integrated", YES, CHECK_USAGE_AND_ERROR );
  }
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
  if ( materi_plasti_phimob )
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
  if ( materi_strainenergy ) {
    check_unknown( "materi_stress", YES, CHECK_USAGE_AND_ERROR );
    check_unknown( "materi_strain_elasti", YES, CHECK_USAGE_AND_ERROR );
  }


    /* read data part */
  input_read_string( echo, str, d, d_is_set );
  input_skip_comment( str );
  while ( strcmp(str,"end_data") ) {

      // if the end_data of an included file is reached, restore the
      // main input stream and continue with the remaining records
    if ( !strcmp(str,"end_data") && include_reading ) {
      include_file_stream.close();
      include_reading = 0;
      include_depth = 0;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      continue;
    }

      // data item name
    input_skip_comment( str );
    if ( echo ) cout << str << " ";
    idat = db_number( str );
    if ( idat<0 ) {
      pri( "\nError in data part." );
      pri( "I do not know ", str );
      exit(TN_EXIT_STATUS);
    }

    // data_ignore (manual Professional 6.401): skip every record with
    // the listed item name. The record loop condition re-tests str, so
    // after consuming the ignored values 'continue' processes the next
    // keyword and "end_data" exits the loop naturally.
    if ( db_active_index( DATA_IGNORE, 0, VERSION_NORMAL ) ) {
      long int iign=0, max_ign=0, ign_item=0, skipping=0, ldum2=0;
      double ddum2[1];
      db_max_index( DATA_IGNORE, max_ign, VERSION_NORMAL, GET );
      for ( iign=0; iign<=max_ign; iign++ ) {
        if ( db_active_index( DATA_IGNORE, iign, VERSION_NORMAL ) ) {
          ign_item = -1;
          db( DATA_IGNORE, iign, &ign_item, ddum2, ldum2,
            VERSION_NORMAL, GET );
          if ( ign_item==-idat || ign_item==idat ) skipping = 1;
        }
      }
      if ( skipping ) {
        if ( echo ) cout << "\n";
        while ( strcmp(str,"end_data") ) {
          input_read_string( echo, str, d, d_is_set );
          input_skip_comment( str );
          if ( d_is_set ) continue;
          if ( db_number(str)>=0 || !strcmp(str,"end_data") ) break;
        }
        continue;
      }
    }

    // import an abaqus mesh (generates tochnog_abaqus.dat)
    if ( idat==INPUT_ABAQUS ) {
      long int input_abaqus_switch=0;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      if ( echo ) cout << " " << str << " ";
      if ( str[0]=='-' ) {
        itmp = db_number( &str[1] );
        if ( itmp<0 ) {
          pri( "\nError in data part." );
          pri( "I do not know ", str );
          exit(TN_EXIT_STATUS);
        }
        input_abaqus_switch = -itmp;
      }
      input_abaqus_switch_global = input_abaqus_switch;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      continue;
    }

    // after all input_abaqus_* sub-options, generate tochnog_abaqus.dat
    if ( idat==INPUT_ABAQUS_CONTINUE ) {
      long int input_abaqus_continue=0;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      if ( echo ) cout << " " << str << " ";
      if ( str[0]=='-' ) {
        itmp = db_number( &str[1] );
        if ( itmp<0 ) {
          pri( "\nError in data part." );
          pri( "I do not know ", str );
          exit(TN_EXIT_STATUS);
        }
        input_abaqus_continue = -itmp;
      }
      if ( input_abaqus_switch_global==-YES && input_abaqus_continue==-YES )
        input_abaqus_read();
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      continue;
    }

    // import a feflow mesh
    if ( idat==INPUT_FEFLOW_MESH ) {
      long int input_feflow_switch=0;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      if ( echo ) cout << " " << str << " ";
      if ( str[0]=='-' ) {
        itmp = db_number( &str[1] );
        if ( itmp<0 ) {
          pri( "\nError in data part." );
          pri( "I do not know ", str );
          exit(TN_EXIT_STATUS);
        }
        input_feflow_switch = -itmp;
      }
      if ( input_feflow_switch==-YES ) input_feflow_read();
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      continue;
    }

    // import a gmsh mesh
    if ( idat==INPUT_GMSH ) {
      long int input_gmsh_switch=0;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      if ( echo ) cout << " " << str << " ";
      if ( str[0]=='-' ) {
        itmp = db_number( &str[1] );
        if ( itmp<0 ) {
          pri( "\nError in data part." );
          pri( "I do not know ", str );
          exit(TN_EXIT_STATUS);
        }
        input_gmsh_switch = -itmp;
      }
      if ( input_gmsh_switch==-YES ) input_gmsh_read();
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      continue;
    }

    // include a data file
    if ( idat==INCLUDE ) {
      if ( include_depth>=1 ) {
        pri( "\nError in data part." );
        pri( "include files cannot contain an include." );
        exit(TN_EXIT_STATUS);
      }
      // read the file name
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      if ( echo ) cout << " " << str << " ";
      include_file_stream.open( str );
      if ( !include_file_stream ) {
        pri( "\nError in data part." );
        pri( "Cannot open include file ", str );
        exit(TN_EXIT_STATUS);
      }
      include_reading = 1;
      include_depth = 1;
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
      continue;
    }

    check( idat, CHECK_USAGE_AND_ERROR );

      // index
    if ( db_no_index(idat ) )
      range[0] = 0;
    else {
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment(str);
      if ( echo ) cout << " " << str << " ";
      if ( !strcmp(str,"-ra") ) {
         range_scan( echo, d, d_is_set, range, range_length );
         // the words of the range with its two -ra markers (Professional
         // 02-08-2026: a 3-entry index range needs mrange 5)
         input_mrange_check( range_length );
      }
      else {
        if ( !string_isinteger(str) ) {
          // Professional spelling "axisymmetric -yes/-no" WITHOUT an
          // index (manual 6.17: the switch makes the whole calculation
          // axi-symmetrical; a per-group group_axisymmetric overrules
          // it). The GNU data model is per-group: the switch is stored
          // on group 0 (the default group of the elements). The -yes/-no
          // token becomes the first data value of the record.
          if ( idat==GROUP_AXISYMMETRIC && str[0]=='-' ) {
            strcpy( saved_first_value, str );
            has_saved_first_value = 1;
            range[0] = 0;
          }
          else {
            pri( "\nError in data part." );
            pri( "Illegal index ", str );
            exit(TN_EXIT_STATUS);
          }
        }
        else
          range[0] = atoi(str);
      }
    }


      // data item values
    length = last_data_value = istr = nstr = in_value_range = 0;
    for ( iv=0; !last_data_value; iv++ ) {
      if ( istr>=nstr ) {
        if ( has_saved_first_value && iv==0 ) {
          // re-inject the index-position token (axisymmetric -yes/-no)
          strcpy( str_total, saved_first_value );
          has_saved_first_value = 0;
        }
        else
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
      else if ( !strcmp(str,"start_if") || !strcmp(str,"start_if_not") ||
                !strcmp(str,"end_if") || !strcmp(str,"end_if_not") )
        last_data_value = 1;
      else if ( !d_is_set ) 
       last_data_value = ( db_number(str)>=0 );
      if ( !last_data_value ) {
        if ( echo ) {
          if ( d_is_set ) cout << " " << d << " ";
          else cout << " " << str << " ";
        }
        if ( idat==POST_CALCUL_LABEL ) {
            // post_calcul_label (manual Professional 6.904) is an output
            // record (print_db.cc writes the item names of post_calcul):
            // its names are read and ignored, so a database file holding
            // the record can be read back
          ival[iv] = 0;
        }
        else if ( db_type(idat)==INTEGER ) {
          if ( !d_is_set && str[0]=='-' &&
               idat==CONTROL_MESH_EXTRUDE_ELEMENT_GROUP_NEW &&
               str[1]!='\0' && string_isinteger( &str[1] ) ) {
            // manual Professional 6.184: a negative new element group
            // drops the extruded elements (other/extrude1: 1 1 -1)
            ival[iv] = -atoi( &str[1] );
          }
          else if ( str[0]=='-' ) {
            itmp = db_number( &str[1] );
            if ( itmp<0 && idat==GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION &&
                 !strcmp(str,"-both") ) {
              // C80.2 (2026-10-16): the C80.1 value -both (shear terms
              // and volumetric part at the centroid point) was removed
              // (never better than -volumetric, too soft in coarse
              // bending).
              pri( "\nError in data part." );
              pri( "group_element_selective_reduced_integration -both was "
                "removed; use -lambda (elastic, bending-dominated) or "
                "-volumetric (general) instead." );
              exit(TN_EXIT_STATUS);
            }
            if ( itmp<0 ) {
              pri( "\nError in data part." );
              pri( "I do not know ", str );
              exit(TN_EXIT_STATUS);
            }
            check( itmp, CHECK_USAGE_AND_ERROR );
            ival[iv] = -itmp;
            if ( itmp==RA ) {
              // mrange: the words of the value range just closed
              if ( in_value_range ) input_mrange_check( value_range_words );
              in_value_range = !in_value_range;
              value_range_words = 0;
            }
            else if ( in_value_range ) value_range_words++;
          }
          else {
            if ( d_is_set )
              ival[iv] = (long int) d;
            else {
                // a number in a value range -ra ... -ra must be an
                // integer (Professional: 9.0 is rejected; F4 T6f)
              if ( ( in_value_range && !string_isinteger(str) ) ||
                   ( !string_isinteger(str) && !string_isdouble(str) ) ) {
                pri( "\n\nError in data part." );
                pri( "Problem reading ", db_name(idat) );
                pri( "I don't know what to do with ", str );
                if ( db_fixed_length(idat) )
                  pri( "Number of data values expected ", db_data_length(idat) );
                exit(TN_EXIT_STATUS);
              }
              ival[iv] = atoi(str);
            }
            if ( in_value_range ) value_range_words++;
          }
        }
        else if ( db_type(idat)==DOUBLE_PRECISION ) {
          if ( d_is_set )
            dval[iv] = d;
          else {
            if ( !string_isdouble(str) ) {
              // force_element_edge_water (manual Professional 6.489):
              // the Professional record is a pure switch
              // (force_edge_water index -yes). A DOUBLE record cannot
              // hold the token through the generic path: store the
              // negative YES/NO enum as the single value, the layout
              // discriminator used by area() (legacy GNU: rho g dirs).
              if ( idat==FORCE_ELEMENT_EDGE_WATER && str[0]=='-' ) {
                long int swit_word = db_number( &str[1] );
                if ( swit_word==YES || swit_word==NO ) {
                  dval[iv] = -(double)swit_word;
                }
                else {
                  pri( "\n\nError in data part." );
                  pri( "Problem reading ", db_name(idat) );
                  pri( "I don't know what to do with ", str );
                  exit(TN_EXIT_STATUS);
                }
              }
              else {
                pri( "\n\nError in data part." );
                pri( "Problem reading ", db_name(idat) );
                pri( "I don't know what to do with ", str );
                if ( db_fixed_length(idat) )
                  pri( "Number of data values expected ", db_data_length(idat) );
                exit(TN_EXIT_STATUS);
              }
            }
            else
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

      // post_element_force (manual Professional 6.927, C84): the record
      // length depends on the space dimensions, 3D 12 values, 2D 6
      // (dir_normal_x dir_normal_y dir_shear0_x dir_shear0_y middle_x
      // middle_y), 1D 2 (dir_normal_x middle_x); the Professional
      // 02-08-2026 rejects any other count in the data part
    if ( idat==POST_ELEMENT_FORCE ) {
      long int pef_expected = ( ndim==3 ? 12 : ( ndim==2 ? 6 : 2 ) );
      if ( length!=pef_expected ) {
        pri( "\n\nError in data part." );
        pri( "Problem reading ", db_name(idat) );
        pri( "Number of data values expected ", pef_expected );
        pri( "Number of data values found ", length );
        exit(TN_EXIT_STATUS);
      }
    }

      // store
    if ( range[0]==-RA ) {
      range_expand( range, integer_range, ldum, range_length );
      input_mrange_check( range_length );
    }
    db_max_index( idat, max, VERSION_NORMAL, GET );
    if ( idat==NODE_DOF ) {
    }
      // post_calcul_label is not stored (output record, see above)
    for ( i=0, ready=(idat==POST_CALCUL_LABEL); !ready; i++) {
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

  // groundflow_phreatic_only -yes: drop the groundflow data of the groups
  // without a multiple phreatic level, before the calculation starts.
  groundflow_phreatic_only_apply();

  // strain_volume_*: table checks and initial volumes of the sets.
  strain_volume_setup();

  // check that elements do not have duplicate nodes (check_element_node)
  {
    long int check_element_node_switch=-YES;
    db( CHECK_ELEMENT_NODE, 0, &check_element_node_switch, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    check_element_node( check_element_node_switch );
  }

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
    // dof_element_dof (manual Professional 6.409, print only): -yes for the
    // dofs an element remembers in its integration points by default
    // (non-continuous field: the element_intpnt_dof restore set of elem.cc,
    // integration-point states of the constitutive update), -no for the
    // continuous fields (principal dofs, displacement and every dof that
    // elem.cc interpolates from the nodes). The record shows the default
    // classification; global_element_dof_apply -no does not change it (as
    // the Professional 02-08-2026, probes f5/v5b a_no).
    long int *dof_element_dof = get_new_int( nuknwn );
    for ( long int iu=0; iu<nuknwn; iu++ ) {
      long int t = -dof_type[iu];
      dof_element_dof[iu] = ( t==MATERI_STRESS || t==MATERI_STRAIN_ELASTI ||
        t==MATERI_STRAIN_PLASTI || t==MATERI_STRAIN_TOTAL ||
        t==MATERI_STRAIN_INTERGRANULAR || t==MATERI_HISTORY_VARIABLES ||
        t==MATERI_PLASTI_CAMCLAY_HISTORY || t==MATERI_PLASTI_HYPO_HISTORY ||
        t==MATERI_PLASTI_KAPPA || t==MATERI_PLASTI_KAPPA_SHEAR ||
        t==MATERI_PLASTI_CAP1_HISTORY || t==MATERI_DAMAGE ||
        t==MATERI_PLASTI_RHO || t==MATERI_MAXWELL_STRESS ||
        t==MATERI_PLASTI_SOFTVAR_LOCAL || t==MATERI_PLASTI_SOFTVAR_NONLOCAL ) ?
        -YES : -NO;
    }
    db( DOF_ELEMENT_DOF, 0, dof_element_dof, ddum, nuknwn, VERSION_NORMAL, PUT );
    delete[] dof_element_dof;
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
    if ( include_reading ) {
      if ( !input_word( include_file_stream, str, 0 ) ) {
        pri( "\nError in data part." );
        pri( "Unexpected end of include file detected." );
        exit(TN_EXIT_STATUS);
      }
    }
    else if ( !input_word( cin, str, 0 ) ) {
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
    if ( !input_word( cin, str, 0 ) ) {
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
      if ( !input_word( cin, str, 0 ) ) {
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
        if ( istring==mstring_limit ) {
          cout << "\nError: maximal number of strings is "
               << ( mstring_given ? mstring_limit : MSTRING ) << "\n";
          pri( "Error detected for define: ", define_words[ndefine] );
          pri( "Maybe the define is not correct." );
          if ( mstring_given )
            pri( "Else, increase mstring in the initialisation part." );
          else
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
    if ( !input_word( cin, str, 0 ) ) {
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
    if ( !input_word( cin, str, 0 ) ) {
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
      if ( !input_word( cin, str, 0 ) ) {
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
  else if ( !strcmp(str,"start_if") || !strcmp(str,"start_if_not") ) {
    // conditional blocks start_if ... end_if and
    // start_if_not ... end_if_not (manual Professional 5.x): the
    // records of the block are applied only when the start_define
    // word is set to true (start_if) or false (start_if_not).
    if ( using_if ) {
      pri( "Error, start_if cannot be nested." );
      exit(TN_EXIT_STATUS);
    }
    using_if = 1;
    using_if_not = ( str[9]=='n' );
    if ( !input_word( cin, str, 0 ) ) {
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
    if ( using_if_not ) apply_if = !apply_if;
    if ( apply_if ) {
      input_read_string( echo, str, d, d_is_set );
      input_skip_comment( str );
    }
    else {
      loop_if:
      if ( !input_word( cin, str, 0 ) ) {
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
      if ( strcmp(str,"end_if") && strcmp(str,"end_if_not") ) {
        goto loop_if;
      }
      else {
        using_if = 0;
        using_if_not = 0;
        input_read_string( echo, str, d, d_is_set );
        input_skip_comment( str );
      }
    }
  }
  else if ( ( !strcmp(str,"end_if") || !strcmp(str,"end_if_not") ) &&
            using_if ) {
    using_if = 0;
    using_if_not = 0;
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
    if ( !input_word( cin, str, 1 ) ) {
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

// Reads simple data records "name index values ... end_data" from in:
// no comments, defines, arithmetic or ranges (the reader of the runtime
// file and of control_input). strict 0 (runtime file): stop silently at
// the first word that does not fit, a record is stored as soon as it is
// complete. strict 1 (control_input, manual Professional 6.139): such a
// word is an error, as in the Professional 02-08-2026, and every record
// is checked as in the data part. mesh_read is set when a node or
// element record was read.
static void input_simple_records( std::istream &in, long int strict,
  long int &mesh_read )

{
  long int iv=0, idat=0, index=0, length=0, itmp=0,
    last_data_value=0, ival[DATA_ITEM_SIZE];
  double dval[DATA_ITEM_SIZE];
  char str[MCHAR];

  strcpy( str, "" );
  if ( !input_word( in, str, 0 ) ) goto simple_records_error;
  string_convert_to_lower_case( str );
  while ( strcmp(str,"end_data" ) ) {
    idat = db_number( str ); 
    if ( idat<0 ) {
      goto simple_records_error;
    }
    if ( strict ) check( idat, CHECK_USAGE_AND_ERROR );
    if ( db_no_index(idat ) ) 
       index = 0;
    else if ( strict ) {
      if ( !input_word( in, str, 0 ) || !string_isinteger( str ) )
        goto simple_records_error;
      index = atoi( str );
    }
    else
       in >> index;
    length = last_data_value = 0;
    for ( iv=0; !last_data_value; iv++ ) {
      if ( !input_word( in, str, 0 ) ) {
        goto simple_records_error;
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
              goto simple_records_error;
            }
            check( itmp, CHECK_USAGE_AND_ERROR );
            ival[iv] = -itmp;
          }
          else {
            if ( !string_isinteger(str) ) {
              goto simple_records_error;
            }
            ival[iv] = atoi(str);
          }
        }
        else if ( db_type(idat)==DOUBLE_PRECISION ) {
          if ( !string_isdouble(str) ) {
            goto simple_records_error;
          }
          dval[iv] = atof(str);
        }
        length = iv + 1;
        if ( length>DATA_ITEM_SIZE ) {
          goto simple_records_error;
        }
      }
    }
    db( idat, index, ival, dval, length, VERSION_NORMAL, PUT );
    if ( idat==NODE || idat==ELEMENT ) mesh_read = 1;
    if ( !strict ) any_runtime = 1;
  }
  return;

  simple_records_error:
  if ( strict ) {
    pri( "\nError for control_input detected." );
    cout << "Error detected near: " << str << "\n";
    exit(TN_EXIT_STATUS);
  }
}

void input_runtime( void ) 

  /* 
    Read the runtime file at the start of each step.
    This allows for changing of data on the fly.
    E.g. use a control_print_gid record if you want 
    to make a plot in some long taking calculation.
  */

{
  long int mesh_read=0;
  char filename[MCHAR];
  ifstream in;

  set_swit(-1,-1,"input_runtime");
  strcpy( filename, data_file_base );
  strcat( filename, ".run" );
  in.open( filename );
  if ( !in ) return;

  input_simple_records( in, 0, mesh_read );
  in.close();
}

// control_input index -yes (manual Professional 6.139): at the start of
// control index icontrol read the data records of the file
// <icontrol>.dat (no comments, defines, arithmetic or ranges; closed by
// end_data, the rest of the file is ignored). Measured on the
// Professional 02-08-2026: a record replaces an existing one with the
// same index, new records are added (also nodes and elements: the mesh
// is extended), a missing file is silently skipped.
void input_control_input( long int icontrol )

{
  long int control_input=-NO, ldum=0, mesh_read=0, inod=0, max_node=-1,
    length=0, idum[1];
  double ddum[1];
  char filename[MCHAR];
  ifstream in;

  db( CONTROL_INPUT, icontrol, &control_input, ddum, ldum, VERSION_NORMAL,
    GET_IF_EXISTS );
  if ( control_input!=-YES ) return;
  long_to_a( icontrol, filename );
  strcat( filename, ".dat" );
  in.open( filename );
  if ( !in ) return;
  input_simple_records( in, 1, mesh_read );
  in.close();
  if ( mesh_read ) {
      // new nodes get zero dofs, as at the start of the calculation (top)
    if ( nuknwn>0 ) {
      double *dwork = get_new_dbl( nuknwn );
      array_set( dwork, 0., nuknwn );
      db_max_index( NODE, max_node, VERSION_NORMAL, GET );
      for ( inod=0; inod<=max_node; inod++ ) {
        if ( !db_active_index( NODE, inod, VERSION_NORMAL ) ) continue;
        if ( !db_active_index( NODE_DOF, inod, VERSION_NORMAL ) )
          db( NODE_DOF, inod, idum, dwork, nuknwn, VERSION_NORMAL, PUT );
          // the start records top() copies before the control loop
        if ( !db_active_index( NODE_START_REFINED, inod, VERSION_NORMAL ) ) {
          length = db_len( NODE, inod, VERSION_NORMAL );
          db( NODE_START_REFINED, inod, idum,
            db_dbl( NODE, inod, VERSION_NORMAL ), length, VERSION_NORMAL, PUT );
        }
        if ( !db_active_index( NODE_DOF_START_REFINED, inod, VERSION_NORMAL ) ) {
          length = db_len( NODE_DOF, inod, VERSION_NORMAL );
          db( NODE_DOF_START_REFINED, inod, idum,
            db_dbl( NODE_DOF, inod, VERSION_NORMAL ), length, VERSION_NORMAL,
            PUT );
        }
      }
      delete[] dwork;
        // new elements get the element records top() initializes
      long int ielem=0, max_elem=-1, zero=0, one=1,
        mnolnuknwn=npointmax*nuknwn, length_nei=1+npointmax*ndim+npointmax+2;
      double *dworkmnol = get_new_dbl( mnolnuknwn );
      double *dworknei = get_new_dbl( length_nei );
      array_set( dworkmnol, 0., mnolnuknwn );
      array_set( dworknei, 0., length_nei );
      dworknei[length_nei-1] = 1;
      db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET );
      for ( ielem=0; ielem<=max_elem; ielem++ ) {
        if ( !db_active_index( ELEMENT, ielem, VERSION_NORMAL ) ) continue;
        if ( !db_active_index( ELEMENT_INTPNT_DOF, ielem, VERSION_NORMAL ) )
          db( ELEMENT_INTPNT_DOF, ielem, idum, dworkmnol, mnolnuknwn,
            VERSION_NORMAL, PUT );
        if ( !db_active_index( ELEMENT_DOF_INITIALISED, ielem,
            VERSION_NORMAL ) )
          db( ELEMENT_DOF_INITIALISED, ielem, &zero, ddum, one,
            VERSION_NORMAL, PUT );
        if ( !db_active_index( NONLOCAL_ELEMENT_INFO, ielem,
            VERSION_NORMAL ) )
          db( NONLOCAL_ELEMENT_INFO, ielem, idum, dworknei, length_nei,
            VERSION_NORMAL, PUT );
      }
      delete[] dworkmnol;
      delete[] dworknei;
    }
    mesh_has_changed( VERSION_NORMAL );
  }
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

void input_gmsh_read( void )

{
  // read the gmsh mesh file tochnog_in.msh (format 2.2)
  // only node, element and element_group are read.
  long int i=0, inode=0, ielem=0, nnode=0, nelem=0, ieltype=0, ntag=0,
    igeom=0, node_tag=0, nnol=0, inol=0, idum[1], 
    elem_data[1+MNOL];
  double ddum[1], xyz[MDIM];
  char str[MCHAR], filename[MCHAR];
  ifstream in;

  strcpy( filename, "tochnog_in.msh" );
  in.open( filename );
  if ( !in ) {
    pri( "Error: cannot open gmsh file ", filename );
    exit(TN_EXIT_STATUS);
  }

  // read MeshFormat section
  while ( input_word( in, str, 1 ) ) {   // skipped words may be long
    if ( !strcmp(str,"$Nodes") ) break;
  }
  in >> nnode;
  for ( inode=1; inode<=nnode; inode++ ) {
    in >> node_tag >> xyz[0] >> xyz[1] >> xyz[2];
    db( NODE, node_tag, idum, xyz, ndim, VERSION_NORMAL, PUT );
  }

  // read Elements section
  while ( input_word( in, str, 1 ) ) {   // skipped words may be long
    if ( !strcmp(str,"$Elements") ) break;
  }
  in >> nelem;
  for ( ielem=1; ielem<=nelem; ielem++ ) {
    in >> igeom >> ieltype >> ntag;
    for ( i=0; i<ntag; i++ ) in >> idum[0];  // skip tags
    // element type -> tochnog element name and number of nodes
    if      ( ieltype==1 ) { elem_data[0] = -BAR2;  nnol = 2; }
    else if ( ieltype==8 ) { elem_data[0] = -BAR3;  nnol = 3; }
    else if ( ieltype==2 ) { elem_data[0] = -TRIA3; nnol = 3; }
    else if ( ieltype==9 ) { elem_data[0] = -TRIA6; nnol = 6; }
    else if ( ieltype==3 ) { elem_data[0] = -QUAD4; nnol = 4; }
    else if ( ieltype==10 ){ elem_data[0] = -QUAD9; nnol = 9; }
    else {
      cout << "Error: gmsh element type " << ieltype << " not supported by tochnog.\n";
      exit(TN_EXIT_STATUS);
    }
    for ( inol=0; inol<nnol; inol++ ) in >> elem_data[1+inol];
    ntag = 1+nnol;
    db( ELEMENT, igeom, elem_data, ddum, ntag, VERSION_NORMAL, PUT );
    // element group from gmsh physical group (first tag)
  }

  in.close();
  cout << "Input from gmsh file " << filename << " read.\n";

}

void input_abaqus_read( void )

{
  // read the abaqus input file abaqus.inp and generate tochnog_abaqus.dat
  // with node and element records. input_abaqus_name limits the element
  // types that are converted. input_abaqus_set limits the element numbers.
  // input_abaqus_group decides whether group_* records are written.
  long int i=0, inode=0, nnol=0, n1=0, n2=0, n3=0, n4=0, n5=0,
    n6=0, n7=0, n8=0, n9=0, elem_id=0, nn=0, *input_abaqus_name=NULL,
    *input_abaqus_set=NULL, nset_=0, in_set=0, j=0,
    input_abaqus_group=-YES, ldum=0;
  double ddum[1], xyz[MDIM];
  char line[MCHAR], *tok=NULL, str2[MCHAR], filename[MCHAR], eltype[MCHAR];
  ifstream in;
  ofstream out;
  long int in_nodes=0, in_elements=0, convert_ok=0, in_elastic=0;
  double young=0., poisson=0.;

  // element types to convert (input_abaqus_name); empty means all
  if ( db_active_index( INPUT_ABAQUS_NAME, 0, VERSION_NORMAL ) ) {
    nn = db_len( INPUT_ABAQUS_NAME, 0, VERSION_NORMAL );
    if ( nn>0 ) {
      input_abaqus_name = get_new_int(nn);
      db( INPUT_ABAQUS_NAME, 0, input_abaqus_name, ddum, nn, VERSION_NORMAL, GET );
    }
  }

  // element numbers to write (input_abaqus_set); empty means all
  if ( db_active_index( INPUT_ABAQUS_SET, 0, VERSION_NORMAL ) ) {
    nset_ = db_len( INPUT_ABAQUS_SET, 0, VERSION_NORMAL );
    if ( nset_>0 ) {
      input_abaqus_set = get_new_int(nset_);
      db( INPUT_ABAQUS_SET, 0, input_abaqus_set, ddum, nset_, VERSION_NORMAL, GET );
    }
  }

  // whether group_* records are written
  db( INPUT_ABAQUS_GROUP, 0, &input_abaqus_group, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );

  strcpy( filename, "abaqus.inp" );
  in.open( filename );
  if ( !in ) {
    pri( "Error: cannot open abaqus file ", filename );
    exit(TN_EXIT_STATUS);
  }
  out.open( "tochnog_abaqus.dat" );

  in_nodes = in_elements = 0;
  while ( in.getline(line,MCHAR) ) {
    string_convert_to_lower_case( line );
    // trim leading spaces
    char *p = line;
    while ( *p==' ' || *p=='\t' ) p++;
    if ( p[0]=='*' ) {
      in_nodes = in_elements = in_elastic = 0;
      // strip trailing comma of keyword
      char *q = p;
      while ( *q && *q!='\n' ) q++;
      *q = 0;
      // keyword
      if      ( !strncmp(p,"*node",5) ) in_nodes = 1;
      else if ( !strncmp(p,"*element",8) ) {
        in_elements = 1;
        strcpy(eltype,"");
        char *eq = strstr(p,"type=");
        if ( eq ) { strcpy(eltype, eq+5); char *c=strchr(eltype,','); if(c)*c=0; }
      }
      else if ( !strncmp(p,"*elastic",8) ) in_elastic = 1;
      continue;
    }
    // data line: split by commas
    if ( in_elastic ) {
      tok = strtok(p,",");
      if ( tok ) young = atof(tok);
      tok = strtok(NULL,",");
      if ( tok ) poisson = atof(tok);
    }
    else if ( in_nodes ) {
      tok = strtok(p,",");
      if ( !tok ) continue;
      inode = atoi(tok);
      xyz[0]=xyz[1]=xyz[2]=0.;
      i=0;
      tok = strtok(NULL,",");
      while ( tok && i<3 ) { xyz[i]=atof(tok); tok=strtok(NULL,","); i++; }
      out << "node  " << inode << "  " << xyz[0];
      if ( ndim>=2 ) out << "  " << xyz[1];
      if ( ndim==3 ) out << "  " << xyz[2];
      out << "\n";
    }
    else if ( in_elements ) {
      tok = strtok(p,",");
      if ( !tok ) continue;
      elem_id = atoi(tok);
      n1=n2=n3=n4=n5=n6=n7=n8=n9=0;
      if      ( !strncmp(eltype,"t2d2",4) || !strncmp(eltype,"t3d2",4) ) {
        strcpy(str2,"-bar2"); nnol=2;
        tok=strtok(NULL,","); n1=atoi(tok?tok:"");
        tok=strtok(NULL,","); n2=atoi(tok?tok:"");
      }
      else if ( !strncmp(eltype,"t2d3",4) || !strncmp(eltype,"t3d3",4) ) {
        strcpy(str2,"-bar3"); nnol=3;
        tok=strtok(NULL,","); n1=atoi(tok?tok:"");
        tok=strtok(NULL,","); n2=atoi(tok?tok:"");
        tok=strtok(NULL,","); n3=atoi(tok?tok:"");
      }
      else if ( !strncmp(eltype,"cps3",4) || !strncmp(eltype,"cpe3",4) ||
                !strncmp(eltype,"cax3",4) || !strncmp(eltype,"s3",2) ) {
        strcpy(str2,"-tria3"); nnol=3;
        tok=strtok(NULL,","); n1=atoi(tok?tok:"");
        tok=strtok(NULL,","); n2=atoi(tok?tok:"");
        tok=strtok(NULL,","); n3=atoi(tok?tok:"");
      }
      else if ( !strncmp(eltype,"cps6",4) || !strncmp(eltype,"cpe6",4) ) {
        strcpy(str2,"-tria6"); nnol=6;
        tok=strtok(NULL,","); n1=atoi(tok?tok:"");
        tok=strtok(NULL,","); n2=atoi(tok?tok:"");
        tok=strtok(NULL,","); n3=atoi(tok?tok:"");
        tok=strtok(NULL,","); n4=atoi(tok?tok:"");
        tok=strtok(NULL,","); n5=atoi(tok?tok:"");
        tok=strtok(NULL,","); n6=atoi(tok?tok:"");
      }
      else if ( !strncmp(eltype,"cps4",4) || !strncmp(eltype,"cpe4",4) ||
                !strncmp(eltype,"cax4",4) || !strncmp(eltype,"m3d4",4) ) {
        strcpy(str2,"-quad4"); nnol=4;
        tok=strtok(NULL,","); n1=atoi(tok?tok:"");
        tok=strtok(NULL,","); n2=atoi(tok?tok:"");
        tok=strtok(NULL,","); n3=atoi(tok?tok:"");
        tok=strtok(NULL,","); n4=atoi(tok?tok:"");
      }
      else if ( !strncmp(eltype,"cps8",4) || !strncmp(eltype,"cpe8",4) ||
                !strncmp(eltype,"cax8",4) || !strncmp(eltype,"m3d8",4) ) {
        strcpy(str2,"-quad9"); nnol=8;
        tok=strtok(NULL,","); n1=atoi(tok?tok:"");
        tok=strtok(NULL,","); n2=atoi(tok?tok:"");
        tok=strtok(NULL,","); n3=atoi(tok?tok:"");
        tok=strtok(NULL,","); n4=atoi(tok?tok:"");
        tok=strtok(NULL,","); n5=atoi(tok?tok:"");
        tok=strtok(NULL,","); n6=atoi(tok?tok:"");
        tok=strtok(NULL,","); n7=atoi(tok?tok:"");
        tok=strtok(NULL,","); n8=atoi(tok?tok:"");
      }
      else {
        cout << "Warning: abaqus element type " << eltype
             << " not converted to tochnog." << "\n";
        continue;
      }
      // filter element types by input_abaqus_name (if specified)
      convert_ok = 1;
      if ( input_abaqus_name ) {
        convert_ok = 0;
        for ( i=0; i<nn; i++ ) {
          if ( str2[0]=='-' && -input_abaqus_name[i]==db_number(&str2[1]) )
            convert_ok = 1;
        }
      }
      // filter element numbers by input_abaqus_set (if specified)
      if ( convert_ok && input_abaqus_set ) {
        in_set = 0;
        for ( j=0; j<nset_; j++ ) {
          if ( input_abaqus_set[j]==elem_id ) in_set = 1;
        }
        if ( !in_set ) convert_ok = 0;
      }
      if ( !convert_ok ) continue;
      out << "element  " << elem_id << "  " << str2;
      out << "  " << n1 << "  " << n2;
      if ( nnol>=3 ) out << "  " << n3;
      if ( nnol>=4 ) out << "  " << n4;
      if ( nnol>=6 ) out << "  " << n5 << "  " << n6;
      if ( nnol>=8 ) out << "  " << n7 << "  " << n8;
      if ( nnol>=9 ) out << "  " << n9;
      out << "\n";
    }
  }

  // write group_* records from the abaqus material (if requested)
  if ( input_abaqus_group==-YES && young>0. ) {
    out << "group_type 0  -materi\n";
    out << "group_materi_elasti_young 0  " << young << "\n";
    if ( poisson>0. )
      out << "group_materi_elasti_poisson 0  " << poisson << "\n";
    out << "group_materi_memory 0  -updated_without_rotation\n";
  }

  out << "end_data\n";
  in.close();
  out.close();
  cout << "Generated tochnog_abaqus.dat from abaqus.inp." << "\n";

}

void input_feflow_read( void )

{
  // read the feflow mesh file feflow.fem (ASCII) and fill node and element.
  // If input_feflow_fem is set to -no, the mesh is read from feflow.dac.
  long int inode=0, ielem=0, nnol=0, n1=0, n2=0, n3=0, n4=0, n5=0,
    n6=0, n7=0, n8=0, n9=0, idum[1], nn=0, i=0, input_feflow_fem=-YES,
    ldum=0;
  double ddum[1], xyz[MDIM];
  char line[MCHAR], *tok=NULL, str2[MCHAR], filename[MCHAR];
  ifstream in;
  long int in_nodes=0, in_elements=0;

  db( INPUT_FEFLOW_FEM, 0, &input_feflow_fem, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( input_feflow_fem==-NO )
    strcpy( filename, "feflow.dac" );
  else
    strcpy( filename, "feflow.fem" );
  in.open( filename );
  if ( !in ) {
    pri( "Error: cannot open feflow file ", filename );
    exit(TN_EXIT_STATUS);
  }

  in_nodes = in_elements = 0;
  while ( in.getline(line,MCHAR) ) {
    string_convert_to_lower_case( line );
    char *p = line;
    while ( *p==' ' || *p=='\t' ) p++;
    if ( *p=='\0' ) continue;

    // section keywords
    if ( strstr(p,"coordinates") || strstr(p,"nodes") ) {
      in_nodes = 1; in_elements = 0; continue;
    }
    if ( strstr(p,"elements") ) {
      in_elements = 1; in_nodes = 0; continue;
    }
    // skip other headers/comments
    if ( strstr(p,"problem") || strstr(p,"dimension") || strstr(p,"version")
         || strstr(p,"feflow") || p[0]=='#' || p[0]=='*' || p[0]=='/' ) continue;

    if ( in_elements ) {
      // element: first token is element number, then node numbers
      strcpy( str2, "" );
      tok = strtok(p," ,\t");
      if ( !tok ) continue;
      ielem = atoi(tok);
      nnol = 0; n1=n2=n3=n4=n5=n6=n7=n8=n9=0;
      tok = strtok(NULL," ,\t");
      while ( tok && nnol<9 ) {
        nnol++;
        if      ( nnol==1 ) n1=atoi(tok);
        else if ( nnol==2 ) n2=atoi(tok);
        else if ( nnol==3 ) n3=atoi(tok);
        else if ( nnol==4 ) n4=atoi(tok);
        else if ( nnol==5 ) n5=atoi(tok);
        else if ( nnol==6 ) n6=atoi(tok);
        else if ( nnol==7 ) n7=atoi(tok);
        else if ( nnol==8 ) n8=atoi(tok);
        else if ( nnol==9 ) n9=atoi(tok);
        tok = strtok(NULL," ,\t");
      }
      if      ( nnol==2 ) strcpy(str2,"-bar2");
      else if ( nnol==3 ) strcpy(str2,"-tria3");
      else if ( nnol==4 ) strcpy(str2,"-quad4");
      else if ( nnol==6 ) strcpy(str2,"-tria6");
      else if ( nnol==8 ) strcpy(str2,"-quad9");
      else continue;
      nn = 1+nnol;
      long int edata[1+9];
      edata[0] = -db_number(&str2[1]);
      edata[1]=n1; edata[2]=n2; edata[3]=n3; edata[4]=n4; edata[5]=n5;
      edata[6]=n6; edata[7]=n7; edata[8]=n8; edata[9]=n9;
      db( ELEMENT, ielem, edata, ddum, nn, VERSION_NORMAL, PUT );
    }
    else if ( in_nodes ) {
      // node: first token is node number, then coordinates
      tok = strtok(p," ,\t");
      if ( !tok ) continue;
      inode = atoi(tok);
      xyz[0]=xyz[1]=xyz[2]=0.;
      i=0;
      tok = strtok(NULL," ,\t");
      while ( tok && i<3 ) { xyz[i]=atof(tok); tok=strtok(NULL," ,\t"); i++; }
      db( NODE, inode, idum, xyz, ndim, VERSION_NORMAL, PUT );
    }
  }

  in.close();
  cout << "Input from feflow file " << filename << " read." << "\n";

}
