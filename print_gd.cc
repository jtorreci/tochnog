/*
    copyright (c) 1998  dennis roddeman
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

void print_giddata( long int task )

{
  long int imat=0, idat=0, ndat=0, igen=0, intv=0;
  char defaul[MCHAR];

  ofstream outprb( "tochnog.prb" );
  ofstream outind( "tochnog.ind" );
  ofstream outmat( "tochnog.mat" );
  ofstream outfirstbas( "tochnog.firstbas" );
  ofstream outcmd( "tochnog.cmd" );

    // bas file
  outfirstbas << "echo *GenData(1)\n";
  outfirstbas << "number_of_space_dimensions *ndime\n";

  task = labs(task);

    // problem data and bas file
  if      ( task==ALL )
    outprb << "21\n";
  else if ( task==GROUND )
    outprb << "6\n";
  else if ( task==MAXFRE )
    outprb << "4\n";
  else if ( task==MAXTIM )
    outprb << "4\n";
  else if ( task==STRESS )
    outprb << "18\n";
  else if ( task==THERMAL )
    outprb << "5\n";
  else if ( task==WAVE )
    outprb << "6\n";
  outprb << "echo#CB#(-yes,-no) -yes\n";
  igen = 1;
  outprb << "number_of_space_dimensions#CB#(1,2,3) 2\n";
  igen++;
  if ( task!=MAXFRE && task!=MAXTIM ) {
    outprb << "derivatives#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "derivatives" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( condif_temperature ) {
    if ( task==THERMAL )
      strcpy( defaul, "On" );
    else
      strcpy( defaul, "Off" );
    outprb << "condif_temperature#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "condif_temperature" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( groundflow_pressure ) {
    if ( task==GROUND )
      strcpy( defaul, "On" );
    else
      strcpy( defaul, "Off" );
    outprb << "groundflow_velocity#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "groundflow_velocity" << "\n";
    outfirstbas << "*endif\n";
    outprb << "groundflow_pressure#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "groundflow_pressure" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( materi_velocity ) {
    if ( task==STRESS )
      strcpy( defaul, "On" );
    else
      strcpy( defaul, "Off" );
    outprb << "materi_damage#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_damage" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_density#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_density" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_displacement#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_displacement" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_history_variables#CB#(0,1,2,3,4,5,6,7,8,9,10,11,12) 0\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "0" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_history_variables" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_maxwell_stress#CB#(0,1,2,3,4,5,6,7,8,9,10,11,12) 0\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "0" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_maxwell_stress" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_plasti_kappa#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_maxwell_stress" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_plasti_rho#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_plasti_rho" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_strain_elasti#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_strain_elasti" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_strain_plasti#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_strain_plasti" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_strain_total#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_strain_total" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_stress#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_stress" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_velocity#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_velocity" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_void_fraction#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_void_fraction" << "\n";
    outfirstbas << "*endif\n";
    outprb << "materi_work#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "materi_work" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( maxwell_er ) {
    if ( task==MAXFRE )
      strcpy( defaul, "On" );
    else
      strcpy( defaul, "Off" );
    outprb << "maxwell_er#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "maxwell_er" << "\n";
    outfirstbas << "*endif\n";
    outprb << "maxwell_ei#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "maxwell_ei" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( maxwell_e ) {
    if ( task==MAXTIM )
      strcpy( defaul, "On" );
    else
      strcpy( defaul, "Off" );
    outprb << "maxwell_e#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "maxwell_e" << "\n";
    outfirstbas << "*endif\n";
    outprb << "maxwell_fe#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "maxwell_fe" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( wave_scalar ) {
    if ( task==WAVE )
      strcpy( defaul, "On" );
    else
      strcpy( defaul, "Off" );
    outprb << "wave_scalar#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "wave_scalar" << "\n";
    outfirstbas << "*endif\n";
    outprb << "wave_fscalar#CB#(Off,On) ";
    outprb << defaul << "\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "wave_fscalar" << "\n";
    outfirstbas << "*endif\n";
  }
  if ( task!=MAXFRE && task!=MAXTIM ) {
    outprb << "residue#CB#(Off,On) Off\n";
    igen++;
    outfirstbas << "*if(strcmp(GenData(" << igen << ")," 
      << '"' << "Off" << '"' << ")!=0)" << "\n";
    outfirstbas << "residue" << "\n";
    outfirstbas << "*endif\n";
  }

  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_partialname( idat, "control" ) && db_external( idat ) &&
        check( idat, CHECK_USAGE ) ) {
      ndat++;
    }
  }
  outprb << "1\n";
  outprb << ndat << "\n";
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_partialname( idat, "control" ) && db_external( idat ) && 
        check( idat, CHECK_USAGE ) ) {
      if      ( idat==CONTROL_TIMESTEP )
        outprb << db_name(idat) << " 1,1\n";
      else if ( idat==CONTROL_TIMESTEP_ITERATIONS && task==MAXFRE )
        outprb << db_name(idat) << " 1\n";
      else if ( idat==CONTROL_TIMESTEP_ITERATIONS && task==MAXTIM )
        outprb << db_name(idat) << " 1\n";
      else if ( idat==CONTROL_TIMESTEP_ITERATIONS && task==THERMAL )
        outprb << db_name(idat) << " 1\n";
      else if ( idat==CONTROL_PRINT )
        outprb << db_name(idat) << " -node_dof\n";
      else
        outprb << db_name(idat) << " -\n";
    }
  }
  outfirstbas << "end_initia\n";

  outfirstbas << "*loop intervals\n";
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_partialname( idat, "control_" ) && db_external( idat ) &&
        check( idat, CHECK_USAGE ) ) {
      intv++;
      outfirstbas << "*if(strcmp(IntvData(" << intv << ")," 
        << '"' << "-" << '"' << ")!=0)" << "\n";
      outfirstbas << db_name(idat) << " *LoopVar" << " *IntvData(" << intv << ")" << "\n";
      outfirstbas << "*endif\n";
    }
  }
  outfirstbas << "*end intervals\n";

    // material file and bas file 

  outmat << "NUMBER: 1 MATERIAL: GROUP0\n";
  outfirstbas << "*loop materials\n";
  outfirstbas << "*Set var material=0\n";
  outfirstbas << "*Set var condif=0\n";
  outfirstbas << "*Set var groundflow=0\n";
  outfirstbas << "*Set var maxwell_frequency=0\n";
  outfirstbas << "*Set var maxwell_time=0\n";
  outfirstbas << "*Set var wave=0\n";
  outfirstbas << "*Set var group=LoopVar-1\n";
  if ( materi_velocity ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_partialname( idat, "group_materi" ) && db_external( idat ) ) {
        imat++;
        outmat << "QUESTION: " << db_name(idat) << "\n";
        outmat << "VALUE: -\n";
        outfirstbas << "*if(strcmp(MatProp(" << imat << ")," 
          << '"' << "-" << '"' << ")!=0)" << "\n";
        outfirstbas << db_name(idat) << " *group" << " *MatProp(" << imat << ")" << "\n";
        outfirstbas << "*Set var material=1\n";
        outfirstbas << "*endif\n";
      }
    }
  }
  if ( condif_temperature ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_partialname( idat, "group_condif" ) && db_external( idat ) ) {
        imat++;
        outmat << "QUESTION: " << db_name(idat) << "\n";
        outmat << "VALUE: -\n";
        outfirstbas << "*if(strcmp(MatProp(" << imat << ")," 
          << '"' << "-" << '"' << ")!=0)" << "\n";
        outfirstbas << db_name(idat) << " *group" << " *MatProp(" << imat << ")" << "\n";
        outfirstbas << "*Set var condif=1\n";
        outfirstbas << "*endif\n";
      }
    }
  }
  if ( groundflow_pressure ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_partialname( idat, "group_groundflow" ) && db_external( idat ) ) {
        imat++;
        outmat << "QUESTION: " << db_name(idat) << "\n";
        outmat << "VALUE: -\n";
        outfirstbas << "*if(strcmp(MatProp(" << imat << ")," 
          << '"' << "-" << '"' << ")!=0)" << "\n";
        outfirstbas << db_name(idat) << " *group" << " *MatProp(" << imat << ")" << "\n";
        outfirstbas << "*Set var groundflow=1\n";
        outfirstbas << "*endif\n";
      }
    }
  }
  if ( maxwell_er ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_partialname( idat, "group_maxwell_frequency" ) && db_external( idat ) ) {
        imat++;
        outmat << "QUESTION: " << db_name(idat) << "\n";
        outmat << "VALUE: -\n";
        outfirstbas << "*if(strcmp(MatProp(" << imat << ")," 
          << '"' << "-" << '"' << ")!=0)" << "\n";
        outfirstbas << db_name(idat) << " *group" << " *MatProp(" << imat << ")" << "\n";
        outfirstbas << "*Set var maxwell_frequency=1\n";
        outfirstbas << "*endif\n";
      }
    }
  }
  if ( maxwell_e ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_partialname( idat, "group_maxwell_time" ) && db_external( idat ) ) {
        imat++;
        outmat << "QUESTION: " << db_name(idat) << "\n";
        outmat << "VALUE: -\n";
        outfirstbas << "*if(strcmp(MatProp(" << imat << ")," 
          << '"' << "-" << '"' << ")!=0)" << "\n";
        outfirstbas << db_name(idat) << " *group" << " *MatProp(" << imat << ")" << "\n";
        outfirstbas << "*Set var maxwell_time=1\n";
        outfirstbas << "*endif\n";
      }
    }
  }
  if ( wave_scalar ) {
    for ( idat=0; idat<MDAT; idat++ ) {
      if ( db_partialname( idat, "group_wave" ) && db_external( idat ) ) {
        imat++;
        outmat << "QUESTION: " << db_name(idat) << "\n";
        outmat << "VALUE: -\n";
        outfirstbas << "*if(strcmp(MatProp(" << imat << ")," 
          << '"' << "-" << '"' << ")!=0)" << "\n";
        outfirstbas << db_name(idat) << " *group" << " *MatProp(" << imat << ")" << "\n";
        outfirstbas << "*Set var wave=1\n";
        outfirstbas << "*endif\n";
      }
    }
  }
  outmat << "QUESTION: group_integration_points\n";
  outmat << "VALUE: -\n";
  outmat << "END MATERIAL\n";
  outfirstbas << "*Set var group=LoopVar-1\n";
  outfirstbas << "group_type *group\n";
  outfirstbas << "*if(material==1&&groundflow==1&&condif==1)\n-materi -groundflow -condif\n";
  outfirstbas << "*elseif(material==1&&groundflow==1)\n-materi -groundflow\n";
  outfirstbas << "*elseif(material==1&&condif==1)\n-materi -condif\n";
  outfirstbas << "*elseif(groundflow==1&&condif==1)\n-groundflow -condif\n";
  outfirstbas << "*elseif(material==1)\n-materi\n";
  outfirstbas << "*elseif(groundflow==1)\n-groundflow\n";
  outfirstbas << "*elseif(condif==1)\n-condif\n";
  outfirstbas << "*elseif(maxwell_frequency==1)\n-maxwell_frequency\n";
  outfirstbas << "*elseif(maxwell_time==1)\n-maxwell_time\n";
  outfirstbas << "*elseif(wave==1)\n-wave\n";
  outfirstbas << "*else\n-none\n*endif\n";
  outfirstbas << "*end materials\n";

    // write data part commands
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_external( idat ) && !db_print_only(idat) &&
         check( idat, CHECK_USAGE ) &&
         ( db_type(idat)==INTEGER || db_type(idat)==DOUBLE_PRECISION ) ) {
      outcmd << db_name(idat) << "\n";
    }
  }

    // write data part commands index info
  for ( idat=0; idat<MDAT; idat++ ) {
    if ( db_external( idat ) && !db_print_only(idat) &&
         check( idat, CHECK_USAGE ) &&
         ( db_type(idat)==INTEGER || db_type(idat)==DOUBLE_PRECISION ) ) {
       if ( db_no_index(idat) )
          outind << "0\n";
       else
          outind << "1\n";
    }
  }

  outprb.close();
  outind.close();
  outfirstbas.close();
  outcmd.close();
  outmat.close();
}
