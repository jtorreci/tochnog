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

void  user_bounda_time( long int index, double time_current, double &val )
 
  /* Purpose: user supplied routine for time table for
       bounda_unknown index.

     You need to skip the error message.

     You need to specify val for all times.
 
     index: Input. Index of the bounda_unknown record for which the time
       value is to be given.
     time_current: Input. Current time.
     val: Output. Value at time (typically specified here with a function).
  */
 
{
  cout << "\nError: routine user_bounda_time not programmed.\n";
  exit(0);
}               

void  user_change_dataitem_time( long int index, double time_current, double &val )

  /* Purpose: user supplied routine for time table for
       change_dataitem index.

     You need to skip the error message.
 
     index: Input. Index of the change_dataitem record for which the time
       value is to be given.
     time_current: Input. Current time.
     val: Output. Value at time (typically specified here with a function).
  */

{
  cout << "\nError: routine user_change_dataitem_time not programmed.\n";
  exit(0);
}

void  user_change_geometry_time( long int index, double time_current, double &val )

  /* Purpose: user supplied routine for time table for
       change_geometry index.

     You need to skip the error message.
 
     index: Input. Index of the change_geometry record for which the time
       value is to be given.
     time_current: Input. Current time.
     val: Output. Value at time (typically specified here with a function).
  */

{
  cout << "\nError: routine user_change_geometry_time not programmed.\n";
  exit(0);
}

void user_plasti( long int task, double user_data[], double new_unknowns[], 
  double old_hisv[], double new_hisv[], double new_stress[], double &f )

  /* Purpose: user supplied routine for calculation of plasticity yield 
       function and flow rule function.

     You need to skip the error message.

     task: Input. 1: Yield function. 2: Flow rule function.
     user_data[]: Input. User supplied data values (if specified
       in datafile).
     new_unknowns: Input. It contains the primary dof values at this
       integration point at the new time point t+dt.
     old_hisv: Input. It contains the old estimates for the 
       materi_history_variables (if initialized in the data file).
     new_hisv: Output. On output it should contain
       the new estimates for the materi history variables (if initialized in the data file).
     new_stress: Input. It contains the most recent estimate for
       the new stresses (if initialized in the data file). 
     f: Output. Yield function if task==1 and flow rule function if task==2.


     Example: f_yield = temp * stresszz * stresszz - user_data[0].
              f_flow = user_data[1] * f_yield 

     Assume that materi_stress is initialized first,
     and condif_temperature is initialized second.
     Then new_unknowns[] contains:
     stressxx, stressxy, stressxz, stressyy, stressyz, stresszz, temp

    double f_field;

    f_yield = new_unknowns[6] * new_stress[5] * new_stress[5] - user_data[0];
    if ( task==1 )
      f = f_yield;
    else {
      assert( task==2 );
      f = user_data[1] * f_yield;
  */

{
  cout << "\nError: routine user_plasti not programmed.\n";
  exit(0);
}


void user_sigma( double user_data[], double new_unknowns[], 
  double inc_epe[], double old_hisv[],
  double new_hisv[], double old_stress[], double new_stress[], 
  double Cuser[3][3][3][3] )

  /* Purpose: user supplied routine for adding incremental stresses 
     (stress increments over dt).

     Note: this routine should only be used when materi_strain_total is
     initialised.

     user_data: Input. User supplied data values (if specified
       in datafile).
     new_unknowns: Input. It contains the primary dof values at this
       integration point at the new time point.
     inc_epe: Input. It contains the most recent estimate for
       the incremental elastic strains, i.e.
       the total strain increment minus the plastic strain increment (inc_epexx,
       inc_epexy, inc_epexz, inc_epeyx, inc_epeyy, inc_epeyz, inc_epezx,
       inc_epezy, inc_epezz).
     new_epe: Input. It contains the most recent estimate for
       the elastic strains (if initialized in the data file).
     old_hisv: Input. On input it contains
       the old estimates for the materi history variables (if initialized in the data file).
     new_hisv: Output. On output it should contain
       the new estimates for the materi history variables (if initialized in the data file).
     old_stress: Input. The old elastic stresses.
       For all number_of_space_dimensions these iold 
       stresses are old_stressxx, old_stressxy, old_stressxz, old_stressyx, old_stressyy, 
       old_stressyz, old_ stresszx, old_stresszy, old_stresszz.
     new_stress: Output. The incremental elastic stresses should 
       be added. For all number_of_space_dimensions these 
       are inc_stressxx, inc_stressxy, inc_stressxz, inc_stressyx, inc_stressyy, 
       inc_stressyz, inc_stresszx, inc_stresszy, inc_stresszz.
     Cuser: Input and Output. The material stiffness matrix needs
       to be added. Cuser[i][j][k][l] = delta stress[i][j] / delta epsilon[[k][l].
       In case of membrane stresses you need to use the membrane stiffness.

  */                             
  
{
}

void user_viscosity( double user_data[], double new_unknowns[], double &visc )

  /* Purpose: user supplied routine for calculation of viscosity.

     user_data: Input. User supplied data values (if specified
       in datafile).
     new_unknowns: Input. It contains the primary dof values at this
       integration point at the new time point.
     visc: Output. Viscosity.

     Example: visc = user_data[0] * ( user_data[1] - temp );

     Assume that materi_strain_total is initialized first,
     and condif_temperature is initialized second.
     Then new_unknowns[] contains:
     eptxx, eptxy, eptxz, eptyy, eptyz, eptzz, temp

     visc = user_data[0] * ( user_data[1] - new_unknowns[6] );
  */

{
  cout << "\nError: routine user_viscosity not programmed.\n";
  exit(0);
}
