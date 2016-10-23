typedef long int integer;
typedef double doublereal;       
typedef short ftnlen;

int umat_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt,
	 drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, predef, 
	dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot, 
	pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc, 
	cmname_len)
doublereal *stress, *statev, *ddsdde, *sse, *spd, *scd, *rpl, *ddsddt, *
	drplde, *drpldt, *stran, *dstran, *time, *dtime, *temp, *dtemp, *
	predef, *dpred;
char *cmname;
integer *ndi, *nshr, *ntens, *nstatv;
doublereal *props;
integer *nprops;
doublereal *coords, *drot, *pnewdt, *celent, *dfgrd0, *dfgrd1;
integer *noel, *npt, *layer, *kspt, *kstep, *kinc;
ftnlen cmname_len;
{

  /* Dummy version of umat routine.

     You can use a umat.f routine as follows:

     - Install f2c on your computer (see http://www.netlib.org/f2c/)
     - f2c umat.f
     - Overwrite this dummy version with your umat.c
     - Outcomment (activate) the F2C statement in the makefile
     - make (compile and link)

     The umat.c will be called for 3d stress states (ntens=6).

     Tochnog MATERI_STRAIN_TOTAL will be mapped to Abaqus STRAN.

     Tochnog MATERI_HISTORY_VARIABLES will be mapped to Abaqus STATEV.

     Tochnog USER_DATA will be mapped to Abaqus PROPS.

     Tochnog iteration count will be mapped to Abaqus KINC.

  */

 return 0;
}
