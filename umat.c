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
  /* Dummy UMAT template: all Abaqus arguments are intentionally unused. */
  (void)stress; (void)statev; (void)ddsdde; (void)sse; (void)spd; (void)scd;
  (void)rpl; (void)ddsddt; (void)drplde; (void)drpldt; (void)stran;
  (void)dstran; (void)time; (void)dtime; (void)temp; (void)dtemp;
  (void)predef; (void)dpred; (void)cmname; (void)ndi; (void)nshr;
  (void)ntens; (void)nstatv; (void)props; (void)nprops; (void)coords;
  (void)drot; (void)pnewdt; (void)celent; (void)dfgrd0; (void)dfgrd1;
  (void)noel; (void)npt; (void)layer; (void)kspt; (void)kstep; (void)kinc;
  (void)cmname_len;

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
