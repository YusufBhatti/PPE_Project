#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

/* Program to create lookup table for ECHAM6/HAM aerosol water uptake
   scheme (Petters & Kreidenweis 'kappa-koehler' theory option)
   
   Compilation:
   -----------
   Requires netcdf header file and library. So if for example the 
   netcdf directory is /opt/netcdf/4.3.0 :

   gcc -o make_lut_kappa_v2 -I/opt/netcdf/4.3.0/include 
          make_lut_kappa_v2.c -lm -lc -L/opt/netcdf/4.3.0/lib -lnetcdf

   Usage:
   -----
   make_lut_kappa_v2 <filename>
   where <filename> is the name of the generated netcdf lookup table file.
   <filename> is optional and defaults to lut_kappa_v2.nc. Any further arguments are ignored.
   
*/

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

// dimension names
#define RD_NAME "DryRadius"
#define T_NAME "Temperature"
#define RH_NAME "RelativeHumidity"
#define K_NAME "kappa"

double rt_solver(double, double, double, double, double, double);
double eqn(double, double, double, double);

int main(int argc, char *argv[])
{
  /* constants */
 
  const double argas = 8.314409;
  const double Mw = 0.018;
  const double rho_w = 1000;

  const double Tmin = 233;
  const double Tmax = 313;
  const double Tstep = 5;

  const double kappamin = 0.04;
  const double kappamax = 1.14;
  const double kappastep = 0.02;

  const double Rdmin = 5E-10;
  const double Rdmax = 5E-6;
  const double Rdstep = 97;

  const double rhmin = 0.30;
  const double rhmax = 0.99;
  const double rhstep_init = 0.03;
  const double rhstep_shrink = 0.96;

  // parameters for numerical solution
  const double gf1 = 0.99;          // always gives negative solution
  double gf2 = 2.5;  
  const double tol = 0.01;

  // default filename and version number
  char def_fn[] = "lut_kappa_v2.nc";
  const char version[] = "2.0";
  const char gf_units[] = "um um-1";

  // netCDF stuff

  int ncfp, dimRd, dimT, dimRH, dimK, varRd, varT, varRH, varK, gf_id;
  size_t lenRd, lenT, lenRH, lenK;
  char* nc_gatt;
  static size_t startix[] = {0,0,0,0};
  size_t cnt[4];
  int gf_dimids[4];
  size_t bufsize;
  int ibufsize;
  int i, j, k, l, dummy;
  int retval;
  
  time_t rawtime;
  char timestr[80];

  // variables
  double Rd;
  double Rdratio;
  double T;
  double rh;
  double rhstep;
  double kappa;
  double sigma;
  double fac_A;
  char *fn;
  double *gf;
  double *dataRd;
  double *dataT;
  double *dataRH;
  double *dataK;

  // executable

  // command line
  if (argc > 2) {
    printf("Usage: %s %s \n", argv[0], "filename (optional)");
    printf("extraneous arguments following %s %s \n", argv[1], "ignored");
  }

  if (argc>=2) 
    fn = argv[1];
  else
    fn = def_fn;
  
  // create netcdf file
  retval = nc_create(fn, NC_CLOBBER, &ncfp);
  if (retval != NC_NOERR) ERR(retval);

  // dimensions: calculate lengths
  /*  lenRd = (size_t) (1+Rdstep);
  lenT = (size_t) ceil((Tmax-Tmin)/Tstep);
  lenRH = (size_t) ceil((rhmax-rhmin)/rhstep);
  lenK = (size_t) ceil((kappamax-kappamin)/kappastep); */

  Rdratio = exp((log(Rdmax)-log(Rdmin))/Rdstep);
      
  // length of coordinate variables - needed for netcdf dimensions 
  // implemented like this for consistency with the data calculation loops below

  // Straightforward cases: for T and kappa, constant additive steps
  i = 0;
  for (T=Tmin; T<=Tmax; T+=Tstep) {
    i++;
  }
  lenT = (size_t) i;

  i = 0;
  for (kappa=kappamin; kappa<=kappamax; kappa+=kappastep) {
    i++;
  }
  lenK = (size_t) i;

  // Dry radius: spacing is linear in ln(Rd) so constant ratio between successive steps
  i = 0;
  for (Rd=Rdmin; Rd<=Rdmax; Rd*=Rdratio) {
    i++;
  }
  lenRd = (size_t) i;
 
  // RH: special case, non-linear spacing
  i = 0;
  rh = rhmin;
  rhstep=rhstep_init;
  while (rh<=rhmax) {
      i++;
      rh = rh + rhstep;
      rhstep = rhstep * rhstep_shrink;
  }
  lenRH = (size_t) i;

  // define netcdf dimensions
  retval = nc_def_dim(ncfp, RD_NAME, lenRd, &dimRd);
  if (retval != NC_NOERR) ERR(retval);
  
  retval = nc_def_dim(ncfp, T_NAME, lenT, &dimT);
  if (retval != NC_NOERR) ERR(retval);
  
  retval= nc_def_dim(ncfp, RH_NAME, lenRH, &dimRH);
  if (retval != NC_NOERR) ERR(retval);
  
  retval = nc_def_dim(ncfp, K_NAME, lenK, &dimK);
  if (retval != NC_NOERR) ERR(retval);

  // define netCDF variables to hold the coordinate values
  retval = nc_def_var(ncfp, RD_NAME, NC_DOUBLE, 1, &dimRd, &varRd);
  if (retval != NC_NOERR) ERR(retval);

  retval = nc_def_var(ncfp, T_NAME, NC_DOUBLE, 1, &dimT, &varT);
  if (retval != NC_NOERR) ERR(retval);

  retval = nc_def_var(ncfp, RH_NAME, NC_DOUBLE, 1, &dimRH, &varRH);
  if (retval != NC_NOERR) ERR(retval);

  retval = nc_def_var(ncfp, K_NAME, NC_DOUBLE, 1, &dimK, &varK);
  if (retval != NC_NOERR) ERR(retval);

  // now define the growth factor 4-D array
  gf_dimids[0] = dimRd;
  gf_dimids[1] = dimT;
  gf_dimids[2] = dimRH;
  gf_dimids[3] = dimK;
  
  // define netcdf variables
  retval = nc_def_var(ncfp, "GF", NC_DOUBLE, 4, gf_dimids, &gf_id);
  if (retval != NC_NOERR) ERR(retval);

  // define global attributes
  time(&rawtime);
  strftime(timestr, 80, "%c", localtime(&rawtime));

  retval = nc_put_att_text(ncfp, NC_GLOBAL, "Created_date", strlen(timestr), timestr);
  if (retval != NC_NOERR) ERR(retval);   

  nc_gatt = getenv("USER");
  retval = nc_put_att_text(ncfp, NC_GLOBAL, "Created_by", strlen(nc_gatt), nc_gatt);
  if (retval != NC_NOERR) ERR(retval);
      
  retval = nc_put_att_text(ncfp, NC_GLOBAL, "Created_using", strlen(argv[0]), argv[0]);
  if (retval != NC_NOERR) ERR(retval);
      
  retval = nc_put_att_text(ncfp, NC_GLOBAL, "Version", strlen(version), version);
  if (retval != NC_NOERR) ERR(retval);

  // define coordinate attributes
  retval = nc_put_att_text(ncfp, varRd, "units", strlen("m"), "m");
  if (retval != NC_NOERR) ERR(retval);

  retval = nc_put_att_text(ncfp, varT, "units", strlen("Kelvin"), "Kelvin");
  if (retval != NC_NOERR) ERR(retval);

  retval = nc_put_att_text(ncfp, varRH, "units", strlen("fractional"), "fractional");
  if (retval != NC_NOERR) ERR(retval);

  retval = nc_put_att_text(ncfp, varK, "units", strlen("dimensionless"), "dimensionless");
  if (retval != NC_NOERR) ERR(retval);

  // define variable attributes
  retval = nc_put_att_text(ncfp, gf_id,  "units", strlen(gf_units), gf_units);
  if (retval != NC_NOERR) ERR(retval);

  // end netcdf definition mode
  retval = nc_enddef(ncfp);
  if (retval != NC_NOERR) ERR(retval);

  // memory allocation
  dataRd = (double*) malloc(lenRd*sizeof(double));
  if (dataRd==NULL) {printf("Allocation error for Rd"); exit(1);}
  dataT = (double*) malloc(lenT*sizeof(double));
  if (dataT==NULL) {printf("Allocation error for T"); exit(1);}
  dataRH = (double*) malloc(lenRH*sizeof(double));
  if (dataRH==NULL) {printf("Allocation error for RH"); exit(1);}
  dataK = (double*) malloc(lenK*sizeof(double));
  if (dataK==NULL) {printf("Allocation error for K"); exit(1);}

  bufsize = lenRd*lenT*lenRH*lenK*sizeof(double);
  gf = (double*) malloc(bufsize);
  if (gf==NULL) {printf("Allocation error for gf"); exit(1);}
  
  cnt[0] = lenRd;
  cnt[1] = lenT;
  cnt[2] = lenRH;
  cnt[3] = lenK;
  ibufsize = (int) bufsize;

  // data for coordinate variables
  // Straightforward cases: for T and kappa, constant additive steps
  i = 0;
  for (T=Tmin; T<=Tmax; T+=Tstep) {
    dataT[i] = T;
    i++;
  }
  retval = nc_put_var_double(ncfp, varT, dataT);
  if (retval != NC_NOERR) {free(dataT); ERR(retval);}
    
  i = 0;
  for (kappa=kappamin; kappa<=kappamax; kappa+=kappastep) {
    dataK[i] = kappa;
    i++;
  }
  retval = nc_put_var_double(ncfp, varK, dataK);
  if (retval != NC_NOERR) {free(dataK); ERR(retval);}

  // Dry radius: spacing is linear in ln(Rd) so constant ratio between successive steps
  i = 0;
  for (Rd=Rdmin; Rd<=Rdmax; Rd*=Rdratio) {
    dataRd[i] = Rd;
    i++;
  }
  retval = nc_put_var_double(ncfp, varRd, dataRd);
  if (retval != NC_NOERR) {free(dataRd); ERR(retval);}

  // RH: special case, non-linear spacing
  i=0;
  rh = rhmin;
  rhstep = rhstep_init;
  while (rh<=rhmax) { 
    dataRH[i] = rh;
    rh = rh + rhstep;
    rhstep = rhstep * rhstep_shrink;
    i++;
  }
  retval = nc_put_var_double(ncfp, varRH, dataRH);
  if (retval != NC_NOERR) {free(dataRH);ERR(retval);}

  i = 0;  
  for (Rd=Rdmin; Rd<=Rdmax; Rd*=Rdratio)           // NBB MULTIPLICATIVE step for logarithmic scale
  {
      for (T=Tmin; T<=Tmax; T+=Tstep)
      {
	  // temperature-dependent surface tension after Seinfeld & Pandis
	  sigma = 0.0761-1.55e-4*(T-273);

	  fac_A = 2*sigma*Mw/(argas*T*Rd*rho_w);
	  
	  rh = rhmin;
	  rhstep = rhstep_init;
	  while (rh<=rhmax) 
	  {
	      rh = rh + rhstep;
	      rhstep = rhstep * rhstep_shrink;

	      for (kappa=kappamin; kappa<=kappamax; kappa+=kappastep) 
	      {
         
		  gf[i] = rt_solver(gf1, gf2, tol, fac_A, rh, kappa);

                  if (gf[i] < 1) gf[i] = 1;
		  gf2 = gf[i]+0.5;

		  i++;

	      }
	  }
      }
  }

  retval = nc_put_vara_double(ncfp, gf_id, startix, cnt, gf);
  if (retval != NC_NOERR) {free(gf); ERR(retval);}

  dummy=nc_close(ncfp);

  free(dataRd);
  free(dataT);
  free(dataRH);
  free(dataK);
  free(gf);

  return 0;
}

/* Bernt's method solver for a scalar-valued function eqn(x)
   of a scalar x */

double rt_solver(double xlow, double xhi, double xtol, double d1, 
		 double d2, double d3)
{
  double za,zb,zc,zd,ze,zeps,fa,fb,fc,ztol1,zxm,zp,zq,zr,zs;

  // executable

  zeps = 4*DBL_EPSILON;
  ztol1 = zeps+1;

  za=xlow;
  zb=xhi;
  fa=eqn(za,d1,d2,d3);
  fb=eqn(zb,d1,d2,d3);

  //     check that f(ax) and f(bx) have different signs

  if (!((fa == 0)|| (fb == 0))) {
    if (fa * (fb/fabs(fb)) > 0) {
      printf("zeroin: f(ax) and f(bx) have the same signs, zeroin returning 0\n");
      printf("zeroin: input was %e %e %e \n", d1, d2, d3);
      return 0;
    }
  }

  zc=za;
  fc=fa;
  zd=zb-za;
  ze=zd;

  if (fabs(fc)< fabs(fb)) {
    za=zb;
    zb=zc;
    zc=za;
    fa=fb;
    fb=fc;
    fc=fa;
  }

  ztol1=2*zeps*fabs(zb)+0.5*xtol;
  zxm = 0.5*(zc-zb);
     
  while ((fabs(zxm) > ztol1) && (fb != 0)) {

    // see IF a bisection is forced

    if ((fabs(ze)>=ztol1) && (fabs(fa)>fabs(fb))) {

      zs=fb/fa;

      if (za != zc) {
             // inverse quadratic interpolation

	zq = fa / fc;
	zr = fb / fc;
	zp=zs*(2*zxm*zq*(zq-zr)-(zb-za)*(zr-1));
	zq=(zq-1)*(zr-1)*(zs-1);
      }

      else {
               // linear interpolation
 
	zp=2*zxm*zs;
	zq=1.0-zs;
      }


      if (zp < 0) {
	zp = -zp ;
      }
      else {
	zq = -zq;
      }

      zs=ze;
      ze=zd;
      if (((2*zp)>(3*zxm*zq-fabs(ztol1*zq)))||
           (zp>fabs(0.5*zs*zq))) {
	zd = zxm ;
	ze = zd;
      }
      else {
	zd=zp/zq;
      }
    }
    else {
      zd=zxm;
      ze=zd;
    }

    za=zb;
    fa=fb;

    if (fabs(zd)<=ztol1) {
      if (zxm < 0) {
	zb = zb-ztol1;
      }
      else {
	zb = zb+ztol1;
      }
    }
    else {
	zb = zb + zd;
    }


    fb=eqn(zb,d1,d2,d3);

    if ((fb*(fc/fabs(fc)))>0) {
      zc=za;
      fc=fa;
      zd=zb-za;
      ze=zd;
    }

    if (fabs(fc)< fabs(fb)) {
      za=zb;
      zb=zc;
      zc=za;
      fa=fb;
      fb=fc;
      fc=fa;
    }

    ztol1=2*zeps*fabs(zb)+0.5*xtol;
    zxm = 0.5*(zc-zb);

  }

  return zb;

}

double eqn(double gf_est, double A, double rh, double k)
{
  double gf3_1;

  // executable

  gf3_1 = pow(gf_est,3) - 1;
  
  return (rh*(exp(-A/gf_est)) - 
	  gf3_1 / (gf3_1 + k)  );

}

