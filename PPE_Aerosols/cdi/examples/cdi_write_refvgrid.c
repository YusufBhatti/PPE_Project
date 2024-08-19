#include <stdio.h>
#include "cdi.h"

#define  nlon   12 // Number of longitudes
#define  nlat    6 // Number of latitudes 
#define  nlev    5 // Number of levels    
#define  nts     3 // Number of time steps

int main(void)
{
  int gridID, zaxisID1, zaxisID2, taxisID;
  int vlistID, varID1, varID2, streamID, tsID;
  int i, nmiss = 0;
  double levs[nlev] = {1, 2, 3, 4, 5};
  double var1[nlon*nlat];
  double var2[nlon*nlat*nlev];


  // Create a grid reference
  gridID = gridCreate(GRID_UNSTRUCTURED, nlon*nlat);
  gridDefNumber(gridID, 123);
  gridDefPosition(gridID, 3);
  gridDefReference(gridID, "http://www.x.y/gridfile.nc");
  gridDefUUID(gridID, "1234569887654321");

  // Create a surface level Z-axis
  zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a general vertical height Z-axis
  zaxisID2 = zaxisCreate(ZAXIS_REFERENCE, nlev);
  zaxisDefLevels(zaxisID2, levs);
  zaxisDefNlevRef(zaxisID2, nlev);
  zaxisDefNumber(zaxisID2, 71);
  //zaxisDefReference(zaxisID2, "http://www.x.y/vgridfile.nc");
  zaxisDefUUID(zaxisID2, "8765432112345678");
 
  // Create a variable list
  vlistID = vlistCreate();

  // Define the variables
  varID1 = vlistDefVar(vlistID, gridID, zaxisID1, TSTEP_INSTANT);
  varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TSTEP_INSTANT);

  // Define the variable names
  vlistDefVarName(vlistID, varID1, "varname1");
  vlistDefVarName(vlistID, varID2, "varname2");

  // Create a Time axis
  taxisID = taxisCreate(TAXIS_ABSOLUTE);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset in netCDF format
  streamID = streamOpenWrite("example.grb2", FILETYPE_GRB2);
  //streamID = streamOpenWrite("example.nc", FILETYPE_NC);
  if ( streamID < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return(1);
    }

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  // Loop over the number of time steps
  for ( tsID = 0; tsID < nts; tsID++ )
    {
      // Set the verification date to 1985-01-01 + tsID
      taxisDefVdate(taxisID, 19850101+tsID);
      // Set the verification time to 12:00:00
      taxisDefVtime(taxisID, 120000);
      // Define the time step
      streamDefTimestep(streamID, tsID);

      // Init var1 and var2
      for ( i = 0; i < nlon*nlat;      i++ ) var1[i] = 1.1;
      for ( i = 0; i < nlon*nlat*nlev; i++ ) var2[i] = 2.2;

      // Write var1 and var2
      streamWriteVar(streamID, varID1, var1, nmiss);
      streamWriteVar(streamID, varID2, var2, nmiss);
    }

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID1);
  zaxisDestroy(zaxisID2);
  gridDestroy(gridID);

  return 0;
}
