#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <netcdf.h>
//#include "ins.h"
#include "greenland.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); abort();}

#define NX 301
#define NY 561

typedef struct NC_DATA_ {
    float x1[NX];
    float y1[NY];
    float usrf[NY][NX];
    float dhdt[NY][NX];
    float surfvelmag[NY][NX];
    float lat[NY][NX];
    float balvelmag[NY][NX];
    float bheatflx[NY][NX];
    float presprcp[NY][NX];
    float lon[NY][NX];
    float surfvely[NY][NX];
    float surfvelx[NY][NX];
    float topg[NY][NX];
    float thk[NY][NX];
    float runoff[NY][NX];
    float smb[NY][NX];
    float airtemp2m[NY][NX];
    float surftemp[NY][NX];
    float *list[100];
} NC_DATA;

static NC_DATA data;
static float *active_data;
static float *x1, *y1, dx, dy;

#define READ_DATA_1D(var) {					\
	if ((retval = nc_inq_varid(ncid, #var, &varid)))	\
	    ERR(retval);					\
	if ((retval = nc_get_var_float(ncid, varid,		\
				       &data.var[0]))) {	\
	    ERR(retval);					\
	} else	{						\
	    printf("   Read %-15s done.\n", #var);		\
	    data.list[data_index_##var] = &data.var[0];		\
	}							\
    }

#define READ_DATA_2D(var) {					\
	if ((retval = nc_inq_varid(ncid, #var, &varid)))	\
	    ERR(retval);					\
	if ((retval = nc_get_var_float(ncid, varid,		\
				       &data.var[0][0]))) {	\
	    ERR(retval);					\
	} else	{						\
	    printf("   Read %-15s done.\n", #var);		\
	    data.list[data_index_##var] = &data.var[0][0];	\
	}							\
    }



void
read_nc_data(char *file_name)
{
    int retval;
    int ncid, varid;

    /* Open the file. */
    if ((retval = nc_open(file_name, NC_NOWRITE, &ncid)))
	ERR(retval);

    READ_DATA_1D(x1); x1 = data.x1;
    READ_DATA_1D(y1); y1 = data.y1;
    dx = (x1[1] - x1[0]);
    dy = (y1[1] - y1[0]);

    READ_DATA_2D(dhdt);
    READ_DATA_2D(surfvelmag);
    READ_DATA_2D(lat);
    READ_DATA_2D(balvelmag);
    READ_DATA_2D(bheatflx);
    READ_DATA_2D(presprcp);
    READ_DATA_2D(lon);
    READ_DATA_2D(usrf);
    READ_DATA_2D(surfvely);
    READ_DATA_2D(surfvelx);
    READ_DATA_2D(topg);
    READ_DATA_2D(thk);
    READ_DATA_2D(runoff);
    READ_DATA_2D(smb);
    READ_DATA_2D(airtemp2m);
    READ_DATA_2D(surftemp);

    if ((retval = nc_close(ncid)))
	ERR(retval);
     
    printf("*** SUCCESS reading netcdf data file: %s.\n", 
	   file_name);
    return;
}

void 
nc_data_set_active(int var_index) 
{
    assert(var_index <= 18 && var_index >= 1);
    active_data = data.list[var_index];
    return;
}


/*
 * Bilinear interpolation of netCDF data.
 * Note: input coord is of computaional region,
 *       output data is scaled.
 * */
double nc_data_scaling = 1.;
double nc_length_scale = 1.;
void
interp_nc_data(double x, double y, double z, double *a)
{
    int i, j;
    double a00, a01, a10, a11;
    double wx, wy;

    x *= nc_length_scale;
    y *= nc_length_scale;
    z *= nc_length_scale;

#if 0
#else
    /* relative position */
    i = (int) x / dx ;
    j = (int) y / dy;

    /* only interior region, no interp near margin */
    assert(i < NX-1);
    assert(j < NY-1);
#endif

    /* Bilinear interpolation */
    a00 = active_data[j*NX + i];
    a01 = active_data[(j+1)*NX + i];
    a10 = active_data[j*NX + (i+1)];
    a11 = active_data[(j+1)*NX + (i+1)];

#if 0
    dx = x[i+1] - x[i];
    dy = y[j+1] - y[j];
#endif


#if 0
    /* check */
    a00 = (x) + 2*(y);
    a10 = (x+dx) + 2*(y);
    a01 = (x) + 2*(y+dy);
    a11 = (x+dx) + 2*(y+dy);
#endif

    wx = (x - (x1[i] - x1[0])) / dx;
    wy = (y - (y1[j] - y1[0])) / dy;
    a[0] = a00 * (1-wx) * (1-wy) 
	+ a01 * (1-wx) * wy
	+ a10 * wx * (1-wy)
	+ a11 * wx * wy;

    a[0] *= nc_data_scaling;
    return;
}
