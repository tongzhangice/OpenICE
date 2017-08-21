#ifndef ICE_NETCDF_

/*
 * Netcdf
 * Variable indecies.
 *
 *  */


#if ICE_BENCH_TEST
#elif ESIMINT_TEST || HEINO_TEST 
#elif TEST_CASE == ICE_GREEN_LAND
enum {
    data_index_dhdt		= 1,
    data_index_surfvelmag	= 2,
    data_index_lat		= 3,
    data_index_balvelmag	= 4,
    data_index_bheatflx		= 5,
    data_index_presprcp		= 6,
    data_index_lon		= 7,
    data_index_usrf		= 8,
    data_index_surfvely		= 9,
    data_index_surfvelx		= 10,
    data_index_topg		= 11,
    data_index_thk		= 12,
    data_index_runoff		= 13,
    data_index_smb		= 14,
    data_index_airtemp2m	= 15,
    data_index_surftemp		= 16,
    data_index_temp		= 17,
    data_index_x1		= 18,
    data_index_y1       	= 19,
};

void read_nc_data(char *file_name);
void nc_data_set_active(int var_index);
void interp_nc_data_3D(double x, double y, double z, int level, double *a);
#define interp_nc_data(x, y, z, a) \
    interp_nc_data_3D(x, y, z, 0, a)
extern double nc_data_scaling;
void func_nc_data(double x, double y, double z, double *a);

#endif

#define ICE_NETCDF_
#endif 
