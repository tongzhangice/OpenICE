#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ins.h"
//static int NX=1556, NY=1357, NODATA_VALUE=-9999, dx=1000, dy=1000;
static int NX=605, NY=271, NODATA_VALUE=-9999, dx=1000, dy=1000;
static double xllcorner=1662*1000, yllcorner=590*1000;

static int len = 8001;
static double x_end = 800e3/*, dx = 100.0*/;
double** read_txt_data(char *file_name)
{
    FILE *fp;
	int i, j, k=0;
	double **data;
    fp = fopen(file_name, "r");
    if (fp==NULL){
        printf("errors when opening 2D geo file!\n");
    }

	data = (double**) calloc(NY, sizeof *data);
	for (j=0;j<NY;j++)
	{
		data[j] = (double*) calloc(NX, sizeof *data[j]);
		for (i=0;i<NX;i++)
		{
            k++;
			fscanf(fp, "%lf", &data[j][i]);
		}
	}
	//free(data[0]);
    //free(data);
	return data;
}

#if 1
void interp_txt_data(double **data, FLOAT x, FLOAT y, FLOAT z, FLOAT *a)
{
	int i, j, ip1, jp1, i1;
	double a00, a01, a10, a11;
	//double wx, wy, val;
    //int NX=1518, NY=1324;
    //double xllcorner=792.5, yllcorner=-125.5;
    int interp_points=4;
    double x1[NX];
    double y1[NY];
    double eps=1e-3;
    
    for (i1=0;i1<NX;i1++)
		x1[i1] = i1*dx + xllcorner;
	for (i1=0;i1<NY;i1++)
		y1[i1] = i1*dy + yllcorner;
	
	i = (int) (x - xllcorner)/dx;
	j = (int) (y - yllcorner)/dy;

    if (i<NX-1&&j<NY-1)
    {
		
        a00 = data[NY-j-2][i];
	    a01 = data[NY-j-1][i];	    
	    a10 = data[NY-j-2][i+1];	    
	    a11 = data[NY-j-1][i+1];
		
		//printf("A00 A01 A10 A11: %lf %lf %lf %lf %d %lf %lf %d %d\n", a00, a01, a10, a11, interp_points, x, y, i, j);
        if (abs(a00-NODATA_VALUE)<eps){		
            a00 = 0;
            interp_points--;
            }
        if (abs(a01-NODATA_VALUE)<eps){		
            a01 = 0;
            interp_points--;
            }
        if (abs(a10-NODATA_VALUE)<eps){		
            a10 = 0;
            interp_points--;
            }
        if (abs(a11-NODATA_VALUE)<eps){
            a11 = 0;
            interp_points--;
            }
        
		//printf("a00 a01 a10 a11: %lf %lf %lf %lf %d %lf %lf\n", a00, a01, a10, a11, interp_points, x, y);
        //if (interp_points==0){
		//	printf("something terrible happened!\n");
		//}	

        if (interp_points>0&&interp_points<4){
            a[0] = (a00+a01+a10+a11)/interp_points;
            }
        
        else if (interp_points==0){
            if (abs(data[NY-j-2][i]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-2][i];
            else if (abs(data[NY-j][i]-NODATA_VALUE)>eps)
                a[0] = data[NY-j][i];
            else if (abs(data[NY-j-1][i-1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-1][i-1];
            else if (abs(data[NY-j-1][i+1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-1][i+1];
            else if (abs(data[NY-j-2][i-1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-2][i-1];
            else if (abs(data[NY-j-2][i+1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-2][i+1];
            else if (abs(data[NY-j][i-1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j][i-1];
            else if (abs(data[NY-j][i+1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j][i+1];
        }
        
        else{
            a[0]=1.0/((x1[i+1]-x1[i])*(y1[j]-y1[j+1]))*(a00*(x1[i+1]-x)*(y1[j]-y)+\
                    a10*(x-x1[i])*(y1[j]-y)+a01*(x1[i+1]-x)*(y-y1[j+1])+a11*(x-x1[i])*(y-y1[j+1]));
            }
    }
    else
    {
        if (j==NY){
            a[0] = data[NY-j][i];}
        else if (i==NX){
            a[0] = data[NY-j-1][i-1];}
        else 
            printf("%d\n", NY-j-1);
            a[0] = data[NY-j-1][i];
    }
	
}
#endif

#if 0
void interp_txt_data(double **data, FLOAT x, FLOAT y, FLOAT z, FLOAT *a)
{
	int i, j, ip1, jp1, i1;
	double a00, a01, a10, a11;
	//double wx, wy, val;
    //int NX=1518, NY=1324;
    //double xllcorner=792.5, yllcorner=-125.5;
    int interp_points=4;
    double x1[NX];
    double y1[NY];
    double eps=1e-3;

    if (y<0)
        y = -y;
    
    for (i1=0;i1<NX;i1++)
		x1[i1] = i1*dx + xllcorner;
	for (i1=0;i1<NY;i1++)
		y1[i1] = i1*dy + yllcorner;
	
	j = (int) (x - xllcorner)/dx;
	i = (int) (y - yllcorner)/dy;

    printf("i, j %d %d %f %f\n", i, j, x, y);

    if (i<NX-1&&j<NY-1)
    {
		
        a00 = data[NY-j-2][i];
	    a01 = data[NY-j-1][i];	    
	    a10 = data[NY-j-2][i+1];	    
	    a11 = data[NY-j-1][i+1];
		
        a00 = data[j][i];
	    a01 = data[NY-j-1][i];	    
	    a10 = data[NY-j-2][i+1];	    
	    a11 = data[NY-j-1][i+1];

		printf("A00 A01 A10 A11: %lf %lf %lf %lf %d %lf %lf %d %d\n", a00, a01, a10, a11, interp_points, x, y, i, j);
        if (abs(a00-NODATA_VALUE)<eps){		
            a00 = 0;
            interp_points--;
            }
        if (abs(a01-NODATA_VALUE)<eps){		
            a01 = 0;
            interp_points--;
            }
        if (abs(a10-NODATA_VALUE)<eps){		
            a10 = 0;
            interp_points--;
            }
        if (abs(a11-NODATA_VALUE)<eps){
            a11 = 0;
            interp_points--;
            }
        
		//printf("a00 a01 a10 a11: %lf %lf %lf %lf %d %lf %lf\n", a00, a01, a10, a11, interp_points, x, y);
        //if (interp_points==0){
		//	printf("something terrible happened!\n");
		//}	

        if (interp_points>0&&interp_points<4){
            a[0] = (a00+a01+a10+a11)/interp_points;
            }
        
        else if (interp_points==0){
            if (abs(data[NY-j-2][i]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-2][i];
            else if (abs(data[NY-j][i]-NODATA_VALUE)>eps)
                a[0] = data[NY-j][i];
            else if (abs(data[NY-j-1][i-1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-1][i-1];
            else if (abs(data[NY-j-1][i+1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-1][i+1];
            else if (abs(data[NY-j-2][i-1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-2][i-1];
            else if (abs(data[NY-j-2][i+1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j-2][i+1];
            else if (abs(data[NY-j][i-1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j][i-1];
            else if (abs(data[NY-j][i+1]-NODATA_VALUE)>eps)
                a[0] = data[NY-j][i+1];
        }
        
        else{
            a[0]=1.0/((x1[i+1]-x1[i])*(y1[j]-y1[j+1]))*(a00*(x1[i+1]-x)*(y1[j]-y)+\
                    a10*(x-x1[i])*(y1[j]-y)+a01*(x1[i+1]-x)*(y-y1[j+1])+a11*(x-x1[i])*(y-y1[j+1]));
            a[0] = (a00+a01+a10+a11)/interp_points;

            //printf("%f\n", a[0]);
            }
    }
    else
    {
        if (j==NY){
            a[0] = data[NY-j][i];}
        else if (i==NX){
            a[0] = data[NY-j-1][i-1];}
        else{
            printf("number %d %d %d\n", NY, j, NY-j-1);
            a[0] = data[NY-j-1][i];}
    }
	
}	
#endif

double* read_txt_data_1D(char *file_name)
{
    FILE *fp;
    double *data;
    int i;
    fp = fopen(file_name, "r");
    
    if (fp==NULL){
        printf("errors when opening 1D geo file!\n");
    }

    data = (double*) calloc(len, sizeof (*data));
    for (i=0; i<len; i++){
        fscanf(fp, "%lf", &data[i]);
        //printf("data:%f\n", data[i]);
    }

    fclose(fp);

    return data;
}


void interp_txt_data_1D(double *data, FLOAT x, FLOAT y, FLOAT z, FLOAT *a)
{
    int i;

    if (x==0){
        a[0] = data[0];
    }
    if (x>0 && x<x_end){
        i = x/dx;
        a[0] = data[i] + (data[i+1]-data[i])*(x-dx*i)/dx;
    }
    if (x==x_end){
        a[0]=data[len-1];
    }


}
