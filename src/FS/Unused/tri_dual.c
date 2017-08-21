#include <stdio.h>
#include <math.h>


#define SOLVER(a, b, c, d, e, f, x, y) {		\
	x = -((b)*(f) - (d)*(e))/((a)*(d) - (b)*(c));	\
	y =  ((a)*(f) - (c)*(e))/((a)*(d) - (b)*(c));	\
    }



void get_tri_center(double *p1, double *p2, double *p3)
{
    /* circumcenter */
    int k;
    double dp1[2], dp2[2];
    double mid1[2], mid2[2];
    double xc, yc;

    for (k = 0; k < 2; k++) {
	dp1[k] = p2[k] - p1[k];
	dp2[k] = p3[k] - p1[k];
    }

    for (k = 0; k < 2; k++) {
	mid1[k] = .5 * (p2[k] + p1[k]);
	mid2[k] = .5 * (p3[k] - p1[k]);
    }

    SOLVER(-dp1[1], dp2[1], dp1[0], -dp2[0], 
	   -mid1[0] + mid2[0], -mid1[1] + mid2[1], 
	   xc, yc);
    
    
    /* cpc = mid1 + s(1) * [ -dp1(2), dp1(1) ]; */
    /* cr = norm ( p(ct(1),:) - cpc ); */

    /* pc(it,:) = cpc; */
    /* r(it,1) = cr; */


    /* test in tirangle */


    /* longest edge */


    
    return;
}


int
main(int argc, char *argv[])
{
    double x[3][2] = {
	{0, 0}, 
	{1, 0}, 
	{0, sqrt(3)}
    };

    get_tri_center(x[0], x[1], x[2]); 


    return 0;
}
