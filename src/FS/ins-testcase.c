/*
 *  Problems defination
 *  */
/* ********************** */
/* *** define equations * */
/* ********************** */

#if TEST_CASE == 1
/* *****************************************************************
 * test case 1:
 * Stokes equations, analytic solution: u:P3, p:P2, unsteady
 *
 * ****************************************************************/
#warning Test case 1: P3P2t
void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
    u[0] = -0.2e0 * (t + 0.5e0) * (double) (y - 2) * (0.1e1 + 0.1e0 * z);
    u[1] = -0.2e1 * (double) (2 - y) * z + t * x - z * (t * x + 0.2e1 - 0.2e1 * t) - t * x * z;
    u[2] = 0.5e0 - z * z;
}

void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p) {
    p[0] = 0.18e1 - 0.6e-1 * t * x * z;
}

void func_gradp(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradp) {
    gradp[0] = -0.6e-1 * t * z;
    gradp[1] = 0;
    gradp[2] = -0.6e-1 * t * x;
}

void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu) {
    gradu[0] = 0;
    gradu[1] = -0.2e0 * (t + 0.5e0) * (0.1e1 + 0.1e0 * z);
    gradu[2] = -0.2e-1 * (t + 0.5e0) * (double) (y - 2);
    gradu[3] = t - 0.2e1 * t * z;
    gradu[4] = 0.2e1 * z;
    gradu[5] = -0.6e1 + (double) (2 * y) - 0.2e1 * t * x + 0.2e1 * t;
    gradu[6] = 0;
    gradu[7] = 0;
    gradu[8] = -0.2e1 * z;
}

void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f) {
    f[0] = -0.2e0 * (double) (y - 2) * (0.1e1 + 0.1e0 * z) - 0.6e-1 * t * z;
    f[1] = (double) x - z * (double) (x - 2) - (double) x * z;
    f[2] = (double) (2 * nu) - 0.6e-1 * t * (double) x;
}

void func_conv(FLOAT x, FLOAT y, FLOAT z, FLOAT *conv) {
    conv[0] = -0.2e0 * (t + 0.5e0) * (0.1e1 + 0.1e0 * z) * (-0.2e1 * (double) (2 - y) * z + t * x - z * (t * x + 0.2e1 - 0.2e1 * t) - t * x * z) - 0.2e-1 * (t + 0.5e0) * (double) (y - 2) * (0.5e0 - z * z);
    conv[1] = -0.2e0 * (t - 0.2e1 * t * z) * (t + 0.5e0) * (double) (y - 2) * (0.1e1 + 0.1e0 * z) + 0.2e1 * z * (-0.2e1 * (double) (2 - y) * z + t * x - z * (t * x + 0.2e1 - 0.2e1 * t) - t * x * z) + (-0.6e1 + (double) (2 * y) - 0.2e1 * t * x + 0.2e1 * t) * (0.5e0 - z * z);
    conv[2] = -0.2e1 * z * (0.5e0 - z * z);
}

void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g1) {
    g1[0] = -0.18e1 + 0.6e-1 * t * x * z;
    g1[1] = -0.2e0 * nu * (t + 0.5e0) * (0.1e1 + 0.1e0 * z);
    g1[2] = -0.2e-1 * nu * (t + 0.5e0) * (double) (y - 2);
}

void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g2) {
    g2[0] = nu * (t - 2 * t * z);
    g2[1] = (double) (2 * nu * z) - 0.18e1 + 0.6e-1 * (double) t * (double) x * (double) z;
    g2[2] = nu * (-6 + 2 * y - 2 * t * x + 2 * t);
}

void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g3) {
    g3[0] = 0;
    g3[1] = 0;
    g3[2] = -(double) (2 * nu * z) - 0.18e1 + 0.6e-1 * (double) t * (double) x * (double) z;
}

void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *T) {
    T[0] = 0.;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0.;
}

#elif TEST_CASE == 2
/* *****************************************************************
 * test case 2:
 * Stokes equations, analytic solution: u:P3, p:P2, steady
 *
 * ****************************************************************/
#warning Test case 2: P3P2
void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
u[0] = -0.10e0 * (double) (y - 2) * (0.1e1 + 0.1e0 * z);
u[1] = -0.2e1 * (double) (2 - y) * z + x - 0.2e1 * x * z;
u[2] = 0.5e0 - z * z;
}

void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p) {
p[0] = 0.18e1 - 0.6e-1 * x * z;
}

void func_gradp(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradp) {
gradp[0] = -0.6e-1 * z;
gradp[1] = 0;
gradp[2] = -0.6e-1 * x;
}

void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu) {
gradu[0] = 0;
gradu[1] = -0.10e0 - 0.10e-1 * z;
gradu[2] = -0.10e-1 * y + 0.20e-1;
gradu[3] = 0.1e1 - 0.2e1 * z;
gradu[4] = 0.2e1 * z;
gradu[5] = -0.4e1 + 0.2e1 * y - (double) (2 * x);
gradu[6] = 0;
gradu[7] = 0;
gradu[8] = -0.2e1 * z;
}

void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f) {
f[0] = -0.6e-1 * z;
f[1] = 0;
f[2] = (double) (2 * nu) - 0.6e-1 * x;
}

void func_conv(FLOAT x, FLOAT y, FLOAT z, FLOAT *conv) {
conv[0] = (-0.10e0 - 0.10e-1 * z) * (-0.2e1 * (double) (2 - y) * z + x - 0.2e1 * x * z) + (-0.10e-1 * (double) y + 0.20e-1) * (0.5e0 - z * z);
conv[1] = -0.10e0 * (0.1e1 - 0.2e1 * z) * (double) (y - 2) * (0.1e1 + 0.1e0 * z) + 0.2e1 * z * (-0.2e1 * (double) (2 - y) * z + x - 0.2e1 * x * z) + (-0.4e1 + (double) (2 * y) - 0.2e1 * x) * (0.5e0 - z * z);
conv[2] = -0.2e1 * z * (0.5e0 - z * z);
}

void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g1) {
g1[0] = -0.18e1 + 0.6e-1 * x * z;
g1[1] = nu * (-0.10e0 - 0.10e-1 * z);
g1[2] = nu * (-0.10e-1 * y + 0.20e-1);
}

void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g2) {
g2[0] = nu * (1 - 2 * z);
g2[1] = (double) (2 * nu * z) - 0.18e1 + 0.6e-1 * (double) x * (double) z;
g2[2] = nu * (-4 + 2 * y - 2 * x);
}

void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g3) {
g3[0] = 0;
g3[1] = 0;
g3[2] = -(double) (2 * nu * z) - 0.18e1 + 0.6e-1 * (double) x * (double) z;
}

void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *T) {
    T[0] = 0.;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0.;
}

#else
#  warning test case error!!!
#endif



/* ------------------------------------------------------------
 *    
 *    B.C. map
 *    
 * ------------------------------------------------------------ */

int
my_bc_map(int bctype)
{
    switch (bctype) {
    default:
	return DIRICHLET;
    }
}
