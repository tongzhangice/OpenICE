/* #include <stdio.h> */
/* #include <math.h> */
/* #include <stdlib.h> */
/* #include <string.h> */
/* #include "ins.h" */


#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482 


/* Basic parameters */
static double L;
static double n = 3;
static double rho;
static double g;
static double spyr;
static double Sx;
static double Sz;
static double D;

static double cx, cy, cbx, cby, ct, A;
static double Z;
static double alpha; 

static double Sp;
static double Su;
//static double Sw;
static double Smu;
static double St;


/* Variables */
//static double t;
#  define PRES_SCALING 1e5


void 
func_init_params(double Lx0, double alpha0)
{
    L = Lx0;			// [m]
    alpha = alpha0;		// [degree]

    /* Basic parameters */
    n = 3.		; // Glen’s parameter
    rho = 910.		; // [kg/mˆ3]
    g = 9.81		; // [m/secˆ2]
    A = 1.e-16		; // [Paˆ(-n)aˆ(-1)]
    spyr = 31556926.	; // [sec/a]
    Sx = L * 1e3	; // [m]
    Sz = 1e3		; // [m]
    D = Sz/Sx	        ; // aspect ratio

    /* Scales */
    Sp    = rho*g*Sz		;	// [Pa=kg/m/secˆ2]
    Su    = pow(2*Sp, n)*Sx*A	;	// [m/a]
    Smu   = 0.5*pow(A, (-1./n)) *
	pow(Su/Sx, (1.-n)/n)	;	// [Pa secˆ2]
    St    = Sx/Su		;	// [a]
    /* Ssig  = Sp/Sz		;	// stress scale [J=Pa/m] */
    /* Ssigb = Sp			;	// bound.stress scale [Pa] */


    phgPrintf("* Delta %e\n", D);
    phgPrintf("* Sx   %e\n", Sx);
    phgPrintf("* Sz   %e\n", Sz);;
    phgPrintf("* Su    %e\n", Su);
    phgPrintf("* St    %e\n", St);
    phgPrintf("* Smu   %e\n", Smu);

    /* Parameters */
    cx = 1e-9 * Su;
    cy = 1e-9 * Su;
    cbx = 0;
    cby = 0;
    ct = 1e-8 / St;
    A = 1e-16;

    Z = 1.;
    L = Lx0;

    return;
}



/*
 * Evaluate order:
 *
 * 0) d, h, s, b
 * 1) u, v, w
 * 2) mu
 * 3) p
 * 4) f 
 *
 *
 * */
static void
func_analytic(double x, double y, double z, int type, double *value)
{
    double s, sx, sy, sxx, sxy, syx, syy, sxxx, sxxy, sxyx, sxyy, syxx, syxy, syyx, syyy,
	b, bx, by, bxx, bxy, byx, byy,
	q, qx, qy, qxx, qxy, qyx, qyy,  
	h, hx, hy, hxx, hxy, hyx, hyy,  hxxx, hxxy, hxyx, hxyy, hyxx, hyxy, hyyx, hyyy, 
	d, dx, dy, dz, dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz,
	dxxx, dxxy, dxxz, dxyx, dxyy, dxyz, dxzx, dxzy, dxzz,
	dyxx, dyxy, dyxz, dyyx, dyyy, dyyz, dyzx, dyzy, dyzz,
	dzxx, dzxy, dzxz, dzyx, dzyy, dzyz, dzzx, dzzy, dzzz;
    double u, ux, uy, uz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz,
	v, vx, vy, vz, vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz, 
	w, wx, wy, wz, wxx, wxy, wxz, wyx, wyy, wyz, wzx, wzy, wzz;	
    double mu, mux, muy, muz, p, px, py, pz;
    double fx, fy, fz;


    /* double _Length_ = 80.; */
    /* double _alpha_ = .5; */
    /* func_init_params(_Length_, _alpha_ * PI / 180.); */


    /* x,y,z, unit: km */
    /* x = atof(argv[1]); */
    /* y = atof(argv[2]); */
    /* z = atof(argv[3]); */
    /* t = atof(argv[4]); */

    /* printf("x %30.15e\n", x); */
    /* printf("y %30.15e\n", y); */
    /* printf("z %30.15e\n", z); */


    /*
     * s, b, h, q, d,           unit: km
     * [sbhqd][xyz],            unit: 1
     * [sbhqd][xyz][xyz],       unit: 1/km
     * [sbhqd][xyz][xyz],       unit: 1/km
     * [sbhqd][xyz][xyz][xyz],  unit: 1/km/km
     * */
    s = -x * tan(alpha) + Z * sin(0.2e1 * PI * x / L) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) / 0.2e1;
    b = -x * tan(alpha) + Z * sin(0.2e1 * PI * x / L) * sin(0.2e1 * PI * y / L) / 0.2e1 - Z;

    if (type == 0) {
	*(value++) = s;
	*(value)   = b;
	//printf("s");
	return;
    }

    h = Z * sin(0.2e1 * PI * x / L) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) / 0.2e1 - Z * sin(0.2e1 * PI * x / L) * sin(0.2e1 * PI * y / L) / 0.2e1 + Z;
    hx = Z * cos(0.2e1 * PI * x / L) * PI / L * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) - Z * cos(0.2e1 * PI * x / L) * PI / L * sin(0.2e1 * PI * y / L);
    hy = Z * sin(0.2e1 * PI * x / L) * cos(0.2e1 * PI * y / L) * PI / L * (0.1e1 - exp(-ct * t)) - Z * sin(0.2e1 * PI * x / L) * cos(0.2e1 * PI * y / L) * PI / L;
    hxx = -0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L);
    hxy = 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) - 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L);
    hyx = 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) - 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L);
    hyy = -0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L);
    hxxx = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L);
    hxxy = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L);
    hxyx = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L);
    hxyy = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L);
    hyxx = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L);
    hyxy = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L);
    hyyx = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L);
    hyyy = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t)) + 0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L);

    sx = -tan(alpha) + Z * cos(0.2e1 * PI * x / L) * PI / L * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    sy = Z * sin(0.2e1 * PI * x / L) * cos(0.2e1 * PI * y / L) * PI / L * (0.1e1 - exp(-ct * t));
    sxx = -0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    sxy = 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    syx = 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    syy = -0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    sxxx = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    sxxy = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    sxyx = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    sxyy = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    syxx = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    syxy = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    syyx = -0.4e1 * Z * cos(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * sin(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));
    syyy = -0.4e1 * Z * sin(0.2e1 * PI * x / L) * pow(PI, 0.3e1) * pow(L, -0.3e1) * cos(0.2e1 * PI * y / L) * (0.1e1 - exp(-ct * t));

    b = -x * tan(alpha) + Z * sin(0.2e1 * PI * x / L) * sin(0.2e1 * PI * y / L) / 0.2e1 - Z;
    bx = -tan(alpha) + Z * cos(0.2e1 * PI * x / L) * PI / L * sin(0.2e1 * PI * y / L);
    by = Z * sin(0.2e1 * PI * x / L) * cos(0.2e1 * PI * y / L) * PI / L;
    bxx = -0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L);
    bxy = 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L);
    byx = 0.2e1 * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L);
    byy = -0.2e1 * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L);

    q = cx * Z * cos(0.2e1 * PI * x / L) * cos(0.2e1 * PI * y / L) * exp(-ct * t) / 0.2e1;
    qx = -cx * Z * sin(0.2e1 * PI * x / L) * PI / L * cos(0.2e1 * PI * y / L) * exp(-ct * t);
    qy = -cx * Z * cos(0.2e1 * PI * x / L) * sin(0.2e1 * PI * y / L) * PI / L * exp(-ct * t);
    qxx = -0.2e1 * cx * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L) * exp(-ct * t);
    qxy = 0.2e1 * cx * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L) * exp(-ct * t);
    qyx = 0.2e1 * cx * Z * sin(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * sin(0.2e1 * PI * y / L) * exp(-ct * t);
    qyy = -0.2e1 * cx * Z * cos(0.2e1 * PI * x / L) * PI * PI * pow(L, -0.2e1) * cos(0.2e1 * PI * y / L) * exp(-ct * t);

    d = (s - z) / h;
    dx = sx / h - (s - z) *  pow( h,  (-2)) * hx;
    dy = sy / h - (s - z) *  pow( h,  (-2)) * hy;
    dz = -1 / h;
    dxx = sxx / h - 2 * sx *  pow( h,  (-2)) * hx + 2 * (s - z) *  pow( h,  (-3)) * hx * hx - (s - z) *  pow( h,  (-2)) * hxx;
    dxy = sxy / h - sx *  pow( h,  (-2)) * hy - sy *  pow( h,  (-2)) * hx + 2 * (s - z) *  pow( h,  (-3)) * hx * hy - (s - z) *  pow( h,  (-2)) * hxy;
    dxz =  pow( h,  (-2)) * hx;
    dyx = sxy / h - sx *  pow( h,  (-2)) * hy - sy *  pow( h,  (-2)) * hx + 2 * (s - z) *  pow( h,  (-3)) * hx * hy - (s - z) *  pow( h,  (-2)) * hxy;
    dyy = syy / h - 2 * sy *  pow( h,  (-2)) * hy + 2 * (s - z) *  pow( h,  (-3)) * hy * hy - (s - z) *  pow( h,  (-2)) * hyy;
    dyz =  pow( h,  (-2)) * hy;
    dzx =  pow( h,  (-2)) * hx;
    dzy =  pow( h,  (-2)) * hy;
    dzz = 0;
    dxxx = sxxx / h - 3 * sxx *  pow( h,  (-2)) * hx + 6 * sx *  pow( h,  (-3)) * hx * hx - 3 * sx *  pow( h,  (-2)) * hxx - 6 * (s - z) *  pow( h,  (-4)) *  pow( hx,  3) + 6 * (s - z) *  pow( h,  (-3)) * hx * hxx - (s - z) *  pow( h,  (-2)) * hxxx;
    dxxy = sxxy / h - sxx *  pow( h,  (-2)) * hy - 2 * sxy *  pow( h,  (-2)) * hx + 4 * sx *  pow( h,  (-3)) * hx * hy - 2 * sx *  pow( h,  (-2)) * hxy + 2 * sy *  pow( h,  (-3)) * hx * hx - 6 * (s - z) *  pow( h,  (-4)) * hx * hx * hy + 4 * (s - z) *  pow( h,  (-3)) * hx * hxy - sy *  pow( h,  (-2)) * hxx + 2 * (s - z) *  pow( h,  (-3)) * hxx * hy - (s - z) *  pow( h,  (-2)) * hxxy;
    dxxz = -2 *  pow( h,  (-3)) * hx * hx +  pow( h,  (-2)) * hxx;
    dxyx = sxxy / h - sxx *  pow( h,  (-2)) * hy - 2 * sxy *  pow( h,  (-2)) * hx + 4 * sx *  pow( h,  (-3)) * hx * hy - 2 * sx *  pow( h,  (-2)) * hxy + 2 * sy *  pow( h,  (-3)) * hx * hx - 6 * (s - z) *  pow( h,  (-4)) * hx * hx * hy + 4 * (s - z) *  pow( h,  (-3)) * hx * hxy - sy *  pow( h,  (-2)) * hxx + 2 * (s - z) *  pow( h,  (-3)) * hxx * hy - (s - z) *  pow( h,  (-2)) * hxxy;
    dxyy = sxyy / h - 2 * sxy *  pow( h,  (-2)) * hy + 2 * sx *  pow( h,  (-3)) * hy * hy - sx *  pow( h,  (-2)) * hyy - syy *  pow( h,  (-2)) * hx + 4 * sy *  pow( h,  (-3)) * hx * hy - 2 * sy *  pow( h,  (-2)) * hxy - 6 * (s - z) *  pow( h,  (-4)) * hx * hy * hy + 4 * (s - z) *  pow( h,  (-3)) * hxy * hy + 2 * (s - z) *  pow( h,  (-3)) * hx * hyy - (s - z) *  pow( h,  (-2)) * hxyy;
    dxyz = -2 *  pow( h,  (-3)) * hx * hy +  pow( h,  (-2)) * hxy;
    dxzx = -2 *  pow( h,  (-3)) * hx * hx +  pow( h,  (-2)) * hxx;
    dxzy = -2 *  pow( h,  (-3)) * hx * hy +  pow( h,  (-2)) * hxy;
    dxzz = 0;
    dyxx = sxxy / h - sxx *  pow( h,  (-2)) * hy - 2 * sxy *  pow( h,  (-2)) * hx + 4 * sx *  pow( h,  (-3)) * hx * hy - 2 * sx *  pow( h,  (-2)) * hxy + 2 * sy *  pow( h,  (-3)) * hx * hx - 6 * (s - z) *  pow( h,  (-4)) * hx * hx * hy + 4 * (s - z) *  pow( h,  (-3)) * hx * hxy - sy *  pow( h,  (-2)) * hxx + 2 * (s - z) *  pow( h,  (-3)) * hxx * hy - (s - z) *  pow( h,  (-2)) * hxxy;
    dyxy = sxyy / h - 2 * sxy *  pow( h,  (-2)) * hy + 2 * sx *  pow( h,  (-3)) * hy * hy - sx *  pow( h,  (-2)) * hyy - syy *  pow( h,  (-2)) * hx + 4 * sy *  pow( h,  (-3)) * hx * hy - 2 * sy *  pow( h,  (-2)) * hxy - 6 * (s - z) *  pow( h,  (-4)) * hx * hy * hy + 4 * (s - z) *  pow( h,  (-3)) * hxy * hy + 2 * (s - z) *  pow( h,  (-3)) * hx * hyy - (s - z) *  pow( h,  (-2)) * hxyy;
    dyxz = -2 *  pow( h,  (-3)) * hx * hy +  pow( h,  (-2)) * hxy;
    dyyx = sxyy / h - 2 * sxy *  pow( h,  (-2)) * hy + 2 * sx *  pow( h,  (-3)) * hy * hy - sx *  pow( h,  (-2)) * hyy - syy *  pow( h,  (-2)) * hx + 4 * sy *  pow( h,  (-3)) * hx * hy - 2 * sy *  pow( h,  (-2)) * hxy - 6 * (s - z) *  pow( h,  (-4)) * hx * hy * hy + 4 * (s - z) *  pow( h,  (-3)) * hxy * hy + 2 * (s - z) *  pow( h,  (-3)) * hx * hyy - (s - z) *  pow( h,  (-2)) * hxyy;
    dyyy = syyy / h - 3 * syy *  pow( h,  (-2)) * hy + 6 * sy *  pow( h,  (-3)) * hy * hy - 3 * sy *  pow( h,  (-2)) * hyy - 6 * (s - z) *  pow( h,  (-4)) *  pow( hy,  3) + 6 * (s - z) *  pow( h,  (-3)) * hy * hyy - (s - z) *  pow( h,  (-2)) * hyyy;
    dyyz = -2 *  pow( h,  (-3)) * hy * hy +  pow( h,  (-2)) * hyy;
    dyzx = -2 *  pow( h,  (-3)) * hx * hy +  pow( h,  (-2)) * hxy;
    dyzy = -2 *  pow( h,  (-3)) * hy * hy +  pow( h,  (-2)) * hyy;
    dyzz = 0;
    dzxx = -2 *  pow( h,  (-3)) * hx * hx +  pow( h,  (-2)) * hxx;
    dzxy = -2 *  pow( h,  (-3)) * hx * hy +  pow( h,  (-2)) * hxy;
    dzxz = 0;
    dzyx = -2 *  pow( h,  (-3)) * hx * hy +  pow( h,  (-2)) * hxy;
    dzyy = -2 *  pow( h,  (-3)) * hy * hy +  pow( h,  (-2)) * hyy;
    dzyz = 0;
    dzzx = 0;
    dzzy = 0;
    dzzz = 0;


    /*
     * u, v,  unit: m/a
     * */
    u = cx * (0.1e1 - pow(d, 0.4e1)) + cbx / h;
    v = (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) / h;
    w = u*h*dx+v*h*dy;

    if (type == 1) {
	*(value++) = u;
	*(value++) = v;
	*(value) = w;
	//printf("u");
	return;
    }

    ux = -0.4e1 * cx * pow(d, 0.3e1) * dx - cbx * pow(h, -0.2e1) * hx;
    uy = -0.4e1 * cx * pow(d, 0.3e1) * dy - cbx * pow(h, -0.2e1) * hy;
    uz = -0.4e1 * cx * pow(d, 0.3e1) * dz;
    uxx = -0.12e2 * cx * d * d * dx * dx - 0.4e1 * cx * pow(d, 0.3e1) * dxx + 0.2e1 * cbx * pow(h, -0.3e1) * hx * hx - cbx * pow(h, -0.2e1) * hxx;
    uxy = -0.12e2 * cx * d * d * dx * dy - 0.4e1 * cx * pow(d, 0.3e1) * dxy + 0.2e1 * cbx * pow(h, -0.3e1) * hx * hy - cbx * pow(h, -0.2e1) * hxy;
    uxz = -0.12e2 * cx * d * d * dx * dz - 0.4e1 * cx * pow(d, 0.3e1) * dxz;
    uyx = -0.12e2 * cx * d * d * dx * dy - 0.4e1 * cx * pow(d, 0.3e1) * dxy + 0.2e1 * cbx * pow(h, -0.3e1) * hx * hy - cbx * pow(h, -0.2e1) * hxy;
    uyy = -0.12e2 * cx * d * d * dy * dy - 0.4e1 * cx * pow(d, 0.3e1) * dyy + 0.2e1 * cbx * pow(h, -0.3e1) * hy * hy - cbx * pow(h, -0.2e1) * hyy;
    uyz = -0.12e2 * cx * d * d * dy * dz - 0.4e1 * cx * pow(d, 0.3e1) * dyz;
    uzx = -0.12e2 * cx * d * d * dx * dz - 0.4e1 * cx * pow(d, 0.3e1) * dxz;
    uzy = -0.12e2 * cx * d * d * dy * dz - 0.4e1 * cx * pow(d, 0.3e1) * dyz;
    uzz = -0.12e2 * cx * d * d * dz * dz - 0.4e1 * cx * pow(d, 0.3e1) * dzz;
    vx = (-0.4e1 * cy * pow(d, 0.3e1) * dx - qx + qx * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dx) / h - (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * pow(h, -0.2e1) * hx;
    vy = (-0.4e1 * cy * pow(d, 0.3e1) * dy - qy + qy * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dy) / h - (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * pow(h, -0.2e1) * hy;
    vz = (-0.4e1 * cy * pow(d, 0.3e1) * dz + 0.4e1 * q * pow(d, 0.3e1) * dz) / h;
    vxx = (-0.12e2 * cy * d * d * dx * dx - 0.4e1 * cy * pow(d, 0.3e1) * dxx - qxx + qxx * pow(d, 0.4e1) + 0.8e1 * qx * pow(d, 0.3e1) * dx + 0.12e2 * q * d * d * dx * dx + 0.4e1 * q * pow(d, 0.3e1) * dxx) / h + (-0.2e1 * (-0.4e1 * cy * pow(d, 0.3e1) * dx - qx + qx * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dx) * hx - (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * hxx) * pow(h, -0.2e1) + 0.2e1 * (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * pow(h, -0.3e1) * hx * hx;
    vxy = (-0.12e2 * cy * d * d * dx * dy - 0.4e1 * cy * pow(d, 0.3e1) * dxy - qxy + qxy * pow(d, 0.4e1) + 0.4e1 * qx * pow(d, 0.3e1) * dy + 0.4e1 * qy * pow(d, 0.3e1) * dx + 0.12e2 * q * d * d * dx * dy + 0.4e1 * q * pow(d, 0.3e1) * dxy) / h + (-(cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * hxy - (-0.4e1 * cy * pow(d, 0.3e1) * dx - qx + qx * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dx) * hy - (-0.4e1 * cy * pow(d, 0.3e1) * dy - qy + qy * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dy) * hx) * pow(h, -0.2e1) + 0.2e1 * (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * pow(h, -0.3e1) * hx * hy;
    vxz = (-0.12e2 * cy * d * d * dx * dz - 0.4e1 * cy * pow(d, 0.3e1) * dxz + 0.4e1 * qx * pow(d, 0.3e1) * dz + 0.12e2 * q * d * d * dx * dz + 0.4e1 * q * pow(d, 0.3e1) * dxz) / h - (-0.4e1 * cy * pow(d, 0.3e1) * dz + 0.4e1 * q * pow(d, 0.3e1) * dz) * pow(h, -0.2e1) * hx;
    vyx = (-0.12e2 * cy * d * d * dx * dy - 0.4e1 * cy * pow(d, 0.3e1) * dxy - qxy + qxy * pow(d, 0.4e1) + 0.4e1 * qx * pow(d, 0.3e1) * dy + 0.4e1 * qy * pow(d, 0.3e1) * dx + 0.12e2 * q * d * d * dx * dy + 0.4e1 * q * pow(d, 0.3e1) * dxy) / h + (-(cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * hxy - (-0.4e1 * cy * pow(d, 0.3e1) * dx - qx + qx * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dx) * hy - (-0.4e1 * cy * pow(d, 0.3e1) * dy - qy + qy * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dy) * hx) * pow(h, -0.2e1) + 0.2e1 * (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * pow(h, -0.3e1) * hx * hy;
    vyy = (-0.12e2 * cy * d * d * dy * dy - 0.4e1 * cy * pow(d, 0.3e1) * dyy - qyy + qyy * pow(d, 0.4e1) + 0.8e1 * qy * pow(d, 0.3e1) * dy + 0.12e2 * q * d * d * dy * dy + 0.4e1 * q * pow(d, 0.3e1) * dyy) / h + (-0.2e1 * (-0.4e1 * cy * pow(d, 0.3e1) * dy - qy + qy * pow(d, 0.4e1) + 0.4e1 * q * pow(d, 0.3e1) * dy) * hy - (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * hyy) * pow(h, -0.2e1) + 0.2e1 * (cy - cy * pow(d, 0.4e1) + cby - q + q * pow(d, 0.4e1)) * pow(h, -0.3e1) * hy * hy;
    vyz = (-0.12e2 * cy * d * d * dy * dz - 0.4e1 * cy * pow(d, 0.3e1) * dyz + 0.4e1 * qy * pow(d, 0.3e1) * dz + 0.12e2 * q * d * d * dy * dz + 0.4e1 * q * pow(d, 0.3e1) * dyz) / h - (-0.4e1 * cy * pow(d, 0.3e1) * dz + 0.4e1 * q * pow(d, 0.3e1) * dz) * pow(h, -0.2e1) * hy;
    vzx = (-0.12e2 * cy * d * d * dx * dz - 0.4e1 * cy * pow(d, 0.3e1) * dxz + 0.4e1 * qx * pow(d, 0.3e1) * dz + 0.12e2 * q * d * d * dx * dz + 0.4e1 * q * pow(d, 0.3e1) * dxz) / h - (-0.4e1 * cy * pow(d, 0.3e1) * dz + 0.4e1 * q * pow(d, 0.3e1) * dz) * pow(h, -0.2e1) * hx;
    vzy = (-0.12e2 * cy * d * d * dy * dz - 0.4e1 * cy * pow(d, 0.3e1) * dyz + 0.4e1 * qy * pow(d, 0.3e1) * dz + 0.12e2 * q * d * d * dy * dz + 0.4e1 * q * pow(d, 0.3e1) * dyz) / h - (-0.4e1 * cy * pow(d, 0.3e1) * dz + 0.4e1 * q * pow(d, 0.3e1) * dz) * pow(h, -0.2e1) * hy;
    vzz = (-0.12e2 * cy * d * d * dz * dz - 0.4e1 * cy * pow(d, 0.3e1) * dzz + 0.12e2 * q * d * d * dz * dz + 0.4e1 * q * pow(d, 0.3e1) * dzz) / h;

    /*
     * w,   unit: m/a
     * */
    wx =  ux *  h *  dx +  u *  hx *  dx +  u *  h *  dxx +  vx *  h *  dy +  v *  hx *  dy +  v *  h *  dxy;
    wy =  uy *  h *  dx +  u *  hy *  dx +  u *  h *  dxy +  vy *  h *  dy +  v *  hy *  dy +  v *  h *  dyy;
    wz =  uz *  h *  dx +  u *  h *  dxz +  vz *  h *  dy +  v *  h *  dyz;
    wxx = uxx * h * dx + 2 * ux * hx * dx + 2 * ux * h * dxx + u * hxx * dx + 2 * u * hx * dxx + u * h * dxxx + vxx * h * dy + 2 * vx * hx * dy + 2 * vx * h * dxy + v * hxx * dy + 2 * v * hx * dxy + v * h * dxxy;
    wxy = uxy * h * dx + ux * hy * dx + ux * h * dxy + uy * hx * dx + u * hxy * dx + u * hx * dxy + uy * h * dxx + u * hy * dxx + u * h * dxxy + vxy * h * dy + vx * hy * dy + vx * h * dyy + vy * hx * dy + v * hxy * dy + v * hx * dyy + vy * h * dxy + v * hy * dxy + v * h * dxyy;
    wxz = uxz * h * dx + ux * h * dxz + uz * hx * dx + u * hx * dxz + uz * h * dxx + u * h * dxxz + vxz * h * dy + vx * h * dyz + vz * hx * dy + v * hx * dyz + vz * h * dxy + v * h * dxyz;
    wyx = uxy * h * dx + ux * hy * dx + ux * h * dxy + uy * hx * dx + u * hxy * dx + u * hx * dxy + uy * h * dxx + u * hy * dxx + u * h * dxxy + vxy * h * dy + vx * hy * dy + vx * h * dyy + vy * hx * dy + v * hxy * dy + v * hx * dyy + vy * h * dxy + v * hy * dxy + v * h * dxyy;
    wyy = uyy * h * dx + 2 * uy * hy * dx + 2 * uy * h * dxy + u * hyy * dx + 2 * u * hy * dxy + u * h * dxyy + vyy * h * dy + 2 * vy * hy * dy + 2 * vy * h * dyy + v * hyy * dy + 2 * v * hy * dyy + v * h * dyyy;
    wyz = uyz * h * dx + uy * h * dxz + uz * hy * dx + u * hy * dxz + uz * h * dxy + u * h * dxyz + vyz * h * dy + vy * h * dyz + vz * hy * dy + v * hy * dyz + vz * h * dyy + v * h * dyyz;
    wzx = uxz * h * dx + ux * h * dxz + uz * hx * dx + u * hx * dxz + uz * h * dxx + u * h * dxxz + vxz * h * dy + vx * h * dyz + vz * hx * dy + v * hx * dyz + vz * h * dxy + v * h * dxyz;
    wzy = uyz * h * dx + uy * h * dxz + uz * hy * dx + u * hy * dxz + uz * h * dxy + u * h * dxyz + vyz * h * dy + vy * h * dyz + vz * hy * dy + v * hy * dyz + vz * h * dyy + v * h * dyyz;
    wzz = uzz * h * dx + 2 * uz * h * dxz + u * h * dxzz + vzz * h * dy + 2 * vz * h * dyz + v * h * dyzz;

    /* [uvw]_[xyz]  unit: 1/a  */
    ux /= 1e3; uy /= 1e3; uz /= 1e3; vx /= 1e3; vy /= 1e3; vz /= 1e3; wx /= 1e3; wy /= 1e3; wz /= 1e3;

    /* [uvw]_[xyz]_[xyz]  unit: 1/a/m  */
    uxx /= 1e6; uxy /= 1e6; uxz /= 1e6; uyx /= 1e6; uyy /= 1e6; uyz /= 1e6; uzx /= 1e6; uzy /= 1e6; uzz /= 1e6; 
    vxx /= 1e6; vxy /= 1e6; vxz /= 1e6; vyx /= 1e6; vyy /= 1e6; vyz /= 1e6; vzx /= 1e6; vzy /= 1e6; vzz /= 1e6; 
    wxx /= 1e6; wxy /= 1e6; wxz /= 1e6; wyx /= 1e6; wyy /= 1e6; wyz /= 1e6; wzx /= 1e6; wzy /= 1e6; wzz /= 1e6;


    /* printf("-----------------------------\n"); */
    /* printf("u %30.15e\n", u); */
    /* printf("v %30.15e\n", v); */
    /* printf("w %30.15e\n", w); */

    /* printf("-----------------------------\n"); */
    /* printf("gu %30.15e\n", ux); */
    /* printf("gu %30.15e\n", uy); */
    /* printf("gu %30.15e\n", uz); */
    /* printf("gu %30.15e\n", vx); */
    /* printf("gu %30.15e\n", vy); */
    /* printf("gu %30.15e\n", vz); */
    /* printf("gu %30.15e\n", wx); */
    /* printf("gu %30.15e\n", wy); */
    /* printf("gu %30.15e\n", wz); */

    /* Check for singularity */
    if (sqrt(ux*ux + uy*uy + uz*uz +
	     vx*vx + vy*vy + vz*vz +
	     wx*wx + wy*wy + wz*wz) < 1e-30) {
	phgInfo(0, "exact singular at %e %e %e\n", x, y, z);
	ux = uy = uz = 0;
	vx = vy = vz = 0;
	wx = wy = wz = 0;
	uy = 1e-16;
    }

    /* mu, unit: Pa */
    mu = .5*pow(A, (-1./n)) * pow(pow(uy + vx, 0.2e1) / 0.4e1 + pow(uz + wx, 0.2e1) / 0.4e1 + pow(vz + wy, 0.2e1) / 0.4e1 - ux * vy - ux * wz - vy * wz, -0.1e1 / 0.3e1);
    mux = pow(mu, 4.) / 6. * (-0.16e2 * A * ((uy + vx) * (uxy + vxx) / 0.2e1 + (uz + wx) * (uxz + wxx) / 0.2e1 + (vz + wy) * (vxz + wxy) / 0.2e1 - uxx * vy - ux * vxy - uxx * wz - ux * wxz - vxy * wz - vy * wxz));
    muy = pow(mu, 4.) / 6. * (-0.16e2 * A * ((uy + vx) * (uyy + vxy) / 0.2e1 + (uz + wx) * (uyz + wxy) / 0.2e1 + (vz + wy) * (vyz + wyy) / 0.2e1 - uxy * vy - ux * vyy - uxy * wz - ux * wyz - vyy * wz - vy * wyz));
    muz = pow(mu, 4.) / 6. * (-0.16e2 * A * ((uy + vx) * (uyz + vxz) / 0.2e1 + (uz + wx) * (uzz + wxz) / 0.2e1 + (vz + wy) * (vzz + wyz) / 0.2e1 - uxz * vy - ux * vyz - uxz * wz - ux * wzz - vyz * wz - vy * wzz));
    //mux /= 1e3;   muy /= 1e3;   muz /= 1e3;

    /*
     * Max (mu) = 1e9 
     *
     *  */
#warning Max mu 1e9
    if (mu > 1e9) {
	mu = 1e9;
	mux = muy = muz = 0;
    }

    if (type == 2) {
	*(value) = mu;
	//printf("mu");
	return;
    }

    /* printf("-----------------------------\n"); */
    /* printf("mu %30.15e\n", 2.*mu); */

    /* p, unit: Pa */
    p = 2 * mu * ux + 2 * mu * vy - rho*g*(s - z) * 1e3;
    px = 2 * mux * ux + 2 * mu * uxx + 2 * mux * vy + 2 * mu * vxy - rho*g*sx;
    py = 2 * muy * ux + 2 * mu * uxy + 2 * muy * vy + 2 * mu * vyy - rho*g*sy;
    pz = 2 * muz * ux + 2 * mu * uxz + 2 * muz * vy + 2 * mu * vyz + rho*g;

    if (type == 3) {
	*(value) = p;
	//printf("p");
	return;
    }

    if (type == 4) {
	/* rhs, unit: Pa/m */
	fx = 2 * mux * ux + 2 * mu * uxx - px + muy * (uy + vx) + mu * (uyy + vxy) + muz * (uz + wx) + mu * (uzz + wxz);
	fy = mux * (uy + vx) + mu * (uxy + vxx) + 2 * muy * vy + 2 * mu * vyy - py + muz * (vz + wy) + mu * (vzz + wyz);
	fz = mux * (uz + wx) + mu * (uxz + wxx) + muy * (vz + wy) + mu * (vyz + wyy) + 2 * muz * wz + 2 * mu * wzz - pz - rho*g;

	*(value++) = fx;
	*(value++) = fy;
	*(value) = fz;
	//printf("f");
	return;
    }


    if (type == 5) {
	double sxx = 2 * mu * ux - p;	
	double sxy = mu * (uy + vx);		
	double sxz = mu * (uz + wx);	
	double syx = mu * (vx + uy);		
	double syy = 2 * mu * vy - p;		
	double syz = mu * (vz + wy);	
	double szx = mu * (wx + uz);	
	double szy = mu * (wy + vz);	
	double szz = 2 * mu * wz - p;		

	double rs = sqrt(1 + sx*sx + sy*sy);
	double rsx = 1/rs * (-sx * sxx - sy * sxy + sxz);
	double rsy = 1/rs * (-sx * syx - sy * syy + syz);
	double rsz = 1/rs * (-sx * szx - sy * szy + szz);
	
	*(value++) = rsx;
	*(value++) = rsy;
	*(value) = rsz;
	//printf("g");
	return;
    }


    /* printf("-----------------------------\n"); */
    /* printf("f %30.15e\n", fx); */
    /* printf("f %30.15e\n", fy); */
    /* printf("f %30.15e\n", fz); */

    return;
}

 
void func_u(double x, double y, double z, double *u) 
/* Input : [m]/LenS
 * Output: [m/a] */
{
    func_analytic(x, y, z, 1, u);
}

void func_p(double x, double y, double z, double *p) 
/* Input : [m]/LenS
 * Output: [Pa]/PresS */
{
    func_analytic(x, y, z, 3, p);
    *p /= PRES_SCALING;
}

void func_gradu(double x, double y, double z, double *gradu)
/* Unused */
{
    return;
}

void func_mu(double x, double y, double z, double *mu)
/* Input : [m]/LenS
 * Output: [?]/1  */
{
    func_analytic(x, y, z, 2, mu);
    *mu *= 2;			/* PHG */
}


void func_f(double x, double y, double z, double *f) 
/* Compensation terms, 
 * Input : [m]/LenS
 * Output: [Pa/m] */
{
    func_analytic(x, y, z, 4, f);
}


void func_g1(double x, double y, double z, double *g1) 
/* Upper surface: 
 * Input : [m]/LenS
 * Output: [Pa] */
{
    func_analytic(x, y, z, 5, g1);
}

void func_g2(double x, double y, double z, double *g1) 
/* Low surface: (unused)
 * Input : [m]/LenS
 * Output: [Pa] */
{
    return;
}

void func_g3(double x, double y, double z, double *g1) 
/* Unused
 * Input : [m]/LenS
 * Output: [Pa] */
{
    return;
}



void func_T(double x, double y, double z, double *T)
{
    T[0] = 0;
}

void func_fT(double x, double y, double z, double *fT)
{
    fT[0] = 0;
}

void func_s(double x, double y, double z, double *s_ptr)
/* Input : [m]
 * Output: [m] */
{
    double v[2];

    phgInfo(0, "   s_ptr: %x\n", s_ptr);
    x /= 1e3;
    y /= 1e3;
    z /= 1e3;
    func_analytic(x, y, z, 0, v);
    *s_ptr = v[0] * 1e3;
}

void func_b(double x, double y, double z, double *b_ptr)
/* Input : [m]
 * Output: [m] */
{
    double v[2];

    x /= 1e3;
    y /= 1e3;
    z /= 1e3;
    func_analytic(x, y, z, 0, v);
    *b_ptr = v[1] * 1e3;
}





double func_vol()
/* Output: m^2, unused */
{
    return 1;
}

void func_a(double x, double y, double *a)
/* Input : [m]/LenS
 * Output: [m/a]
 */
{
    x /= 1e3;
    y /= 1e3;

    *a = .5 * Z * ct * sin(2*PI*x / L) * sin(2*PI*y / L) * exp(-ct*t);
    *a *= 1e3;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *q) 
{
    func_a(x, y, q);
}


