#define PI 3.141592653589793
#define e 2.718281828459045
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum);
float * creaPunti(int Npnt, float (*f)(float), long int *dumm);
float anglePHI(float phi);
float angleTHETA(float costheta);
float radialC(float x);
float radialS(float x);
float radialIsoThermal(float x);
float areaCerchio(float r);
float volSfera(float r);
float density(float *punti, float (*A)(float), float dist1, float dist2, int Npnt);
void readParamFloat(char *name);
void readParamInt(char *name);
void rejection(float ymin, float ymax, int numpunti, float intervallo, long *dumm, float rho0, char *name, float(*f)(float, float, float), float *point, float Rs, float Rmin);
float NFW(float rho0, float r, float rs);
float findRHO0(float mass, float Rmax, float c);
float totalNFW_mass(float rho0, float Rmax, float rs);
void makeHist(char *name, int Nbin);