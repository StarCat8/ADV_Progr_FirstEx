#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allvars.h"

#define PI 3.141592653589793
#define e 2.718281828459045
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum){
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    if (*idum < 0 || iff == 0)
        {
        iff=1;
        mj=labs(MSEED-labs(*idum));
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) 
            {
            ii=(21*i) % 55; 
            ma[ii]=mk;       
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
            }
        for (k=1;k<=4;k++)  
            for (i=1;i<=55;i++) 
                {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
                }
        inext=0;  
        inextp=31;
        *idum=1;
        }
    if (++inext == 56) inext=1;   
    if (++inextp == 56) inextp=1; 
    mj=ma[inext]-ma[inextp];      
    if (mj < MZ) mj += MBIG;      
    ma[inext]=mj;                 
    return mj*FAC;                
    }

float * creaPunti(int Npnt, float (*f)(float), long int *dumm){
    int i;
    float w;
    float *z=malloc(sizeof(float) * Npnt);
    for(i=0; i<Npnt; i++){
        z[i] = f(ran3(dumm));
        //printf("\n%f", z[i]);
    }
    return z;
    }

float anglePHI(float phi){
    phi=2*PI*phi;
    return phi;
    }
float angleTHETA(float costheta){
    float theta;
    costheta=2*costheta-1;
    theta = acos(costheta);
    return theta;
    }
float radialC(float x){
    x = pow(x, 0.5);
    return x;
    }

float radialS(float x){
    x = pow(x, (1.*1/3)); //F(r) = r^3
    return x;
    }

float radialIsoThermal(float x){
    return x;
}

float areaCerchio(float r){
    float Area = PI*r*r;
    return Area;
    }

float volSfera(float r){
        float volume = 4/3*PI*r*r*r;
        return volume;
    }

float density(float *punti, float (*A)(float), float dist1, float dist2, int Npnt){
    int i, j=0;
    float a, rho;
    a=(A(dist2)-A(dist1));
    for(i=0; i<Npnt; i++){
            if(*(punti + i)>=dist1 && *(punti + i)<=dist2){
                j++;
            }
        rho = j/a;
        }
    return rho;
    }
void readParamFloat(char *name){
    char buf[256];
    float x;
    FILE *ptrToParam;
    ptrToParam = fopen(name, "r"); 
    //printf("%s", name);
    while(fscanf(ptrToParam, "%s %f\n", buf, &x)!=EOF){
        if(!strcmp(buf, "circle_radius")){
            circle_radius = x;
        }
        if(!strcmp(buf, "uniform_sphere_radius")){
            uniform_sphere_radius = x;
        }
        if(!strcmp(buf, "IsoThermal_sphere_radius")){
            isothermal_sphere_radius = x;
        }
        if(!strcmp(buf, "total_mass")){
            total_mass = x;
        }
        if(!strcmp(buf, "concentration")){
            concentration = x;
        }
        if(!strcmp(buf, "Rmax")){
            Rmax = x;
        }
        if(!strcmp(buf, "r_min")){
            r_min = x;
        }
        }
        fclose(ptrToParam);
    }
void readParamInt(char *name){
    char buf[256];
    int x;
    FILE *ptrToParam;
    ptrToParam = fopen(name, "r");  
    while(fscanf(ptrToParam, "%s %d\n", buf, &x)!=EOF){
        if(!strcmp(buf, "Npoints")){
            Npnt = x;
        }
    }
    fclose(ptrToParam);
    }

void rejection(float ymin, float ymax, int numpunti, float intervallo, long *dumm, float rho0, char *name, float(*f)(float, float, float), float *point, float Rs, float Rmin){
	FILE *fileT;
	fileT = fopen(name, "w");
	float xtry, ytry, y; //
	int i = 0;
	while(i < numpunti){ //
		xtry = 1.*(intervallo * ran3(dumm)) + Rmin;
		ytry = 1.*((ymax-ymin)*ran3(dumm));
		y = f(rho0, xtry, Rs); //
        if(ytry < y){
            fprintf(fileT, "%f %f\n", log10(xtry), log10(ytry));
            point[i] = xtry;
            i++;
            }
        }
	fclose(fileT);
	}

float NFW(float rho0, float r, float rs){
    float x;
    x = rho0 /(( ( r / rs ) * ( 1 + r / rs ) ) *  ( 1 + r / rs ));
    return x;
}

float findRHO0(float mass, float Rmax, float c){
    float x;
    x = mass / (4 * PI * (Rmax/c) * (Rmax/c) * (Rmax/c) * (log(( (Rmax/c) + Rmax ) / (Rmax/c) ) - Rmax / ((Rmax/c) + Rmax ) ));
    return x;
}

float totalNFW_mass(float rho0, float Rmax, float rs){
    float x;
    x = 4 * PI * rho0 * rs * rs * rs * (log(( rs + Rmax ) / rs ) - Rmax / (rs + Rmax ) );
    return x;
}

void makeHist(char *name, int Nbin){ //crea istogramma da file formattato
    FILE *Hp, *ptr;
    int *HISTO = malloc(sizeof(int) * Nbin);
    int index, i;
    i=0;
    while(i<Nbin){
        HISTO[i]=0;
        i++;
    }
    float MIN, MAX, buf;
    ptr = fopen(name, "r");
    fscanf(ptr, "%f %*f %*f %*f %*d", &buf);
    MIN = 1.e33;
    MAX = -1.e33;
    //trova minimo e massimo
    while(fscanf(ptr, "%f %*f %*f %*f %*d", &buf) == 1){
        if(buf<MIN){
            MIN=buf;
        }
        if(buf>MAX){
            MAX=buf;
        }
    }
    fclose(ptr);
    ptr = fopen(name, "r");
    while(fscanf(ptr, "%f %*f %*f %*f %*d\n", &buf) == 1){
        buf = (buf - MIN)/MAX * Nbin;
        index = (int)buf;
        HISTO[index]=HISTO[index]+1;
    }
    i=0;
    Hp = fopen("NFW_data/Hist.txt", "w+");
    while(i<Nbin){
        fprintf(Hp, "%f %d\n", i*MAX/Nbin + MIN, HISTO[i]);
        i=i+1;
    }
    fclose(Hp);
    fclose(ptr);
    free(HISTO);
    }
