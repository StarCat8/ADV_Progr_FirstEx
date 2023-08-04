//CREA DISTRIBUZIONE CERCHIO

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "build/funzioni.h"
#include "build/allvars.h"

float (*f)(float);
float (*A)(float);
float (*F)(float, ... );

int main(int argc, char *argv[]){
    if(argc != 2) {
        printf("Errore nel numero di argomenti passati. E' richiesto solamente il numero di punti da generare\n");
        return 1;
    }
    int i, j, k;
    long int uu;
    float RHO0, RHOMAX, massPerPoint;
    long int *dumm = &uu;
    float x, phi, theta, RHO, RHO_cumulative, r, dr, X, Y, Z;
    float *ang1, *ang2, *rad;
    char paramFile[] = "param.txt";
    readParamFloat(paramFile);
    readParamInt(paramFile);
    int Npnt = atoi(argv[1]);
    float floatNpnt = (float)Npnt;
    massPerPoint = total_mass/floatNpnt;
    printf("la %f %f %f %d %fal", uniform_sphere_radius, isothermal_sphere_radius, circle_radius, Npnt, x);
    RHO0 = findRHO0(total_mass, Rmax, concentration);


    /*crea cerchio */
        FILE *distrC;
        char cartellaCIR[] = "CIR_data";
        mkdir(cartellaCIR, 0777);
        distrC = fopen("CIR_data/distrCerchio.txt", "w+");
        FILE *denProfileC;
        FILE *denProfileC_cumulative;
        FILE *denProfileCLOG;
        denProfileC = fopen("CIR_data/densProCerchio.txt", "w+");
        denProfileC_cumulative = fopen("CIR_data/densProCerchio_cumulative.txt", "w+");
        denProfileCLOG = fopen("CIR_data/densProCerchioLOG.txt", "w+");
        /*angolo*/
        f=&anglePHI;
        ang1 = creaPunti(Npnt, f, dumm);
        /*raggio*/
        f=&radialC;
        rad = creaPunti(Npnt, f, dumm);
        /*print to file*/
        for ( i = 0; i < Npnt; i++ ) {
            X=*(rad + i) * circle_radius * cos(*(ang1 + i));
            Y=*(rad + i) * circle_radius * sin(*(ang1 + i));
            fprintf(distrC, "%f %f 1\n", X, Y);
        }
        fclose(distrC);
        A=&areaCerchio;
        r=0; 
        dr=0.01;
        RHO_cumulative = 0;
        while(r<=1){
            RHO = density(rad, A, r, (r+dr), Npnt) * massPerPoint;
            RHO_cumulative += RHO * (A(r+dr)-A(r));
            fprintf(denProfileC, "%f %f 1\n", r * circle_radius, RHO);
            fprintf(denProfileC_cumulative, "%f %f 1\n", r * circle_radius, RHO_cumulative);
            fprintf(denProfileCLOG, "%f %f 1\n", r * circle_radius, log10(RHO));
            r=r+dr;
        }
        fclose(denProfileC);
        fclose(denProfileC_cumulative);
        fclose(denProfileCLOG);
        free(ang1);
        free(rad);
    

    /*crea sfera uniforme*/
        FILE *distrS;
        char cartellaUNI[] = "UNI_data";
        mkdir(cartellaUNI, 0777);
        distrS = fopen("UNI_data/distrSfera.txt", "w+");
        FILE *denProfileS, *denProfileS_cumulative;
        FILE *denProfileSLOG;
        denProfileS = fopen("UNI_data/densProSfera.txt", "w+");
        denProfileS_cumulative = fopen("UNI_data/densProSfera_cumulative.txt", "w+");
        denProfileSLOG = fopen("UNI_data/densProSferaLOG.txt", "w+");
        /*angoli*/
        f=&anglePHI;
        ang1 = creaPunti(Npnt, f, dumm);
        f=&angleTHETA;
        ang2 = creaPunti(Npnt, f, dumm);
        /*raggio*/
        f=&radialS;
        rad = creaPunti(Npnt, f, dumm);
        /*print to file*/
        for ( i = 0; i < Npnt; i++ ) {
            X=*(rad + i) * uniform_sphere_radius * cos(*(ang1 + i)) * sin(*(ang2+i));
            Y=*(rad + i) * uniform_sphere_radius * sin(*(ang1 + i))*sin(*(ang2+i));
            Z=*(rad + i) * uniform_sphere_radius * cos(*(ang2+i));
            fprintf(distrS, "%f %f %f 1\n", X, Y, Z);
        }
        fclose(distrS);
        A=&volSfera;
        r=0;
        dr=0.01;
        RHO_cumulative = 0;
        while(r<=1){
            RHO = density(rad, A, r, (r+dr), Npnt) * massPerPoint;
            RHO_cumulative += RHO * (A(r+dr)-A(r));
            fprintf(denProfileS_cumulative, "%f %f 1\n", r * uniform_sphere_radius, RHO_cumulative);
            fprintf(denProfileS, "%f %f 1\n", r * uniform_sphere_radius, RHO);
            fprintf(denProfileSLOG, "%f %f 1\n", r * uniform_sphere_radius, log10(RHO));
            r=r+dr; 
        }
        fclose(denProfileS);
        fclose(denProfileS_cumulative);
        fclose(denProfileSLOG);
        free(ang1);
        free(ang2);
        free(rad);

    /*crea sfera isoterma*/
        FILE *distrIsoTherm;
        char cartellaISO[] = "ISO_data";
        mkdir(cartellaISO, 0777);
        distrIsoTherm = fopen("ISO_data/distrSferaISO.txt", "w+");
        FILE *denProfileIsoTherm, *denProfileIsoTherm_cumulative;
        FILE *denProfileIsoThermLOG;
        denProfileIsoTherm = fopen("ISO_data/densProSferaISO.txt", "w+");
        denProfileIsoTherm_cumulative = fopen("ISO_data/densProSferaISO_cumulative.txt", "w+");
        denProfileIsoThermLOG = fopen("ISO_data/densProSferaISOLOG.txt", "w+");
        /*angoli*/
        f=&anglePHI;
        ang1 = creaPunti(Npnt, f, dumm);
        f=&angleTHETA;
        ang2 = creaPunti(Npnt, f, dumm);
        /*raggio*/
        f=&radialIsoThermal;
        rad = creaPunti(Npnt, f, dumm);
        /*print to file*/
        for ( i = 0; i < Npnt; i++ ) {
            X=*(rad + i) * isothermal_sphere_radius * cos(*(ang1 + i)) * sin(*(ang2+i));
            Y=*(rad + i) * isothermal_sphere_radius * sin(*(ang1 + i))*sin(*(ang2+i));
            Z=*(rad + i) * isothermal_sphere_radius * cos(*(ang2+i));
            fprintf(distrIsoTherm, "%f %f %f 1\n", X, Y, Z);
        }
        fclose(distrIsoTherm);
        A=&volSfera;
        r=0;
        dr=0.01;
        RHO_cumulative = 0;
        while(r<=1){
            RHO = density(rad, A, r, (r+dr), Npnt) * massPerPoint;
            RHO_cumulative += RHO * (A(r+dr)-A(r));
            fprintf(denProfileIsoTherm, "%f %f 1\n", r * isothermal_sphere_radius, RHO);
            fprintf(denProfileIsoTherm_cumulative, "%f %f 1\n", r * isothermal_sphere_radius, RHO_cumulative);
            fprintf(denProfileIsoThermLOG, "%f %f 1\n", log10(r * isothermal_sphere_radius), log10(RHO));
            r=r+dr;
        }
        fclose(denProfileIsoTherm);
        fclose(denProfileIsoThermLOG);
        free(ang1);
        free(ang2);
        free(rad);
        
    /*Navarro Frenk and White*/
    FILE *distProfileNFW;
    FILE *denProfileNFW;
    FILE *denProfileNFW_cumulative;
    FILE *denProfileNFWLOG;
    char cartellaNFW[] = "NFW_data";
    mkdir(cartellaNFW, 0777);
    denProfileNFW = fopen("NFW_data/densProNFW.txt", "w+");
    denProfileNFW_cumulative = fopen("NFW_data/densProNFW_cumulative.txt", "w+");
    denProfileNFWLOG = fopen("NFW_data/densProNFWLOGz.txt", "w+");
    distProfileNFW = fopen("NFW_data/distribProNFW.txt", "w+");
    float *NFWpoints = malloc(sizeof(float) * Npnt);       // allocate memory for the array of points of the distribution
    char NFWdist[] = "NFW_data/NFWdistribc.txt";
    RHOMAX = NFW( RHO0 , r_min, Rmax/concentration);
    printf("\n ->%f<- \n", RHOMAX);
    F = &NFW;
    rejection(1, RHOMAX, Npnt, Rmax, dumm, RHO0, NFWdist, F, NFWpoints, Rmax/concentration, r_min);
    makeHist(NFWdist, 20);
    f=&anglePHI;
    ang1 = creaPunti(Npnt, f, dumm);
    f=&angleTHETA;
    ang2 = creaPunti(Npnt, f, dumm);
    for ( i = 0; i < Npnt; i++ ) {
        X=*(NFWpoints + i)  * cos(*(ang1 + i)) * sin(*(ang2+i));
        Y=*(NFWpoints + i) *  sin(*(ang1 + i))*sin(*(ang2+i));
        Z=*(NFWpoints + i) *  cos(*(ang2+i));
        fprintf(distProfileNFW, "%f %f %f 1\n", X, Y, Z);
    }
        A=&volSfera;
        r=1;
        dr=Rmax/1000;
        RHO_cumulative = 0;
        while(r<=Rmax){
            RHO = density(NFWpoints, A, r, (r+dr), Npnt) * massPerPoint;
            RHO_cumulative += RHO * (A(r+dr)-A(r));
            if(RHO > 0){
                fprintf(denProfileNFW, "%f %f 1\n", r , RHO);
                fprintf(denProfileNFW_cumulative, "%f %f 1\n", r , RHO_cumulative);
                fprintf(denProfileNFWLOG, "%f %f 1\n", log10(r), log10(RHO));
            }
            r=r+dr;
        }
        fclose(denProfileNFW);
        fclose(denProfileNFW_cumulative);
        fclose(denProfileNFWLOG);
        free(ang1);
        free(ang2);
        free(NFWpoints);
    
    return 0;
}


