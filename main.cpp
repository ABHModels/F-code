#include <stdlib.h>
#include <iostream>
#include <algorithm>

using namespace std;

#include "def.h"
#include "functions.cpp"

#define imax 100
#define pmax 62

int main(int argc, char *argv[])
{
	long double spin2;
	long double iobs;
	long double isco, xin;
	
	long double robs_i, robs_f; 
	long double pstep;
	
	long double pobs, pold, plar, psma, pdiff, wei, ptol=1.e-05;
	long double robscur, robs[imax], robstemp[imax];
	
	int aux_e=0, aim; 
	long double aux0, aux1, aux2;
	long double auxr1, auxr2, auxxm, auxxl, auxz, auxp1, auxp2, auxp3, auxpp, auxz1;
	long double auxr[imax];
	
	long double traced[5], tracestat[pmax];
	
	long double gerem[imax]; double val;
	long double rmissold, rmiss, germtol = 1.0e-05;
	int fin, count;
	
	long double gmm[2], gfac[pmax], gmin=10.0, gmax=0.0, gmino, gmaxo, gstar[pmax], gcur, gdiff, gtol = 1.e-05;
	long double trf, jac;
	long double valley1, valley2;
    long double F[imax], F_const, temp, sum, varl, varr, varm;
    long double f_col, Mbh, Mdd, D, A1, A2;
    long double eta, z[imax], theta[imax];

	int i, j, k, m;
	int ii;
    
    int iii, jjj;
	char filename_o[128];
	
	FILE *foutput, *finput;
	
	long double pgdata[100][2];
   	long double spin_vs_defpar[30][30][2];
    long double spinarray[]={-0.998, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.5804822, 0.6381477, 0.6800805, 0.7142378, 0.7435757, 0.7695684, 0.7930734, 0.8146393, 0.8346417, 0.8533508, 0.8709681, 0.8876485, 0.9035142, 0.9186634, 0.9331764, 0.9471198, 0.9605495, 0.9735132, 0.9860516, 0.9982};
    
    finput = fopen("grids/grid_Mdl.dat", "r");
    for(i = 0; i < 30; i++)
        for(j = 0; j < 30; j++)
            fscanf(finput, "%Lf %Lf", &spin_vs_defpar[i][j][0], &spin_vs_defpar[i][j][1]);
    fclose(finput);
	
    printf("%Lf\n", spin_vs_defpar[29][29][1]);
	/* ----- Set free parameters ----- */
    
    for(iii = 0; iii < 30; iii++) {
        for(jjj = 0; jjj < 30; jjj++) {
	
            spin = spin_vs_defpar[iii][jjj][0];      /* spin parameter */
            Mdl = spin_vs_defpar[iii][jjj][1];
            /* Deformation Parameters */
                 
            epsi3 = 0.0;
            //a13 = spin_vs_defpar[iii][jjj][1];
            a13 = 0.0;
            a22 = 0.0;//spin_vs_defpar[iii][jjj][1];
            a52 = 0.0;
            
            spin2 = spin*spin;
            
            /* ----- Set inner radius of the disk ----- */
            //isco = atof(argv[6]);
            find_isco(15, isco);
            eta = 1.0 - Energy_func(isco);
            
            printf("ISCO = %Lf, eta = %Lf\n", isco, eta);
            /* ----- Set computational parameters ----- */
            robs_i = isco;
            robs_f = 1000000.;
            /*
            GAULEG
            COMPUTE THE EMISSION RADII
            */
            auxr1 = 1.0/sqrt(robs_f);
            auxr2 = 1.0/sqrt(robs_i);
            m = (imax + 1.)/2.;
            auxxm = 0.5*(auxr2 + auxr1);
            auxxl = 0.5*(auxr2 - auxr1);
            
            for (i = 1; i <= m; i++) {
                auxz = cos( Pi*(i - 0.25)/(imax + 0.5) );
                do {
                    auxp1 = 1.0;
                    auxp2 = 0.0;
                    for (j = 1; j <= imax; j++) {
                        auxp3 = auxp2;
                        auxp2 = auxp1;
                        auxp1 = ((2.0*j - 1.0)*auxz*auxp2 - (j - 1.0)*auxp3)/j;
                    }
                    auxpp = imax*(auxz*auxp1 - auxp2)/(auxz*auxz - 1.0);
                    auxz1 = auxz;
                    auxz  = auxz1 - auxp1/auxpp;
                } while ( fabs(auxz - auxz1) > 3.0e-14 );
                auxr[i - 1]    = auxxm - auxxl*auxz;
                auxr[imax - i] = auxxm + auxxl*auxz;
            }

            for (i = 0; i <= imax - 1; i++)
                robs[imax - i - 1] = 1.0/(auxr[i]*auxr[i]);
            printf("xin = %Lf\n", robs[0]);
            
            /*for (i = 0; i <= imax - 1; i++)
            {
                z[i] = 3.0*Mdl*(1.0 - sqrt(isco/robs[i]))/eta;
                theta[i] = atan2(robs[i], z[i]);
                printf("theta = %Lf\n", theta[i]);
            }*/
            
            sprintf(filename_o,"data_Mdl/DimlessF_a%.05Le.Mdl_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,Mdl,a13,a22,a52);
            //printf("%s\n",filename_o);
            foutput = fopen(filename_o,"w");
            //fprintf(foutput, "%.10Le %.10Le\n", isco, 0.0);
            //computation between ISCO and robs[0];
            {
                temp = Energy_func(robs[0]) - Omega_func(robs[0])*Lz_func(robs[0]);
                temp = temp*temp*sqrt(-G_func(robs[0]));
                F_const = - derOmega_func(robs[0])/temp;
                
                sum = 0.0;
                
                for(ii = 0; ii <= imax - 1; ii++) {
                    robstemp[ii] = isco + ii*(robs[0] - isco)/(imax - 1);
                }
                for(j = 0; j < imax - 1; j++) {
                    varl = (Energy_func(robstemp[j]) - Omega_func(robstemp[j])*Lz_func(robstemp[j]))*derLz_func(robstemp[j]);
                    varr = (Energy_func(robstemp[j+1]) - Omega_func(robstemp[j+1])*Lz_func(robstemp[j+1]))*derLz_func(robstemp[j+1]);
                    varm = 0.5*(varl + varr);
                    sum = sum + varm*(robstemp[j+1] - robstemp[j]);
                }
                F[0] = F_const*sum;
            }
            
            //computation between robs[0] and robs[99];
            for (i = 1; i <= imax - 1; i++) {
                temp = Energy_func(robs[i]) - Omega_func(robs[i])*Lz_func(robs[i]);
                temp = temp*temp*sqrt(-G_func(robs[i]));
                F_const = - derOmega_func(robs[i])/temp;
                //printf("%d: Omega_func = %Le\n", i, Omega_func(robs[i]));
                for(ii = 0; ii <= imax - 1; ii++) {
                    robstemp[ii] = robs[i - 1] + ii*(robs[i] - robs[i - 1])/(imax - 1);
                }
                for(j = 0; j < imax - 1; j++) {
                    varl = (Energy_func(robstemp[j]) - Omega_func(robstemp[j])*Lz_func(robstemp[j]))*derLz_func(robstemp[j]);
                    varr = (Energy_func(robstemp[j+1]) - Omega_func(robstemp[j+1])*Lz_func(robstemp[j+1]))*derLz_func(robstemp[j+1]);
                    varm = 0.5*(varl + varr);
                    sum = sum + varm*(robstemp[j+1] - robstemp[j]);
                }
                F[i] = F_const*sum;
            }
            
            //fprintf(foutput, "{");
            for(i = 0; i <= imax - 1; i++)
            {
                if(F[imax - i - 1] < 0)
                    printf("ERROR!!!\t%s: %d: F[%d] = %Le\n", filename_o, imax - i - 1, imax - i - 1, F[imax - i -1]);
            }
                 
            for (i = 0; i <= imax - 1; i++)
                fprintf(foutput, "%.10Le %.10Le\n", robs[imax - i - 1], F[imax - i - 1]);
            fclose(foutput);
        
        }
    }
    
	return 0;
    }
