#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

//const long double Pi  = 3.141592653589793;
//const long double Piby2 = 1.5707963267948966192;

const long double Pi = acos(-1.0);
//const long double sigma = 5.67e-05;

long double xobs, yobs;
long double epsi3, a13, a22, a52;
long double spin, Mdl;

void metric(long double r, long double th, long double mn[][4]);
void metric_rderivatives(long double r, long double th, long double dmn[][4]);
void find_isco(long double z1, long double& r_isco);

long double Omega_func(long double r);
long double Energy_func(long double r);
long double Lz_func(long double r);
long double G_func(long double r);
long double derOmega_func(long double r);
long double derLz_func(long double r);
void constant_func(long double f_col, long double Mbh, long double Mdd, long double D, long double& A1, long double& A2);
long double integrate_func(long double g_min, long double g_max, long double g_star, long double trf, long double cosem, long double Mbh, long double lflag, long double energy, int index, long double robs[], long double F[], long double A1, long double A2);

#endif
