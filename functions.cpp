void metric(long double r, long double th, long double mn[][4])
{
    long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
    long double g_tt, g_rr, g_thth, g_pp, g_tp;
    
    t1 = cos(th);
    t2 = spin * spin;
    t3 = r * r;
    t4 = pow(t3, 0.2e1);
    t5 = t3 * t4;
    t1 = t2 * pow(t1, 0.2e1);
    t6 = (t1 + t3) * r + epsi3;
    t7 = a22 + t3;
    t8 = sin(th);
    t9 = t3 + t2;
    t8 = pow(t8, 0.2e1);
    t10 = r * t3 + a13;
    t11 = -t2 * r * t7 * t8 + t10 * t9;
    t12 = -0.2e1 * r + t9;
    t13 = a52 + t3;
    t14 = 0.1e1 / r;
    t15 = t2 * a22;
    t16 = 0.1e1 / t12;
    t13 = 0.1e1 / t13;
    t11 = pow(t11, -0.2e1);
    
    g_tt = t6 * (t2 * pow(t7, 0.2e1) * t8 + 0.2e1 * r * t4 - t4 * t9) * r * t11;
    g_rr = t6 * r * t16 * t13;
    g_thth = t14 * epsi3 + t1 + t3;
    g_pp = (pow(t9, 0.2e1) * pow(t10, 0.2e1) - t2 * t5 * t12 * t8) * t6 * t8 * t11 * t14;
    g_tp = -(0.2e1 * t5 + t3 * (a13 * (a22 + t2) + ((r * a22 + a13) * r + t15) * r) + t15 * a13) * spin * t6 * t8 * t11;
    
    mn[0][0] = g_tt;
    mn[0][3] = g_tp;
    mn[1][1] = g_rr;
    mn[2][2] = g_thth;
    mn[3][0] = mn[0][3];
    mn[3][3] = g_pp;
}

void metric_rderivatives(long double r, long double th, long double dmn[][4])
{
    long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23;
    long double dgttdr, dgtpdr, dgppdr;
    
    t1 = cos(th);
    t2 = spin * spin;
    t3 = r * r;
    t4 = pow(t3, 0.2e1);
    t5 = t3 * t4;
    t1 = t2 * pow(t1, 0.2e1);
    t6 = (t1 + t3) * r + epsi3;
    t7 = a22 + t3;
    t8 = sin(th);
    t9 = t3 + t2;
    t10 = -0.2e1 * r + t9;
    t8 = pow(t8, 0.2e1);
    t11 = r * t3 + a13;
    t12 = t9 * t11;
    t13 = -t2 * r * t7 * t8 + t12;
    t1 = t1 / 0.3e1 + t3;
    t14 = 0.2e1 / 0.3e1 * t2;
    t15 = 0.3e1 / 0.5e1 * t2;
    t16 = (t15 + t3) * r + 0.2e1 / 0.5e1 * a13;
    t15 = r * t16 - t15 * (a22 / 0.3e1 + t3) * t8;
    t13 = 0.1e1 / t13;
    t17 = pow(t13, 0.2e1);
    t18 = 0.3e1 * t1;
    t19 = a22 + t2;
    t20 = r * a22;
    t21 = t2 * a22;
    t22 = 0.1e1 / 0.2e1;
    t23 = t17 * t8;
    t9 = -t2 * t5 * t10 * t8 + pow(t9, 0.2e1) * pow(t11, 0.2e1);
    t11 = 0.1e1 / r;
    
    dgttdr = t17 * ((t6 * (0.10e2 * r * t15 * t13 - 0.1e1) - t18 * r) * (-t2 * pow(t7, 0.2e1) * t8 + t10 * t4) - 0.6e1 * t6 * (t3 * ((-0.5e1 / 0.3e1 + r) * r + t14) - t14 * t7 * t8) * t3);
    dgtpdr = t23 * spin * ((0.20e2 * t6 * t15 * t13 - 0.6e1 * t1) * (t22 * (t3 * (a13 * t19 + ((a13 + t20) * r + t21) * r) + t21 * a13) + t5) - 0.12e2 * t6 * ((t19 / 0.6e1 + t3 / 0.3e1) * a13 + t4 + t20 * (0.5e1 / 0.12e2 * t3 + t2 / 0.4e1)) * r);
    dgppdr = 0.10e2 * t23 * t6 * (-t15 * t9 * t11 * t13 + t12 * t16 - 0.4e1 / 0.5e1 * t4 * ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t2) * t2 * t8) + t23 * t9 * t11 * (-t11 * t6 + t18);
    
    dmn[0][0] = dgttdr;
    dmn[0][3] = dgtpdr;
    dmn[3][0] = dmn[0][3];
    dmn[3][3] = dgppdr;
}

long double Omega_func(long double r)
{
    long double g001, g031, g331;
    long double dg[4][4];
    long double Omega;
    long double th;
    
    th = Pi/2.0;
    metric_rderivatives(r, th, dg);
    g001 = dg[0][0];
    g031 = dg[0][3];
    g331 = dg[3][3];
    
    Omega = (-g031 + sqrt(g031*g031 - g001*g331))/g331;
    return Omega;
}

long double Energy_func(long double r)
{
    long double Omega, aux, E;
    long double g00, g03, g33;
    long double g[4][4];
    long double th;
    
    th = Pi/2.0;
    
    metric(r, th, g);
    g00 = g[0][0];
    g03 = g[0][3];
    g33 = g[3][3];
    
    Omega = Omega_func(r);
    aux = sqrt(-g00 - 2.0*g03*Omega - g33*Omega*Omega);
    E = - (g00 + g03*Omega)/aux;
    return E;
}

long double Lz_func(long double r)
{
    long double Omega, aux, Lz;
    long double g00, g03, g33;
    long double g[4][4];
    long double th;
    
    th = Pi/2.0;
    
    metric(r, th, g);
    g00 = g[0][0];
    g03 = g[0][3];
    g33 = g[3][3];
    
    Omega = Omega_func(r);
    aux = sqrt(-g00 - 2.0*g03*Omega - g33*Omega*Omega);
    Lz = (g03 + g33*Omega)/aux;
    return Lz;
}

long double G_func(long double r)
{
    long double g[4][4];
    long double G;
    long double th;
    
    th = Pi/2.0;
    metric(r, th, g);
    G = -((g[0][3]*g[0][3])/g[3][3] - g[0][0])*g[1][1]*g[3][3];
    return G;
}

long double derOmega_func(long double r)
{
    long double dr, temp, derOmega;
    long double th;
    
    th = Pi/2.0;
    
    dr = 0.0001*r;
    temp = r + dr;
    dr = temp - r;
    
    derOmega = 0.5*(Omega_func(r + dr) - Omega_func(r - dr))/dr;
    return derOmega;
}

long double derLz_func(long double r)
{
    long double dr, temp, derLz;
    long double th;
    
    th = Pi/2.0;
    
    dr = 0.0001*r;
    temp = r + dr;
    dr = temp - r;
    
    derLz = 0.5*(Lz_func(r + dr) - Lz_func(r - dr))/dr;
    return derLz;
}

void find_isco(long double z1, long double& isco)
{
    
    int i, j, casenum=1, stop=0, count=0;
    long double detol = 1.0e-8;
    long double rll,rul,rinit=z1,rnew,rold,rstep=1.e-5;
    long double mn[4][4],dmn[4][4];
    long double mnp[4][4],mnm[4][4],dmnp[4][4],dmnm[4][4];
    long double ep,em,eold,enew,omp,omm,omold,omnew;
    long double epsq,emsq;
    long double deold,denew;
    long double dr=1.0e-5;
    long double sqspin=spin*spin;
    
    if(spin>0.)
        rll = 1.+sqrt(1.-sqspin);
    else if(spin<0.)
        rll = 1.+sqrt(1.-sqspin);
    else
        rll = 1.;
    
    if(a13>-5.7){
        rul = rinit; rold = rul;
        metric(rold+dr,Pi/2.,mnp);
        metric(rold-dr,Pi/2.,mnm);
        metric_rderivatives(rold+dr,Pi/2.,dmnp);
        metric_rderivatives(rold-dr,Pi/2.,dmnm);
        omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3];
        omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3];
        ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(-mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp);
        em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(-mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm);
        deold = 0.5*(ep-em)/dr;
        
        do{
            count++;
            if(count>40){
                printf("No convergence after %i iterations.\n",count);
                break;
            }
            
            rnew = (rll+rul)/2.;
            metric(rnew+dr,Pi/2.,mnp);
            metric(rnew-dr,Pi/2.,mnm);
            metric_rderivatives(rnew+dr,Pi/2.,dmnp);
            metric_rderivatives(rnew-dr,Pi/2.,dmnm);
            omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3];
            omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3];
            ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(-mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp);
            em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(-mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm);
            denew = 0.5*(ep-em)/dr;
            
            if(fabs(denew)<fabs(detol)) {
                //printf("denew = %Le, deold = %Le\n",denew,deold);
                stop = 1;
            }
            else if((denew*deold)>0.0) {
                if(rnew<rold)
                    rul = rnew;
                else if(rnew>rold)
                    rll = rnew;
                else
                    printf("rold=rnew? rold = %Le, rnew = %Le",rold,rnew);
            }
            else if((denew*deold)<0.0) {
                if(rnew<rold)
                    rll = rnew;
                else if(rnew>rold)
                    rul = rnew;
                else
                    printf("rold=rnew? rold = %Le, rnew = %Le",rold,rnew);
            }
            else
                printf("Compare enew and eold. eold = %Le, enew = %Le\n",deold,denew);
            
            
            rold = rnew;
            omold = omnew;
            deold = denew;
            
        }while(stop==0);
    }
    else {
        rold = 4.0; // ISCO radius for these special cases is never above this, so we start here. -SN
        
        do{
            rnew = rold-rstep;
            
            if(rnew<rll){
                printf("Couldn't find ISCO? rnew = %Le, rll = %Le\n",rnew,rll);
                break;
            }
            
            metric(rnew+dr,Pi/2.,mnp);
            metric(rnew-dr,Pi/2.,mnm);
            metric_rderivatives(rnew+dr,Pi/2.,dmnp);
            metric_rderivatives(rnew-dr,Pi/2.,dmnm);
            omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3];
            omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3];
            epsq = -mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp;
            emsq = -mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm;
            
            if(epsq>0. && emsq>0.){
                ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(epsq);
                em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(emsq);
                denew = 0.5*(ep-em)/dr;
                
                if(fabs(denew)<fabs(detol)) {
                    //printf("denew = %Le, deold = %Le\n",denew,deold);
                    stop = 1;
                    //break;
                }
                else {
                    //printf("%Le\n",rnew);
                    rold = rnew;
                }
            }
            else{
                printf("epsq = %Le, emsq = %Le, rnew = %Le\n",epsq,emsq,rnew);
                stop = 1;
            }
            
        }while(stop==0);
    }
    
    
    if(stop==1)
        isco=rnew;
}

void vertical_stability(long double z1, long double isco, long double iobs_deg)
{
    int i;
    long double r, dr, th, dth;
    long double sqspin=spin*spin;
    long double en, elle, omega;
    long double mn[4][4], dmn[4][4];
    long double tdot, Veff[3], sqfreq_th=10.;
    
    char filename_o[128];
    FILE *fout;
    
    r = 10.*isco;
    dr = 0.1*isco;
    dth = 0.01;
    
    sprintf(filename_o,"vertistab_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
    fout = fopen(filename_o,"w");
    
    do{
        th = Pi/2.;
        metric(r,th,mn);
        metric_rderivatives(r,th,dmn);
        omega = (-dmn[0][3]+sqrt(dmn[0][3]*dmn[0][3]-dmn[0][0]*dmn[3][3]))/dmn[3][3];
        en = - (mn[0][0]+mn[0][3]*omega)/sqrt(-mn[0][0]-2.*mn[0][3]*omega-mn[3][3]*omega*omega);
        elle =  (mn[0][3]+mn[3][3]*omega)/sqrt(-mn[0][0]-2.*mn[0][3]*omega-mn[3][3]*omega*omega);
        
        
        for(i=0;i<3;i++) {
            th = Pi/2. - 1.*dth + (long double)i*dth;
            metric(r,th,mn);
            Veff[i] = (en*en*mn[3][3] + 2.*en*elle*mn[0][3] + elle*elle*mn[0][0])/(mn[0][3]*mn[0][3] - mn[0][0]*mn[3][3]);
        }
        
        th = Pi/2.;
        metric(r,th,mn);
        tdot = (en*mn[3][3] + elle*mn[0][3])/(mn[0][3]*mn[0][3] - mn[0][0]*mn[3][3]);
        sqfreq_th = -0.5*(Veff[2]+Veff[0]-2*Veff[1])/(dth*dth*mn[2][2]*tdot*tdot);
        
        fprintf(fout,"%.10Le %.10Le\n",r,sqrt(sqfreq_th));
        
        r-=dr;
        
    }while(sqfreq_th>0. && r>isco);
    
    fclose(fout);
}


