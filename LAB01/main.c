#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 500
#define beta 0.001
#define gamma 0.1
#define t_max 100.0
#define dt 0.1
#define u_0 1
#define TOL 1e-6

double alpha=beta*N-gamma;

double trapezoid_method(double u_old, double u_new)
{

    double u=u_old+dt*(alpha*u_old-beta*u_old*u_old+alpha*u_new-beta*u_new*u_new)/2.0;
    return u;
}
double** a_butcher()
{
    double **a=malloc(2*sizeof(double*));
    for (int i=0; i<2; i++)
    {
        a[i]=malloc(2*sizeof(double));
    }
    a[0][0]=0.25;
    a[0][1]=0.25-sqrt(3)/6.0;
    a[1][0]=0.25+sqrt(3)/6.0;
    a[1][1]=0.25;
    return a;
}
int predictor(double *U_new, double u_old)
{
    double U_old[2];
    double m[2][2];
    double F[2];
    double epsilon;
    for (int i=0; i<2; i++)
    {
        U_old[i]=u_old;
        U_new[i]=u_old;
    }
    double **a=a_butcher();
    int mi=0;
    do
    {
        mi++;
        F[0]=U_old[0]-u_old-dt*(a[0][0]*(alpha*U_old[0]-beta*U_old[0]*U_old[0])+a[0][1]*(alpha*U_old[1]-beta*U_old[1]*U_old[1]));
        F[1]=U_old[1]-u_old-dt*(a[1][0]*(alpha*U_old[0]-beta*U_old[0]*U_old[0])+a[1][1]*(alpha*U_old[1]-beta*U_old[1]*U_old[1]));

        m[0][0]=1.0-dt*a[0][0]*(alpha-2.0*beta*U_old[0]);
        m[0][1]=-dt*a[0][1]*(alpha-2.0*beta*U_old[1]);
        m[1][0]=-dt*a[1][0]*(alpha-2.0*beta*U_old[0]);
        m[1][1]=1.0-dt*a[1][1]*(alpha-2.0*beta*U_old[1]);
        U_new[0]+=(F[1]*m[0][1]-F[0]*m[1][1])/(m[0][0]*m[1][1]-m[0][1]*m[1][0]);
        U_new[1]+=(F[0]*m[1][0]-F[1]*m[0][0])/(m[0][0]*m[1][1]-m[0][1]*m[1][0]);
        U_old[0]=U_new[0];
        U_old[1]=U_new[1];
    }while(mi<20);
    for (int i=0; i<2; i++)
    {
        free(a[i]);
    }
    free(a);

    return mi;
}
double corector(double *U, double u_old)
{
    double u_new = u_old+dt*(alpha*U[0]-beta*U[0]*U[0]+alpha*U[1]-beta*U[1]*U[1])/2.0;
    return u_new;
}
int main()
{
    FILE *out_picard=fopen("out_picard.dat", "w");
    FILE *out_newton=fopen("out_newton.dat", "w");
    FILE *out_RK2=fopen("out_RK2.dat", "w");
    double *U=malloc(2*sizeof(double));
    double u_new, u_old, u_holder, epsilon, t;
    double un_new, un_old, un_holder;
    double uRK2_new, uRK2_old;
    int mi=0;
    int mi_n=0;
    int mi_RK2;
    int IT_MAX=t_max/dt;
    u_old=u_0;
    u_new=u_0;
    un_new=u_0;
    un_old=u_0;
    uRK2_new=u_0;
    uRK2_old=u_0;
    fprintf(out_picard, "%lf\t%lf\t%lf\t%d\n", 0.0, u_old, N-u_old, mi);
    fprintf(out_newton, "%lf\t%lf\t%lf\t%d\n", 0.0, un_old, N-un_old, mi);
    fprintf(out_RK2, "%lf\t%lf\t%lf\t%d\n", 0.0, uRK2_old, N-uRK2_old, mi);
    for (int i=0; i<IT_MAX; i++)
    {
        t=i*dt;
        do //Picard
        {
            mi++;
            u_holder=trapezoid_method(u_old, u_new);
            epsilon=fabs(u_holder-u_new);
            u_new=u_holder;
        }while(epsilon>TOL&&mi<=20);
        u_old=u_new;
        fprintf(out_picard, "%lf\t%lf\t%lf\t%d\n", t+dt, u_old, N-u_old, mi);
        mi=0;
        do //Newton
        {
            mi_n++;
            un_holder=un_new-(un_new-trapezoid_method(un_old, un_new))/(1.0-dt*(alpha-2.0*beta*un_new)/2.0);
            epsilon=fabs(un_holder-un_new);
            un_new=un_holder;
        }while(epsilon>TOL);
        un_old=un_new;
        fprintf(out_newton, "%lf\t%lf\t%lf\t%d\n", t+dt, un_old, N-un_old, mi_n);
        mi_n=0;
        mi_RK2=predictor(U, uRK2_old);
        uRK2_new=corector(U, uRK2_old);
        fprintf(out_RK2, "%lf\t%lf\t%lf\t%d\n", t+dt, uRK2_new, N-uRK2_new, mi_RK2);
        uRK2_old=uRK2_new;
    }
    fclose(out_picard);
    fclose(out_newton);
    fclose(out_RK2);
    free(U);
    return 0;
}