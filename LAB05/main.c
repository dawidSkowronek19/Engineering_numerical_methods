#include <stdio.h>
#include <stdlib.h>
#define n_x 200
#define n_y 90
#define Delta 0.01
#define i1 50
#define j1 55
#define rho 1
#define mi 1
#define IT_MAX 20000

void Gauss_Siedel(double **psi, double **zeta, double Omega)
{
    for (int i=1; i<n_x; i++)
    {
        for (int j=1; j<n_y; j++)
        {
            if (i>i1 || j>j1)
            {
                psi[i][j]=(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-Delta*Delta*zeta[i][j])/4.0;
                zeta[i][j]=(zeta[i+1][j]+zeta[i-1][j]+zeta[i][j+1]+zeta[i][j-1])/4.0-Omega*rho*((psi[i][j+1]-psi[i][j-1])*(zeta[i+1][j]-zeta[i-1][j])-(psi[i+1][j]-psi[i-1][j])*(zeta[i][j+1]-zeta[i][j-1]))/(16.0*mi);
            }
        }
    }
}
double Q_wy(double Q_we)
{
    return Q_we*(n_y*n_y*n_y-j1*j1*j1-3.0*j1*n_y*n_y+3.0*j1*j1*n_y)/(n_y*n_y*n_y);
}
void Boundary_conditions_psi(double **psi, double Q_we)
{
    for (int j=0; j<=n_y; j++)
    {
        if (j>=j1)
        {
            psi[0][j]=Q_we*((Delta*Delta*Delta*j*j*j)/3.0-(Delta*Delta*j*j)/2.0*(Delta*j1+Delta*n_y)+Delta*j*Delta*j1*Delta*n_y)/(2.0*mi);
        }
        psi[n_x][j]=Q_wy(Q_we)*((Delta*Delta*Delta*j*j*j)/3.0 - (Delta*Delta*j*j*Delta*n_y)/2.0)/(2.0*mi) + Q_we*Delta*Delta*j1*j1*(-Delta*j1+3.0*Delta*n_y)/(12.0*mi);
    }
    for (int i=1; i<n_x; i++)
    {
        psi[i][n_y]=psi[0][n_y];
        if (i<=i1)
        {
            psi[i][j1]=psi[0][j1];
        }
        if(i>=i1)
        {
            psi[i][0]=psi[0][j1];
        }
    }
    for (int j=1; j<=j1; j++)
    {
        psi[i1][j]=psi[0][j1];
    }
}
void Boundary_conditions_zeta(double **psi, double **zeta, double Q_we)
{
    for (int j=0; j<=n_y; j++)
    {
        if (j>=j1)
        {
            zeta[0][j]=Q_we*(2.0*Delta*j-Delta*j1-Delta*n_y)/(2.0*mi);   
        }
        zeta[n_x][j]=Q_wy(Q_we)*(2.0*Delta*j-Delta*n_y)/(2.0*mi);
    }
    for (int i=1; i<n_x; i++)
    {
        zeta[i][n_y]=2.0*(psi[i][n_y-1]-psi[i][n_y])/(Delta*Delta);
        if (i>=i1+1)
        {
            zeta[i][0]=2.0*(psi[i][1]-psi[i][0])/(Delta*Delta);
        }
        if(i<=i1)
        {
            zeta[i][j1]=2.0*(psi[i][j1+1]-psi[i][j1])/(Delta*Delta);
        }
    }
    for (int j=1; j<j1; j++)
    {
        zeta[i1][j]=2.0*(psi[i1+1][j]-psi[i1][j])/(Delta*Delta);
    }
    zeta[i1][j1]=(zeta[i1-1][j1]+zeta[i1][j1-1])/2.0;

}
double Gamma_eror(double **psi, double **zeta)
{
    int j2=j1+2;
    double Sum=0.0;
    for (int i=1; i<n_x; i++)
    {
        Sum+=psi[i+1][j2]+psi[i-1][j2]+psi[i][j2+1]+psi[i][j2-1]-4.0*psi[i][j2]-Delta*Delta*zeta[i][j2];
    }
    return Sum;
}
void u_v(double **psi, double Q_we)
{
    char file_name_u[100];
    char file_name_v[100];
    snprintf(file_name_u, sizeof(file_name_u), "out_u_Q=%.0lf.dat", Q_we);
    snprintf(file_name_v, sizeof(file_name_v), "out_v_Q=%.0lf.dat", Q_we);
    FILE *out_u=fopen(file_name_u, "w");
    FILE *out_v=fopen(file_name_v, "w");
    double **u=malloc((n_x+1)*sizeof(double*));
    double **v=malloc((n_x+1)*sizeof(double*));
    for (int i=0; i<=n_x; i++)
    {
        u[i]=malloc((n_y+1)*sizeof(double));
        v[i]=malloc((n_y+1)*sizeof(double));
    }
    for (int i=1; i<n_x; i++)
    {
        for (int j=1; j<n_y; j++)
        {
            u[i][j]=0.0;
            v[i][j]=0.0;
            if (i>i1 || j>j1)
            {
                u[i][j]=(psi[i][j+1]-psi[i][j])/(Delta);
                v[i][j]=-(psi[i+1][j]-psi[i][j])/(Delta);
            }
            fprintf(out_u, "%lf\t%lf\t%lf\n", i*Delta, j*Delta, u[i][j]);
            fprintf(out_v, "%lf\t%lf\t%lf\n", i*Delta, j*Delta, v[i][j]);
        }
        fprintf(out_u, "\n");
        fprintf(out_v, "\n");
    }
    for (int i=0; i<=n_x; i++)
    {
        free(u[i]);
        free(v[i]);
    }
    free(u);
    free(v);
    fclose(out_u);
    fclose(out_v);
}
int main()
{
    double Q_we[3]={-1000.0, -4000.0, 4000.0};
    char file_name_psi[100], file_name_zeta[100], file_name_error[100];
    for (int h=0; h<3; h++)
    {
        snprintf(file_name_psi, sizeof(file_name_psi), "out_psi_Q=%.0lf.dat", Q_we[h]);
        snprintf(file_name_zeta, sizeof(file_name_zeta), "out_zeta_Q=%.0lf.dat", Q_we[h]);
        snprintf(file_name_error, sizeof(file_name_error), "out_error_Q=%.0lf.dat", Q_we[h]);


        FILE *out_psi=fopen(file_name_psi, "w");
        FILE *out_zeta=fopen(file_name_zeta, "w");
        FILE *out_error=fopen(file_name_error, "w");
        double Omega, error;
        double **zeta=malloc((n_x+1)*sizeof(double*));
        double **psi=malloc((n_x+1)*sizeof(double*));
        for (int i=0; i<=n_x; i++)
        {
            zeta[i]=malloc((n_y+1)*sizeof(double));
            psi[i]=malloc((n_y+1)*sizeof(double));
        }
        for (int i=0; i<=n_x; i++)
        {
            for (int j =0; j <= n_y; j++)
            {
                psi[i][j]=0.0;
                zeta[i][j]=0.0;
            }
        }
        Boundary_conditions_psi(psi,Q_we[h]);
        for (int IT=1; IT<=IT_MAX; IT++)
        {
            if(IT<2000)
                Omega=0.0;
            else 
                Omega=1.0;
            Gauss_Siedel(psi, zeta,Omega);
            Boundary_conditions_zeta(psi,zeta,Q_we[h]);
            error=Gamma_eror(psi,zeta);
            fprintf(out_error, "%d\t%lf\n", IT, error);
        }

        for (int i=0; i<=n_x; i++)
        {
            for (int j=0; j<=n_y; j++)
            {
                fprintf(out_psi, "%lf\t%lf\t%lf\n", i*Delta, j*Delta, psi[i][j]);
                fprintf(out_zeta, "%lf\t%lf\t%lf\n", i*Delta, j*Delta, zeta[i][j]);
            }
            fprintf(out_psi, "\n");
            fprintf(out_zeta, "\n");
        }
        u_v(psi, Q_we[h]);
        for (int i=0; i<=n_x; i++)
        {
            free(zeta[i]);
            free(psi[i]);
        }
        free(zeta);
        free(psi);
        fclose (out_error);
        fclose(out_psi);
        fclose(out_zeta);
    }
    return 0;
}