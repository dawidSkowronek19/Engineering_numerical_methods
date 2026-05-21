#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.142

#define n_x 400
#define n_y 90
#define i1 200
#define i2 210
#define j1 50
#define delta 0.01
#define D 0
#define t_max 0.6
const double sigma = 10.0*delta;
const double x_A = 0.45;
const double y_A = 0.45;

typedef struct coordinate{double x, y;}coordinate;

double starting_coundition(double x, double y)
{
    return 1.0/(2.0*PI*sigma*sigma)*exp((-(x-x_A)*(x-x_A)-(y-y_A)*(y-y_A))/(2.0*sigma*sigma));
}

coordinate** velocity_field(double **psi)
{
    coordinate **v=malloc((n_x+1)*sizeof(coordinate*));
    for (int i=0; i<=n_x; i++)
    {
        v[i]=malloc((n_y+1)*sizeof(coordinate));
    }
    for (int i=1; i<n_x; i++)
    {
        for (int j=1; j<n_y; j++)
        {
            v[i][j].x=(psi[i][j+1]-psi[i][j-1])/(2.0*delta);
            v[i][j].y=-(psi[i+1][j]-psi[i-1][j])/(2.0*delta);
        }
    }
    
    for (int i = i1; i<=i2; i++)
    {
        for (int j=0; j<=j1; j++)
        {
            v[i][j].x=0.0;
            v[i][j].y=0.0;
        }
    }
    
    for (int i=1; i<n_x; i++)
    {
        v[i][0].x=0.0;
        v[i][0].y=0.0;
        v[i][n_y].x=0.0;
        v[i][n_y].y=0.0;
    }

    for (int j=0; j<=n_y; j++)
    {
        v[0][j].x=v[1][j].x;
        v[n_x][j].x=v[n_x-1][j].x;
        v[0][j].y=v[1][j].y;
        v[n_x][j].y=v[n_x-1][j].y;
    }
    return v;
}
double v_max(coordinate **v)
{
    double max=sqrt(v[0][0].x*v[0][0].x+v[0][0].y*v[0][0].y);
    double norm;
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            norm=sqrt(v[i][j].x*v[i][j].x+v[i][j].y*v[i][j].y);
            if (norm >max)
                max=norm;
        }
    }
    return max;
}
void relaxation(double delta_t, double **u,  coordinate **v)
{
    double **u_old=malloc((n_x+1)*sizeof(double*));
    for (int i=0; i<=n_x; i++)
    {
        u_old[i]=malloc((n_y+1)*sizeof(double));
    }

    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            u_old[i][j]=u[i][j];
        }
    }
    
    for (int iteration_numb=0; iteration_numb<20; iteration_numb++)
    {
        for (int i=0; i<=n_x; i++)
        {
            for (int j=1; j<n_y; j++)
            {
                if (i>=i1&&i<=i2&&j<=j1)
                {
                    continue;
                }
                else if (i == 0)
                {
                    u[i][j]=(1.0/(1.0+2.0*D*delta_t/(delta*delta)))*(u_old[i][j]-0.5*delta_t*v[i][j].x*((u[i+1][j]-u[n_x][j]+
                        u_old[i+1][j]-u_old[n_x][j])/(2.0*delta))-0.5*delta_t*v[i][j].y*((u[i][j+1]-u[i][j-1]+
                            u_old[i][j+1]-u_old[i][j-1])/(2.0*delta))+0.5*delta_t*D*((u[i+1][j]+u[n_x][j]+u[i][j+1]+u[i][j-1])
                            /(delta*delta)+(u_old[i+1][j]+u_old[n_x][j]+u_old[i][j+1]+u_old[i][j-1]-
                            4.0*u_old[i][j])/(delta*delta)));
                }
                else if (i ==n_x)
                {
                    u[i][j]=(1.0/(1.0+2.0*D*delta_t/(delta*delta)))*(u_old[i][j]-0.5*delta_t*v[i][j].x*((u[0][j]-u[i-1][j]+
                        u_old[0][j]-u_old[i-1][j])/(2.0*delta))-0.5*delta_t*v[i][j].y*((u[i][j+1]-u[i][j-1]+
                            u_old[i][j+1]-u_old[i][j-1])/(2.0*delta))+0.5*delta_t*D*((u[0][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
                            /(delta*delta)+(u_old[0][j]+u_old[i-1][j]+u_old[i][j+1]+u_old[i][j-1]-
                            4.0*u_old[i][j])/(delta*delta)));
                }
                else
                {
                    u[i][j]=(1.0/(1.0+2.0*D*delta_t/(delta*delta)))*(u_old[i][j]-0.5*delta_t*v[i][j].x*((u[i+1][j]-u[i-1][j]+
                        u_old[i+1][j]-u_old[i-1][j])/(2.0*delta))-0.5*delta_t*v[i][j].y*((u[i][j+1]-u[i][j-1]+
                            u_old[i][j+1]-u_old[i][j-1])/(2.0*delta))+0.5*delta_t*D*((u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
                            /(delta*delta)+(u_old[i+1][j]+u_old[i-1][j]+u_old[i][j+1]+u_old[i][j-1]-
                            4.0*u_old[i][j])/(delta*delta)));
                }
            }
        }
    }

    for (int i=0; i<=n_x; i++)
    {
        free(u_old[i]);
    }
    free(u_old);
}
double density_integral(double **u)
{
    double SUM=0.0;
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            SUM+=delta*delta*u[i][j];
        }
    }
    return SUM;
}
double average_x_position(double **u)
{
    double SUM=0.0;
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            SUM+=delta*delta*u[i][j]*i*delta;
        }
    }
    return SUM;
}
void u_print(double **u, int it)
{
    char name[100];
    if (it%5==0)
    {
        snprintf(name, sizeof(name), "./u/t=%d.txt", it);
        FILE *out=fopen(name, "a");
        for (int i=0; i<=n_x; i++)
        {
            for (int j=0; j<=n_y; j++)
            {
                fprintf(out, "%lf\t%lf\t%lf\n", i*delta, j*delta, u[i][j]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
}
int main()
{
    double **psi=malloc((n_x+1)*sizeof(double*));
    coordinate **v;
    double **u=malloc((n_x+1)*sizeof(double*));
    for (int i=0; i<=n_x; i++)
    {
        psi[i]=malloc((n_y+1)*sizeof(double));
        u[i]=malloc((n_y+1)*sizeof(double));
    }
    FILE *out_psi=fopen("strumien.dat", "r");
    
    int x, y;
    double val;
    int line = 0;

    while (fscanf(out_psi, "%d %d %lf", &x, &y, &val) == 3)
    {
        psi[x][y]=val;
        line++;
    }
    
    v=velocity_field(psi);
    FILE *out_vx=fopen("out_vx.txt", "w");
    FILE *out_vy=fopen("out_vy.txt", "w");
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            fprintf(out_vx, "%lf\t%lf\t%lf\n", i*delta, j*delta, v[i][j].x);
            fprintf(out_vy, "%lf\t%lf\t%lf\n", i*delta, j*delta, v[i][j].y);
        }
        fprintf(out_vx, "\n");
        fprintf(out_vy, "\n");
    }
    fclose(out_vx);
    fclose(out_vy);
    double t=0.0;
    double delta_t = delta/(4.0*v_max(v));
    double iteration_numb=t_max/delta_t;
    FILE *out_rho=fopen("out_rho.txt", "w");
    FILE *out_averageX=fopen("out_averageX.txt","w");
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            u[i][j]=starting_coundition(i*delta, j*delta);
        }
    }
    fprintf(out_rho, "%lf\t%lf\n", t, density_integral(u));
    fprintf(out_averageX, "%lf\t%lf\n", t, average_x_position(u));
    u_print(u, 0);
    for (int i=0; i<iteration_numb; i++)
    {
        t=i*delta_t;
        relaxation(delta_t, u, v);
        u_print(u, i+1);
        fprintf(out_rho, "%lf\t%lf\n", t+delta_t, density_integral(u));
        fprintf(out_averageX, "%lf\t%lf\n", t+delta_t, average_x_position(u));

    }
    for (int i=0; i<=n_x; i++)
    {
        free(psi[i]);
        free(v[i]);
        free(u[i]);
    }
    free(psi);
    free(v);
    free(u);
    fclose(out_psi);
    fclose(out_rho);
    fclose(out_averageX);
    return 0;
}