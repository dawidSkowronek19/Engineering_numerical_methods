#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


#define n_x 128
#define n_y 128
#define Delta 0.2
#define TOL 1e-8
# define PI 3.1415926535897932384626433832795028841971
void boundary_conditions(double **V)
{
    double x_max=n_x*Delta;
    double y_max=n_y*Delta;
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            V[i][j]=0.0;
        }
    }
    for (int i=0; i<=n_x; i++)
    {
        V[i][n_y]=-sin(2.0*PI*i*Delta/x_max);
        V[i][0]=sin(2.0*PI*i*Delta/x_max);
    }
    for (int j=0; j<=n_y;j++)
    {
        V[0][j]=sin(PI*Delta*j/y_max);
        V[n_x][j]=sin(PI*Delta*j/y_max);
    }
}
void relax(int k, double **V)
{
    for (int i=k; i<=n_x-k; i+=k)
    {
        for (int j=k; j<=n_y-k; j+=k)
        {
            V[i][j]=(V[i+k][j]+V[i-k][j]+V[i][j+k]+V[i][j-k])/4.0;
        }
    }
}
void interpolation(double **V, int k)
{
    for (int i=0; i<=n_x-k; i+=k)
    {
        for (int j=0; j<=n_y-k; j+=k)
        {
            V[i+k/2][j+k/2]=(V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k])/4.0;
            if (i!=n_x-k)
            {
                V[i+k][j+k/2]=(V[i+k][j]+V[i+k][j+k])/2.0;
            }
            if (j!=n_y-k)
            {
                V[i+k/2][j+k]=(V[i][j+k]+V[i+k][j+k])/2.0;
            }
            if (j!=0)
            {
                V[i+k/2][j]=(V[i][j]+V[i+k][j])/2.0;
            }
            if (i!=0)
            {
                V[i][j+k/2]=(V[i][j]+V[i][j+k])/2.0;
            }
        }
    }
}
double STOP(double **V, int k)
{
    double S=0.0;
    for (int i=0; i<=n_x-k; i+=k)
    {
        for (int j=0; j<=n_y-k; j+=k)
        {
            S+=(Delta*Delta*k*k/2.0)*((V[i+k][j]-V[i][j]+V[i+k][j+k]-V[i][j+k])*(V[i+k][j]-V[i][j]+V[i+k][j+k]-V[i][j+k])+(V[i][j+k]-V[i][j]+V[i+k][j+k]-V[i+k][j])*(V[i][j+k]-V[i][j]+V[i+k][j+k]-V[i+k][j]))/(4.0*k*k*Delta*Delta);
        }
    }
    return S;
}
int main()
{
    bool first_iteration=true;
    int k;
    int it=0;
    double S_old;
    double S_new;
    double epsilon;
    double **V=malloc((n_x+1)*sizeof(double*));
    for (int i=0; i<=n_x; i++)
    {
        V[i]=malloc((n_y+1)*sizeof(double));
    }

    boundary_conditions(V);

    for (int iteration=0; iteration<5; iteration++)
    {
        k=pow(2,4-iteration);
        char filename1[64];
        char filename2[64];

        sprintf(filename1, "error_k=%d.out", k);
        sprintf(filename2, "potential_k=%d.out", k);

        FILE *error = fopen(filename1, "w");
        FILE *potential = fopen(filename2, "w");

        do
        {
            it++;
            if (first_iteration==true)
            {
                relax(k,V);
                S_old=STOP(V, k);
                fprintf(error, "%d\t%lf\n", it,S_old);
                epsilon=1.0;
                first_iteration=false;
            }
            else
            {
                relax(k,V);
                S_new=STOP(V, k);
                fprintf(error, "%d\t%lf\n", it,S_new);
                epsilon=fabs((S_new-S_old)/S_old);
                S_old=S_new;
            }
        }while (epsilon>=TOL);
        if (k!=1)
        {
            interpolation(V,k);
        }
        first_iteration=true;
        for (int i=0; i<=n_x; i+=k)
        {
            for (int j=0; j<=n_y; j+=k)
            {
                fprintf(potential, "%lf\t%lf\t%lf\n", i*Delta, j*Delta, V[i][j]);
            }
            fprintf(potential, "\n");
        }
        fclose(error);
        fclose(potential);
    }

    for (int i=0; i<=n_x; i++)
    {
        free(V[i]);
    }
    free(V);
    return 0;
}