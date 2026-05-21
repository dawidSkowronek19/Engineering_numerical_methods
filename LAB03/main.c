#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define epislon 1.01
#define Delta 0.1
#define n_x 150
#define n_y 150
#define V1 -1.0
#define V2 -1.0
#define TOL 1e-8
#define m 50
#define t_max 500.0

double rho1(double x, double y)
{
    double x_max=Delta*n_x;
    double y_max=Delta*n_y;
    double sigma_x=0.1*x_max;
    double sigma_y=0.1*y_max;

    return -exp(-(x-0.35*x_max)*(x-0.35*x_max)/(sigma_x*sigma_x)-(y-0.5*y_max)*(y-0.5*y_max)/(sigma_y*sigma_y));
}
double rho2(double x, double y)
{
    double x_max=Delta*n_x;
    double y_max=Delta*n_y;
    double sigma_x=0.1*x_max;
    double sigma_y=0.1*y_max;

    return -exp(-(x-0.65*x_max)*(x-0.65*x_max)/(sigma_x*sigma_x)-(y-0.5*y_max)*(y-0.5*y_max)/(sigma_y*sigma_y));
}
double rho3(double x, double y)
{
    double x_max=Delta*n_x;
    double y_max=Delta*n_y;
    double sigma_x=0.1*x_max;
    double sigma_y=0.1*y_max;

    return -exp(-(x-0.45*x_max)*(x-0.45*x_max)/(sigma_x*sigma_x)-(y-0.75*y_max)*(y-0.75*y_max)/(sigma_y*sigma_y));
}

double rho(double x, double y)
{
    return 3.5*rho1(x,y)+2.5*rho2(x,y)+4.0*rho3(x,y);
}
void matrix_V_0(double **matrix)
{
    for (int i=0; i<=n_x; i++)
    {
        for (int j=0; j<=n_y; j++)
        {
            if (j==0)
            {
                matrix[i][j]=V1;
            }
            else if (j==n_y)
            {
                matrix[i][j]=V2;
            }
            else
            {
                matrix[i][j]=0.0;
            }
            //printf("%lf ", matrix[i][j]);
        }
        //printf("\n");
        
    }
}
double ITERATION_STOP(double **V)
{
    double S=0.0;
    for (int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_y; j++)
        {
            S+=Delta*Delta*(0.5*(V[i+1][j]-V[i][j])/Delta + 0.5*(V[i][j+1]-V[i][j])/Delta -rho(i*Delta, j*Delta)*V[i][j]);
        }
    }
    return S;
}
void relax_global(double **matrix_V_n, double **matrix_V_s, double *S, double omega)
{
    
    for (int i=1; i<n_x; i++)
    {
        for (int j=1; j<n_y; j++)
        {
            matrix_V_n[i][j]=0.25*(matrix_V_s[i+1][j]+matrix_V_s[i-1][j]+matrix_V_s[i][j-1]+matrix_V_s[i][j+1]+Delta*Delta*rho(i*Delta, j*Delta)/epislon);

        }
    }
    for (int j=0; j<=n_y; j++)
    {
        matrix_V_n[0][j]=matrix_V_n[1][j];
        matrix_V_n[n_x][j]=matrix_V_n[n_x-1][j];
    }
    for (int i=0; i<=n_x; i++)
    {
        matrix_V_n[i][0]=V1;
        matrix_V_n[i][n_y]=V2;
    }
    S[0]=ITERATION_STOP(matrix_V_s);
    S[1]=ITERATION_STOP(matrix_V_n);
    for(int i=0; i<=n_x; i++)
    {
        for(int j=0; j<=n_y; j++)
        {
            matrix_V_s[i][j]=(1-omega)*matrix_V_s[i][j]+omega*matrix_V_n[i][j];
        }
    }
}
void gradient(double **V, double **grad_Vx, double **grad_Vy)
{
    for (int i=1; i<n_x; i++)
    {
        for (int j=1; j<n_y; j++)
        {
            grad_Vx[i][j]=(V[i+1][j]-V[i-1][j])/(2.0*Delta);
            grad_Vy[i][j]=(V[i][j+1]-V[i][j-1])/(2.0*Delta);
        }
    }
    for (int i=0; i<=n_x; i++)
    {
        grad_Vx[i][0]=grad_Vx[i][1];
        grad_Vx[i][n_y]=grad_Vx[i][n_y-1];
        grad_Vy[i][0]=grad_Vy[i][1];
        grad_Vy[i][n_y]=grad_Vy[i][n_y-1];

    }
    for (int j=0; j<n_y; j++)
    {
        grad_Vx[0][j]=grad_Vx[1][j];
        grad_Vx[n_x][j]=grad_Vx[n_x-1][j];
        grad_Vy[0][j]=grad_Vy[1][j];
        grad_Vy[n_x][j]=grad_Vy[n_x-1][j];
    }
    grad_Vx[0][0]=grad_Vx[1][1];
    grad_Vx[n_x][0]=grad_Vx[n_x-1][1];
    grad_Vx[n_x][n_y]=grad_Vx[n_x-1][n_y-1];
    grad_Vx[0][n_y]=grad_Vx[1][n_y-1];

    grad_Vy[0][0]=grad_Vy[1][1];
    grad_Vy[n_x][0]=grad_Vy[n_x-1][1];
    grad_Vy[n_x][n_y]=grad_Vy[n_x-1][n_y-1];
    grad_Vy[0][n_y]=grad_Vy[1][n_y-1];

}
void ion_postion(double delta_t, double *r, double *v, double **grad_Vx, double **grad_Vy)
{
    FILE *position=fopen("position.txt", "w");
    
    double error;
    int IT_MAX=0;
    double t=0;
    double IT_NUMB=t_max/delta_t;
    double r_holder[2];
    double v_holder[2];

    double r_old[2];
    double v_old[2];

    bool STOP=false;

    int Vx,Vy,Vx_old,Vy_old;

    for (int i=0; i<2; i++)
    {
        r_holder[i]=r[i];
        v_holder[i]=v[i];
        r_old[i]=r[i];
        v_old[i]=v[i];
    }

    for (int i=0; i<=IT_NUMB; i++)
    {
        fprintf(position, "%lf\t%lf\t%lf\n", t, r_old[0], r_old[1]);
        t=i*delta_t;
        do
        { 
            
            IT_MAX++;
            Vx_old = (int)round(r_old[0] / Delta);
            Vy_old = (int)round(r_old[1] / Delta);
            Vx = (int)round(r_holder[0] / Delta);
            Vy = (int)round(r_holder[1] / Delta);
            if (Vx>n_x||Vy>n_y||Vx<0||Vy<0)
            {
                STOP=true;
            }
            if (STOP==false)
            {
                r[0]=r_old[0]+delta_t*(v_old[0]+v_holder[0])/2.0;
                r[1]=r_old[1]+delta_t*(v_old[1]+v_holder[1])/2.0;

                v[0]=v_old[0]-delta_t*(grad_Vx[Vx_old][Vy_old]+grad_Vx[Vx][Vy])/(2.0*m);
                v[1]=v_old[1]-delta_t*(grad_Vy[Vx_old][Vy_old]+grad_Vy[Vx][Vy])/(2.0*m);
            error=(fabs(r[0]-r_holder[0])+fabs(r[1]-r_holder[1])+fabs(v[0]-v_holder[0])+fabs(v[1]-v_holder[1]))/4.0;
            }
            for (int j=0; j<2; j++)
            {
                r_holder[j]=r[j];
                v_holder[j]=v[j];
            }
        }while(error>=TOL&&STOP==false&&IT_MAX<100);
        for (int j=0; j<2; j++)
        {
            r_old[j]=r[j];
            v_old[j]=v[j];
        }
        IT_MAX=0;

    }
    fclose(position);
}
int main()
{
    FILE *mop=fopen("map_of_potential.txt","w");
    
    double **matrix_V_n = malloc((n_x+1) * sizeof(double *));
    for (int j = 0; j < n_x+1; j++) {
        matrix_V_n[j] = malloc((n_y+1) * sizeof(double));
    }
    double **matrix_V_s = malloc((n_x+1) * sizeof(double *));
    for (int j = 0; j < n_x+1; j++) {
        matrix_V_s[j] = malloc((n_y+1) * sizeof(double));
    }
    double **matrix_grad_Vx = malloc((n_x+1) * sizeof(double *));
    for (int j = 0; j < n_x+1; j++) {
        matrix_grad_Vx[j] = malloc((n_y+1) * sizeof(double));
    }
    double **matrix_grad_Vy = malloc((n_x+1) * sizeof(double *));
    for (int j = 0; j < n_x+1; j++) {
        matrix_grad_Vy[j] = malloc((n_y+1) * sizeof(double));
    }
    
    double *S=malloc(2*sizeof(double));
    double *r=malloc(2*sizeof(double));
    double *v=malloc(2*sizeof(double));
    r[0]=0.45*n_x*Delta;
    r[1]=0.5*n_y*Delta;
    v[0]=-0.003*n_x*Delta;
    v[1]=0.01*n_y*Delta;
    double omega=0.5;
    double error;
    double S_new, S_old;
    int IT_NUMB=0;
    
    matrix_V_0(matrix_V_s);
    printf("POTENTIAL_CALC\n");
    do
    {
        IT_NUMB++;
        relax_global(matrix_V_n, matrix_V_s, S, omega);
        error=fabs((S[1]-S[0])/(S[0]+1e-16));
        printf("%d\n", IT_NUMB);
    }while(error>TOL&&IT_NUMB<1000);
    
    
    for (int i=0; i<=n_x; i++)
    {
        for(int j=0; j<=n_y; j++)
        {
            fprintf(mop, "%lf\t%lf\t%lf\n", Delta*i, Delta*j, matrix_V_s[i][j]);
        }
        fprintf(mop, "\n");
    }
    printf("GRAD_CALC\n");
    gradient(matrix_V_n,matrix_grad_Vx, matrix_grad_Vy);
    printf("POSITION_CALC\n");
    ion_postion(0.01, r, v, matrix_grad_Vx, matrix_grad_Vy);
    
    
    for (int j = 0; j < n_x+1; j++) {
        free(matrix_V_n[j]);
        free(matrix_V_s[j]);
        free(matrix_grad_Vx[j]);
        free(matrix_grad_Vy[j]);
    }
    free(matrix_V_n);
    free(matrix_V_s);
    free(matrix_grad_Vx);
    free(matrix_grad_Vy);
    free(S);
    free(r);
    free(v);
    fclose(mop);
    return 0;
}