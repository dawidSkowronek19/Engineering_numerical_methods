#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double PI =M_PI;
#define M 1.0
#define t_max 2.0
#define TOL 1e-6
#define S 0.75
#define x_0 0.0
#define y_0 0.0
const double G = 4.0*PI*PI;

typedef struct position{double x,y;}position;

position Force(position r)
{
    position a;
    double radius=sqrt((r.x-x_0)*(r.x-x_0)+(r.y-y_0)*(r.y-y_0));
    a.x=-G*M*(r.x-x_0)/(radius*radius*radius);
    a.y=-G*M*(r.y-y_0)/(radius*radius*radius);

    return a;
}

void RK4(position *r_f, position *v_f, double dt)
{
    position r_tmp, r,v;
    position a;
    r=r_f[0];
    v=v_f[0];

    position k_r[4];
    position k_v[4];
    
    r_tmp.x=r.x;
    r_tmp.y=r.y;
    a=Force(r_tmp);

    k_r[0].x=v.x;
    k_r[0].y=v.y;
    k_v[0].x=a.x;
    k_v[0].y=a.y;

    r_tmp.x=r.x+dt*k_r[0].x/2.0;
    r_tmp.y=r.y+dt*k_r[0].y/2.0;
    a=Force(r_tmp);

    k_r[1].x=v.x+dt*k_v[0].x/2.0;
    k_r[1].y=v.y+dt*k_v[0].y/2.0;
    k_v[1].x=a.x;
    k_v[1].y=a.y;

    r_tmp.x=r.x+dt*k_r[1].x/2.0;
    r_tmp.y=r.y+dt*k_r[1].y/2.0;
    a=Force(r_tmp);

    k_r[2].x=v.x+dt*k_v[1].x/2.0;
    k_r[2].y=v.y+dt*k_v[1].y/2.0;
    k_v[2].x=a.x;
    k_v[2].y=a.y;

    r_tmp.x=r.x+dt*k_r[2].x;
    r_tmp.y=r.y+dt*k_r[2].y;
    a=Force(r_tmp);

    k_r[3].x=v.x+dt*k_v[2].x;
    k_r[3].y=v.y+dt*k_v[2].y;
    k_v[3].x=a.x;
    k_v[3].y=a.y;

    r_f[0].x+=dt*(k_r[0].x+2.0*k_r[1].x+2.0*k_r[2].x+k_r[3].x)/6.0;
    r_f[0].y+=dt*(k_r[0].y+2.0*k_r[1].y+2.0*k_r[2].y+k_r[3].y)/6.0;
    v_f[0].x+=dt*(k_v[0].x+2.0*k_v[1].x+2.0*k_v[2].x+k_v[3].x)/6.0;
    v_f[0].y+=dt*(k_v[0].y+2.0*k_v[1].y+2.0*k_v[2].y+k_v[3].y)/6.0;

}
double max_value(position r_d, position v_d, position r_o, position v_o)
{
    double p=4.0;
    double denom=pow(2,p)-1.0;
    double E[4];
    E[0]=(r_d.x-r_o.x)/denom;
    E[1]=(r_d.y-r_o.y)/denom;
    E[2]=(v_d.x-v_o.x)/denom;
    E[3]=(v_d.y-v_o.y)/denom;
    double max=fabs(E[0]);
    for (int i=1; i<4; i++)
    {
        if (fabs(E[i])>max)
        {
            max=fabs(E[i]);
        }
    }
    return max;
}
double step_control(double max,double dt_old)
{
    double p=4.0;
    double dt_new;
    dt_new=pow(S*TOL/max,1.0/(p+1))*dt_old;
    return dt_new;
    
}
double Energy(position r, position v)
{
    return (v.x*v.x+v.y*v.y)/2.0 - G*M/sqrt(r.x*r.x+r.y*r.y);
}
double momentum(position r, position v)
{
    return fabs(r.x*v.y-r.y*v.x);
}
int main()
{
    position r_RK4_d, v_RK4_d, r_RK4_o, v_RK4_o, r_holder, v_holder;
    double dt=0.01;
    double t=0.0;
    double max;
    r_holder.x= 1.0;
    r_holder.y=0.0;
    v_holder.x=0.0;
    v_holder.y=sqrt(G);
    r_RK4_d=r_holder;
    v_RK4_d=v_holder;
    r_RK4_o=r_holder;
    v_RK4_o=v_holder;

    FILE *out_RK4=fopen("out_RK4.dat", "w");
    fprintf(out_RK4, "%e\t%e\t%e\t%e\t%e\t%e\n", t, r_RK4_d.x, r_RK4_d.y, dt,Energy(r_RK4_d, v_RK4_d), momentum(r_RK4_d, v_RK4_d));
    do
    {
        r_RK4_d=r_holder;
        v_RK4_d=v_holder;
        r_RK4_o=r_holder;
        v_RK4_o=v_holder;
        
        RK4(&r_RK4_d, &v_RK4_d, dt);
        RK4(&r_RK4_d, &v_RK4_d, dt);
        RK4(&r_RK4_o, &v_RK4_o, 2.0*dt);
        max=max_value(r_RK4_d, v_RK4_d, r_RK4_o, v_RK4_o);
        if (max<TOL)
        {
            t+=2*dt;
            fprintf(out_RK4, "%e\t%e\t%e\t%.12e\t%.12e\t%.12e\n", t, r_RK4_d.x, r_RK4_d.y, dt, Energy(r_RK4_d, v_RK4_d), momentum(r_RK4_d, v_RK4_d)); // t, r.x, r.y, dt, E, L
            r_holder.x= r_RK4_d.x;
            r_holder.y=r_RK4_d.y;
            v_holder.x=v_RK4_d.x;
            v_holder.y=v_RK4_d.y;
        }
        dt=step_control(max, dt);
    }while(t<t_max);

    fclose(out_RK4);
    return 0;
}