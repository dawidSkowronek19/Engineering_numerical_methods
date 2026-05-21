#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

const int nx=150;
const double Delta=0.1;
const double dt=0.05;
const double alpha=1.0;
const double beta=1.0;
const double xF=2.0;
const double t_max=50.0;
const double sigma=0.5;
const double xA=7.5;

double aF(double t, double x)
{
    if (fabs(x-xA)<(Delta/2.0))
    {
        return cos(50.0*t/t_max);
    }
    else
    {
        return 0;
    }
}
void boundaryConditions(std::vector<double> &u, std::vector<double> &v)
{
    u.at(0)=0.0;
    u.at(nx)=0.0;
    v.at(0)=0.0;
    v.at(nx)=0.0;
}
void startConditions(std::vector<double> &u, std::vector <double> &v)
{
    for (int i=0; i<=nx; i++)
    {
        u.at(i)=exp(-(i*Delta-xA)*(i*Delta-xA)/(2.0*sigma*sigma));
        //u.at(i)=sin(i*Delta)*(i*Delta);
        v.at(i)=0.0;
    }

}

void verlet(std::vector<double> &u0, std::vector<double> &u, std::vector<double> &vHolder, std::vector<double> &v, double t)
{
    std::vector<double> a;
    std::vector<double> utmp;
    a.resize(nx+1);
    utmp.resize(nx+1);
    utmp=u;
    for (int i=1; i<nx; i++) //aN
    {
        a.at(i)=(u.at(i+1)-2*u.at(i)+u.at(i-1))/(Delta*Delta)-beta*(u.at(i)-u0.at(i))/dt +alpha*aF(t, Delta*i);
    }
    for (int i=1; i<nx; i++)// v(N+0.5) u(N+1)
    {
        vHolder.at(i)=v.at(i)+dt*a.at(i)/2.0;
        u.at(i)=u0.at(i)+dt*vHolder.at(i);
    }
    u0=utmp;
    for (int i=1; i<nx; i++)// a(N+1)
    {
        a.at(i)=(u.at(i+1)-2*u.at(i)+u.at(i-1))/(Delta*Delta)-beta*(u.at(i)-u0.at(i))/dt +alpha*aF(t, Delta*i);
    }
        for (int i=0; i<nx; i++)// v(N+0.5) u(N+1)
    {
        v.at(i)=vHolder.at(i)+dt*a.at(i)/2.0;
    }
    u0=u;

}

double Energy(double t, std::vector<double> u, std::vector<double> v)
{
    double Sum=Delta*(((u.at(1)-u.at(0))/Delta)*(u.at(1)-u.at(0))/Delta+(u.at(nx)-u.at(nx-1))/Delta*(u.at(nx)-u.at(nx-1))/Delta)/4.0;
    for (int i=1; i<nx; i++)
    {
        Sum+=Delta*(v.at(i)*v.at(i)+(u.at(i+1)-u.at(i-1))/(2.0*Delta)*(u.at(i+1)-u.at(i-1))/(2.0*Delta))/2.0;
    }
    return Sum;
}
int main()
{
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> u0;
    std::vector<double> vHolder;
    double t=0.0;
    u.resize(nx+1);
    v.resize(nx+1);
    u0.resize(nx+1);
    vHolder.resize(nx+1);

    startConditions(u,v);
    boundaryConditions(u,v);
    u0=u;
    vHolder=v;
    double E;
    //std::ofstream fileU("./data/u.dat");
    std::ofstream fileE("./data4/E.dat");
    while(t<t_max)
    {
        t+=dt;
        std::ofstream fileU("./data4/u.dat"+std::to_string(t));
        verlet(u0,u, vHolder, v, t);
        E=Energy(t,u,v);
        for (int i=0; i<=nx; i++)
        {
            fileU << i*Delta <<" "<< u.at(i)<<"\n";
        }
        fileE << t <<" "<< E<<"\n";
        

    }
    return 0;
}