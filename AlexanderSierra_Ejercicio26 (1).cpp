#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

void Metodo_Explicito(int x_init,int x_end,double delta_t,double omega,string filename);
void Metodo_Explicito(int x_init,int x_fin,double delta_t,double omega,string filename)
    {
    double y = 1,x = x_init, z_0 = 0,z_n = 0;
    ofstream outfile;
    outfile.open(filename);
    while(x<x_fin)
        {
        z_n=z_0;
        z_0=z_0-delta_t*omega*omega*y;
        y=y+delta_t*z_n;
        x=x+delta_t;
        outfile<<x<<" "<<y<<" "<<z_0<<endl;
        }
    outfile.close();
    }

void Runge_kutta4th(int x_init,int x_end,double delta_t,double omega,string filename);
void Runge_kutta4th(int x_init,int x_end,double delta_t,double omega,string filename)
    {
    ofstream outfile;
    outfile.open(filename);
    
    float y=1,x=x_init,z_0=0,z_n=0;
    float y_next = 0,z_next = 0;
    float f_z0=0,f_y0=0;
    float f_z1=0,f_y1=0;
    float f_z2=0,f_y2=0;
    float f_z3=0,f_y3=0;
    float f_zmean=0,f_ymean=0;

    while(x<x_end)
        {
        z_n=z_0;
        f_z0=-delta_t*omega*omega*y;
        f_y0=delta_t*z_n;
        y_next=y+f_y0*1/2;
        z_next=z_0+f_z0*1/2;
       
        f_z1=-delta_t*omega*omega*y_next;
        f_y1=delta_t*z_next;
        y_next=y+f_y1*1/2;
        z_next=z_0+f_z1*1/2;
           
        f_z2=-delta_t*omega*omega*y_next;
        f_y2=delta_t*z_next;
        y_next=y+f_y2;
        z_next=z_0+f_z2;
       
        f_z3=-delta_t*omega*omega*y_next;
        f_y3=delta_t*z_next;
        
        f_zmean=f_z0*1/6+f_z1*1/3+f_z2*1/3+f_z3*1/6;
        f_ymean=f_y0*1/6+f_y1*1/3+f_y2*1/3+f_y3*1/6;
        z_0=z_0+f_zmean;
        y=y+f_ymean;
        x=x+delta_t;
        outfile<<x<<" "<<y<<" "<<z_0<<endl;
        }
    outfile.close();
    }


int main(){
   
    float omega = 1.0;
    Metodo_Explicito(0.0, 10000.0, omega/2, omega, "euler.dat");
    Runge_kutta4th(0.0, 10000.0, omega/2, omega, "rk4.dat");
    //leap_frog(0.0, 10000.0, omega/2, omega, "leapfrog.dat");
   
    return 0;  
}


