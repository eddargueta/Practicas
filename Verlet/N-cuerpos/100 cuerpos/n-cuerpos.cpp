// n-cuerpos.cpp
// Aplicacion del metodo de velocidad de Verlet 
// para el movimiento gravitacional.
// Para correr con varios nucleos y optimizacion (-O2):
// 1. Compilar con g++ -O2 -fopenmp -o n-cuerpos.out n-cuerpos.cpp.
// 2. Escribir en la terminal OMP_NUM_THREADS=X (X: numero de nucleos).
// 3. ./n-cuerpos.out.

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <fstream>
#include <chrono>
#include <random>
#include <omp.h>

using namespace std;

random_device rd;
default_random_engine generator(rd()); // rd() random seed.
uniform_real_distribution<long double> random_interval(-1.0, 1.0);

// Variables globales.
double *xp, *yp, *xpn, *ypn; // Posicion en x, y.
double *vx, *vy, *vxn, *vyn; // Velocidad en x, y.
double *ax, *ay, *axn, *ayn; // Aceleracion en x, y.
double *masa; // Masa de cada particula.
double E, P, L; 
int n_cuerpos; // Numero de cuerpos.

// Constantes.
const double G  = 6.6743e-11;
const double ua = 1.4960e11; //1.50e11;

// Inicializar masas.
void init_masa()
{
    /*masa[0] = 7.349e22; // Luna.
    masa[1] = 5.972e24; // Tierra.
    masa[2] = 1.989e30; // Sol.*/
    
    //Las masas inician con el mismo valor.
    for (int i = 0; i < n_cuerpos; i++)
        masa[i] = 10e18;
}

// Inicializar posiciones.
void init_posicion()
{
    /*xp[0] = 3.8440e8 + ua; //Luna. 
    yp[0] = 0.0;
    xp[1] = ua; // Tierra.  
    yp[1] = 0.0;
    xp[2] = 0.0; // Sol.
    yp[2] = 0.0;*/ 
    
    // Valores aleatorios entre [-1, 1] UA en un cuadrado xy.
    for (int i = 0; i < n_cuerpos; i++)
    {
        xp[i] = random_interval(generator) * ua;
        yp[i] = random_interval(generator) * ua;
    }
}

// Inicializar velocidades.
void init_velocidad()
{
    /*vx[0] = 0.0; // Luna.
    vy[0] = 3.0802e4;
    vx[1] = 0.0; // Tierra.
    vy[1] = 2.9780e4;
    vx[2] = 0.0; // Sol.
    vy[2] = 0.0;*/
       
    for (int i = 0; i < n_cuerpos; i++)
    {
        double r = sqrt(xp[i]*xp[i] + yp[i]*yp[i]);
        double v = sqrt(G * M_PI * n_cuerpos * masa[i] * r)/(2*ua);

        vx[i] = v * (-yp[i]/r + random_interval(generator));
        vy[i] = v * (xp[i]/r + random_interval(generator));
    } 
}

void energia_momentos()
{
    E = 0.0; 
    P = 0.0; 
    L = 0.0;
  
    double K = 0.0;
    double U = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < n_cuerpos; i++) 
    {
        double r = sqrt(xp[i]*xp[i] + yp[i]*yp[i]);
        double v = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
        
        K += 0.5*masa[i]*v*v;
        P += masa[i]*v;
        L += masa[i]*r*v;
        
        #pragma omp parallel for
        for (int j = 0; j < i; j++)
        {
            double dx  = xp[i] - xp[j];
            double dy  = yp[i] - yp[j];
            double rij = sqrt(dx*dx + dy*dy); 	
            
            U += G*masa[i]*masa[j]/rij;
        }
    }
    
    E = K - U;
}

void aceleracion(double *px, double *py, double *xa, double *ya, 
    const double t, double *ac)
{
    #pragma omp parallel for
    for (int i = 0; i < n_cuerpos; i++)
    {
        xa[i] = 0;
        ya[i] = 0;

        #pragma omp parallel for
        for (int j = 0; j < n_cuerpos; j++)
        {
            if (i == j)
                continue;
            
            double dx = px[i] - px[j];
            double dy = py[i] - py[j];
            double r3 = pow(dx*dx + dy*dy, 1.5);
            
            xa[i] -= (G*masa[j]*dx)/r3;
            ya[i] -= (G*masa[j]*dy)/r3;
        }   
    }
  
    // Llena los arreglos de a(t)/a(t+dt).
    for(int i = 0; i < n_cuerpos; i++)
    {
        ac[i            ] = xa[i];
        ac[i + n_cuerpos] = ya[i];
    }
}

void act_posicion(double *y_nueva, const double t, const double h)
{
    // Calcula posicion r(t+dt).
    for (int i = 0; i < n_cuerpos; i++)
    {
        xpn[i] = xp[i] + vx[i]*h + 0.5*ax[i]*h*h;
        ypn[i] = yp[i] + vy[i]*h + 0.5*ay[i]*h*h;
    }
    
    // Llena arrgeglo r(t+dt).
    for (int i = 0; i < n_cuerpos; i++)
    {
    	y_nueva[i            ] = xpn[i];
        y_nueva[i + n_cuerpos] = ypn[i];
    }  
}

void act_velocidad(double *y_nueva, const double t, const double h)
{
    // Calcula velocidad v(t+dt).
    for (int i = 0; i < n_cuerpos; i++)
    {
        vxn[i] = vx[i] + 0.5*(ax[i] + axn[i])*h;
        vyn[i] = vy[i] + 0.5*(ay[i] + ayn[i])*h;
    }
    
    // Llena arreglo v(t+dt).
    for (int i = 0; i < n_cuerpos; i++)
    {
    	y_nueva[i + 2*n_cuerpos] = vxn[i];
        y_nueva[i + 3*n_cuerpos] = vyn[i];
    }          
}

void verificar_colisiones(double t)
{
    double dist_lim = 1.0e7;
    
    #pragma omp parallel for
    for (int i = 0; i < n_cuerpos; i++)
    {
        if (masa[i] != 0.0)
        {
            #pragma omp parallel for
            for (int j = 0; j < i; j++)
            {
                if (masa[j] != 0.0)
                {
                    double dx = xp[i] - xp[j];
                    double dy = yp[i] - yp[j];
                    double distancia = sqrt(pow(dx, 2) + pow(dy, 2));

                    if (distancia < dist_lim)
                    {
                        double nueva_masa = masa[i] + masa[j];
                        vx[i] = (masa[i] * vx[i] + masa[j] * vx[j]) / nueva_masa;
                        vy[i] = (masa[i] * vy[i] + masa[j] * vy[j]) / nueva_masa;
                        masa[i] = nueva_masa;

                        // Particula j sigue la misma trayectoria que particula i 
                        // pero sin masa.
                        xp[j] = (xp[i] + xp[j]) / 2;
                        yp[j] = (yp[i] + yp[j]) / 2;
                        vx[j] = vx[i];
                        vy[j] = vy[i];
                        ax[j] = ax[i];
                        ay[j] = ay[i];
                        masa[j] = 0.0;
                        cout << "Colision " << i << " " << j << " en t = " << t << endl;
                    }   
                }   
            }   
        }        
    } 
}

void salidaSolucion(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    
    for (int i = 0; i < n_cuerpos; i++)
        of << scientific << setprecision(9) << "\t" << xp[i];
        
    for (int i = 0; i < n_cuerpos; i++)
        of << scientific << setprecision(9) << "\t" << yp[i];
    
    of << endl;
}

/*void salidaVelocidad(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    
    for (int i = 0; i < n_cuerpos; i++)
        of << scientific << setprecision(9) << "\t" << vx[i];
        
    for (int i = 0; i < n_cuerpos; i++)
        of << scientific << setprecision(9) << "\t" << vy[i];
    
    of << endl;
}*/

void salidaEnergia(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    of << scientific << setprecision(9) << "\t" << E << endl;
}

void salidaMomentumLineal(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    of << scientific << setprecision(9) << "\t" << P << endl;
}

void salidaMomentumAngular(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    of << scientific << setprecision(9) << "\t" << L << endl;
}

/*void salidaPCM(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    of << scientific << setprecision(9) << "\t" << Pcm << endl;
}

void salidaMasa(double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;

    for (int i = 0; i < n_cuerpos; i++)
        of << "\t" << masa[i];

    of << endl;
}*/

int main()
{
    // Parametros.
    const double h_step  = 24*3600; // Tamanio de paso (1 dia).
    const long int Niter = 1.825e6; //1825000; // Numero de iteraciones (5000 anios).  
    const int out_cada   = 1000;    // Out_cada iteraciones.
    
    n_cuerpos = 100; // Numero de cuerpos.

    // Otras variables.
    const int n_ec = n_cuerpos * 4; // Numero ecuaciones.
    double tiempo = 0.0;

    // Archivos de salida.
    ofstream of_posicion("100-posicion.dat", ios::out);
    // ofstream of_velocidad("velocidad.dat", ios::out);
    ofstream of_energia("100-energia.dat", ios::out);
    ofstream of_momentumLineal("100-momLineal.dat", ios::out);
    ofstream of_momentumAngular("100-momAngular.dat", ios::out);
    //ofstream of_masa("masa.dat", ios::out);
    //ofstream of_centroMasa("centroMasa.dat", ios::out);

    // Reserva espacio para
    double *y	     = new double[n_ec]; // Posiciones y velocidades r(t), v(t),
    double *y_nueva  = new double[n_ec]; // r(t+dt), v(t+dt).
    double *a        = new double[2*n_cuerpos]; // Aceleraciones a(t),
    double *a_nueva  = new double[2*n_cuerpos]; // a(t+dt).
    
    // Reserva espacio para posicion, velocidad, 
    // aceleracion y masa.
    xp = y;
    yp = y +   n_cuerpos;
    vx = y + 2*n_cuerpos;
    vy = y + 3*n_cuerpos;
    ax = a;
    ay = a + n_cuerpos;
    xpn = y_nueva;
    ypn = y_nueva +   n_cuerpos;
    vxn = y_nueva + 2*n_cuerpos;
    vyn = y_nueva + 3*n_cuerpos;
    axn = a_nueva;
    ayn = a_nueva + n_cuerpos;
    masa = new double[n_cuerpos];
  
    for (int i = 0; i < n_ec; i++)
        y_nueva[i] = 0.0;

    // Inicializa masas.
    init_masa();

    // Inicializa posiciones.
    init_posicion();

    // Inicializa velocidades. 
    init_velocidad();
    
    // Inicializa aceleraciones.
    aceleracion(xp, yp, ax, ay, tiempo, a);
    
    energia_momentos();

    verificar_colisiones(tiempo);

    // Ciclo principal.
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Niter; i++)
    {
        if (i%1000 == 0)
        {
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            cout << "iteracion: " << i << ", tardo en calcularse: " << duration.count() << "ms" << endl;
            start = stop;
        }
        
        xp = y;
        yp = y +  n_cuerpos;
        vx = y + 2*n_cuerpos;
        vy = y + 3*n_cuerpos;
        ax = a;
        ay = a + n_cuerpos;
        xpn = y_nueva;
        ypn = y_nueva +   n_cuerpos;
        vxn = y_nueva + 2*n_cuerpos;
        vyn = y_nueva + 3*n_cuerpos;
        axn = a_nueva;
        ayn = a_nueva + n_cuerpos;
    
    
        if (i % out_cada == 0)            
        {
            salidaSolucion(tiempo, of_posicion);
            salidaEnergia(tiempo, of_energia);
            //salidaVelocidad(tiempo, of_velocidad);
            salidaMomentumLineal(tiempo, of_momentumLineal);
            salidaMomentumAngular(tiempo, of_momentumAngular);
            //salidaPCM(tiempo, of_centroMasa);
            //salidaMasa(tiempo, of_masa);
        }
        
        // Algoritmo de velocidad de Verlet:
        act_posicion(y_nueva, tiempo + h_step, h_step); // r(t+dt).
        aceleracion(xpn, ypn, axn, ayn, tiempo + h_step, a_nueva); // a(t+dt).
        act_velocidad(y_nueva, tiempo + h_step, h_step); // v(t+dt).
            
        energia_momentos();
        verificar_colisiones(tiempo);

        for (int i = 0; i < n_ec; i++)
        {
            y[i] = y_nueva[i];
            a[i] = a_nueva[i];
        }

        tiempo += h_step; // ti+1.    
    }

    // Libera el espacio de memoria reservado.
    delete[] y, y_nueva, a, a_nueva; 

    return 0;
}


