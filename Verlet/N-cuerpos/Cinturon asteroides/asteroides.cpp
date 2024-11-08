// asteroides.cpp
// Aplicacion del metodo de velocidad de Verlet 
// para el movimiento gravitacional.

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <fstream>
#include <chrono>
#include <random>

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
const double masa_asteroide = 1.3877e16; // 1.3877e21/100e3.

// Inicializar masas.
void init_masa()
{
    masa[0] = 1.9890e30; // Sol.
    masa[1] = 3.3020e23; // Mercurio.
    masa[2] = 4.8700e24; // Venus.
    masa[3] = 5.9720e24; // Tierra.
    masa[4] = 6.4185e23; // Marte.
    masa[5] = 1.8990e27; // Jupiter.
    // Asteroides.
    masa[6]  = 9.4300e20; // Ceres.
    masa[7]  = 2.2000e20; // Palas.
    masa[8]  = 2.7100e20; // Vesta.
    masa[9]  = 9.0300e17; // Higia.
    masa[10] = 3.0000e19; // Juno. 
    
    double masa_min = 0.01*masa_asteroide;
    double masa_max = masa_asteroide;
    
    srand(time(0));
    
    // Asteroides pequenios.
    for (int i = 11; i < n_cuerpos; i++) 
        masa[i] = masa_min + (rand()/RAND_MAX)*(masa_max - masa_min);    

}

// Inicializar posiciones.
void init_posicion()
{
    xp[0]  = 0.0;       // Sol.
    yp[0]  = 0.0; 
    xp[1]  = 0.3870*ua; // Mercurio.
    yp[1]  = 0.0;
    xp[2]  = 0.7233*ua; // Venus.
    yp[2]  = 0.0;
    xp[3]  = ua;        // Tierra.  
    yp[3]  = 0.0;
    xp[4]  = 1.5237*ua; // Marte.
    yp[4]  = 0.0; 
    xp[5]  = 5.2034*ua; // Jupiter.
    yp[5]  = 0.0;
    
    // Asteorides grandes.
    xp[6]  = 2.7664*ua; // Ceres. 
    yp[6]  = 0.0; 
    xp[7]  = 2.7720*ua; // Palas. 
    yp[7]  = 0.0;
    xp[8]  = 2.3615*ua; // Vesta. 
    yp[8]  = 0.0;
    xp[9]  = 3.1395*ua; // Higia. 
    yp[9]  = 0.0;
    xp[10] = 2.6705*ua; // Juno.
    yp[10] = 0.0;
    
    // Asteorides pequenios. 
    // Valores aleatorios entre [-3.65, 3.65] UA en un cuadrado xy.
    for (int i = 11; i < n_cuerpos; i++)
    {
        while (true)
        {
            double xi = random_interval(generator) * (3.65*ua);
	    double yi = random_interval(generator) * (3.65*ua);

            double r = sqrt(pow(xi, 2) + pow(yi, 2));

            if (r >= 2.06*ua && r <= 3.65*ua)
            {
                xp[i] = xi;
                yp[i] = yi;

                break;
            }   
        }    
    }
}

// Inicializar velocidades.
void init_velocidad()
{
    vx[0]  = 0.0;      // Sol.
    vy[0]  = 0.0;
    vx[1]  = 0.0;      // Mercurio.
    vy[1]  = 4.7872e4;
    vx[2]  = 0.0;      // Venus.
    vy[2]  = 3.5000e4;
    vx[3]  = 0.0;      // Tierra.
    vy[3]  = 2.9780e4;
    vx[4]  = 0.0;      // Marte.
    vy[4]  = 2.4077e4;
    vx[5]  = 0.0;      // Jupiter.
    vy[5]  = 1.3070e4;
    
    // Asteorides grandes.
    vx[6]  = 0.0;      // Ceres. 
    vy[6]  = 1.7882e4;
    vx[7]  = 0.0;      // Palas. 
    vy[7]  = 1.7650e4;
    vx[8]  = 0.0;      // Vesta. 
    vy[8]  = 1.9340e4;
    vx[9]  = 0.0;      // Higia. 
    vy[9]  = 1.6760e4;
    vx[10] = 0.0;      // Juno. 
    vy[10] = 1.7930e4;
    
    // Asteroides pequenios.
    for (int i = 11; i < n_cuerpos; i++)
    {
        double r = sqrt(pow(xp[i], 2) + pow(yp[i], 2));
        vx[i] = sqrt(G*masa[0]/r)*(-yp[i]/r);
        vy[i] = sqrt(G*masa[0]/r)*(xp[i]/r);
    }    
}

void energia_momentos()
{
    E = 0.0; 
    P = 0.0; 
    L = 0.0;
    
    double K = 0.0;
    double U = 0.0;
    
    for (int i = 0; i < n_cuerpos; i++) 
    {
        double r = sqrt(xp[i]*xp[i] + yp[i]*yp[i]);
        double v = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
        
        K += 0.5*masa[i]*v*v;
        P += masa[i]*v;
        L += masa[i]*r*v;
        
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
    for (int i = 0; i < n_cuerpos; i++)
    {
        xa[i] = 0;
        ya[i] = 0;

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
    const double h_step  = 8.64e3; // Tamanio de paso (1 dia).
    const long int Niter = 1.825e4; // Numero de iteraciones (50 anios).  
    const int out_cada   = 1;   // Out_cada iteraciones.
    
    n_cuerpos = 100; // Numero de cuerpos.

    // Otras variables.
    const int n_ec = n_cuerpos * 4; // Numero ecuaciones.
    double tiempo = 0.0;

    // Archivos de salida.
    ofstream of_posicion("ast-posicion.dat", ios::out);
    // ofstream of_velocidad("velocidad.dat", ios::out);
    ofstream of_energia("ast-energia.dat", ios::out);
    ofstream of_momentumLineal("ast-momLineal.dat", ios::out);
    ofstream of_momentumAngular("ast-momAngular.dat", ios::out);
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


