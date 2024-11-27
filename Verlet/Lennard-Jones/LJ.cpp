// LJ.cpp
// Aplicacion del metodo de velocidad de Verlet 
// para el potencial de Lennard-Jones en un gas
// monoatomico (argon).
// Para correr con varios nucleos y optimizacion (-O2):
// 1. Compilar con g++ -O2 -fopenmp -o LJ.out LJ.cpp.
// 2. Escribir en la terminal OMP_NUM_THREADS=X (X: numero de nucleos).
// 3. ./LJ.out.

/*
Temperaturas de referencia T0-->T*:
T*1 = 0°C = 273.15K = 2.28,
T*2 = 25°C = 298.15K = 2.48,
T*3 = 76.85°C = 350K = 2.92,
T*4 = 326.85°C = 600K = 5.00.
*/

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
double E, K, U; // Energias.
// double P, L, Pcm; // Energia y momenta del sistema.
int n_cuerpos; // Numero de cuerpos.

// Constantes.
const double l  = 45.0;//15.0; // Lado del cuadrado (unidades sigma).
const double rc = 2.5;  // 2.5*sigma radio de corte.
const double T0 = 0.1;  // Temperatura de referencia.
const double kB = 1.0;  // Cte. de Boltzmann. 

// Inicializar masas.
void init_masa()
{
    // Las masas inician con el mismo valor.
    for (int i = 0; i < n_cuerpos; i++)
        masa[i] = 1.0; // Masa unitaria.
}

// Inicializar posiciones.
void init_posicion()
{
    srand(time(0));

    // Valores aleatorios entre [-1, 1]l en un cuadrado xy.
    for (int i = 0; i < n_cuerpos; i++)
    {
        xp[i] = ((double)rand()/RAND_MAX)*l;
        yp[i] = ((double)rand()/RAND_MAX)*l;
        //xp[i] = random_interval(generator) * l;
        //yp[i] = random_interval(generator) * l;
    }    
}

// Inicializar velocidades.
void init_velocidad()
{
    K = 0.0;
    
    // Velocidades medias.
    double vxm = 0.0; 
    double vym = 0.0;
    double Km  = 0.0; //se puede borrar.
    
    for (int i = 0; i < n_cuerpos; i++)
    {
        vx[i] = random_interval(generator) * 0.5;
        vy[i] = random_interval(generator) * 0.5;
        
        // Energia cinetica inicial.
        K += 0.5*masa[i]*(vx[i]*vx[i] + vy[i]*vy[i]);
         
        vxm += vx[i];
        vym += vy[i];  
    }
    
    // Implementar fs para 1 < T* < 5:
    /*Km = K/n_cuerpos; 
    
    double T = Km/kB; 

    double fs = sqrt(T0/T); 
   
    vxm = vxm/n_cuerpos;
    vym = vym/n_cuerpos;*/
    
    // Asegura que la velocidad del centro de masa sea 0.
    for (int i = 0; i < n_cuerpos; i++)
    {
    	vx[i] = vx[i] - vxm;
        vy[i] = vy[i] - vym;
        
        //vx[i] = (vx[i] - vxm)*fs; // Para 1 < T* < 5.
        //vy[i] = (vy[i] - vym)*fs;  
    }    
}

void aceleracion(double *px, double *py, double *xa, double *ya, 
    const double t, double *ac)
{
    U = 0.0;
    
    double rc2 = rc*rc;
    double Uc = 4*((1/pow(rc2, 6)) - (1/pow(rc2, 3)));
    
    for (int i = 0; i < n_cuerpos; i++)
    {
        xa[i] = 0.0;
        ya[i] = 0.0;
    }
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n_cuerpos-1; i++)
    {
    	for (int j = i+1; j < n_cuerpos; j++)
    	{
	        double dx = px[i] - px[j];
            double dy = py[i] - py[j];
            
            // Condiciones de contorno periodicas:
            // calcula las distancias entre i y la j-imagen mas cercana.
            dx -= l*round(dx/l);
            dy -= l*round(dy/l);
            
            double r2 = dx*dx + dy*dy; // |ri-rj|^2
            
            // Calula las aceleraciones entre las i's y las j's que esten 
            // dentro del radio corte al cuadrado (rc2).
            if (r2 < rc2)
            {
                double factor = 48*(1/pow(r2,7) - 0.5*(1/pow(r2,4)));
            	xa[i] += factor*dx; 
            	ya[i] += factor*dy;
            	xa[j] -= factor*dx; 
            	ya[j] -= factor*dy;
            	
            	// Actualiza la energia potencial U.
            	U += 4*(1/pow(r2, 6) - 1/pow(r2, 3)) - Uc;
            }     	
    	}
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
    
    // Verifica que las particulas se manten-
    // gan dentro de la caja.
    for (int i = 0; i < n_cuerpos; i++)
    {
        if (xpn[i] > l)
            xpn[i] -= l;
        else if (xpn[i] <= 0)
            xpn[i] += l;
            
        if (ypn[i] > l)
            ypn[i] -= l;
        else if (ypn[i] <= 0)
            ypn[i] += l;        
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
    K = 0.0;
    
    // Calcula velocidad v(t+dt) y energia cinetica.
    for (int i = 0; i < n_cuerpos; i++)
    {
        vxn[i] = vx[i] + 0.5*(ax[i] + axn[i])*h;
        vyn[i] = vy[i] + 0.5*(ay[i] + ayn[i])*h;
        
        K += 0.5*masa[i]*(vxn[i]*vxn[i] + vyn[i]*vyn[i]);
    }
    
    // Llena arreglo v(t+dt).
    for (int i = 0; i < n_cuerpos; i++)
    {
    	y_nueva[i + 2*n_cuerpos] = vxn[i];
        y_nueva[i + 3*n_cuerpos] = vyn[i];
    }            
}

void escalar_velocidad(double *y_nueva)
{ 
   double Km = K/n_cuerpos;
    
   double T = Km/kB; // Tempertura instantanea.
    
   double fs = sqrt(T0/T); // Para 2 grados de libertad.
   
   // Escala v(t) acorde a la temperatura.
   for (int i = 0; i < 2*n_cuerpos; i++)
       y_nueva[i + 2*n_cuerpos] *= fs;
}

void energia()
{
    E = 0.0;
    E = (K + U)/n_cuerpos;
    //E = K + U;
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

void salidaVelocidad(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    
    for (int i = 0; i < n_cuerpos; i++)
        of << scientific << setprecision(9) << "\t" << vx[i];
        
    for (int i = 0; i < n_cuerpos; i++)
        of << scientific << setprecision(9) << "\t" << vy[i];
    
    of << endl;
}

void salidaEnergia(const double tiempo, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;
    of << scientific << setprecision(9) << "\t" << E << endl;
}

/*
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

void salidaPCM(const double tiempo, ofstream &of)
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
}
*/

int main()
{
    // Parametros.
    const double h_step  = 0.005; // Tamanio de paso (adimensional).
    const long int Niter = 1000000; // Numero de iteraciones. 
    const int out_cada   = 100;   // Out_cada iteraciones.
    
    n_cuerpos = 500; // Numero de cuerpos.

    // Otras variables.
    const int n_ec = n_cuerpos * 4; // Numero ecuaciones.
    double tiempo = 0.0;

    // Archivos de salida.
    ofstream of_posicion("T01-posicion.dat", ios::out);
    ofstream of_velocidad("T01-velocidad.dat", ios::out);
    ofstream of_energia("T01-energia.dat", ios::out);
    //ofstream of_momentumLineal("momLineal-prueba.dat", ios::out);
    //ofstream of_momentumAngular("momAngular-prueba.dat", ios::out);
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

    // Inicializa velocidades y la 
    // energia cinetica inicial.
    init_velocidad();
    
    // Inicializa aceleraciones y la
    // energia potencial inicial.
    aceleracion(xp, yp, ax, ay, tiempo, a);
 
    energia();
    
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
            salidaVelocidad(tiempo, of_velocidad);
            //salidaMomentumLineal(tiempo, of_momentumLineal);
            //salidaMomentumAngular(tiempo, of_momentumAngular);
            //salidaPCM(tiempo, of_centroMasa);
            //salidaMasa(tiempo, of_masa);
        }
        
        // Algoritmo de velocidad de Verlet:
        act_posicion(y_nueva, tiempo + h_step, h_step); // r(t+dt).
        aceleracion(xpn, ypn, axn, ayn, tiempo + h_step, a_nueva); // a(t+dt).
        act_velocidad(y_nueva, tiempo + h_step, h_step); // v(t+dt).
            
        escalar_velocidad(y_nueva);
        energia();
         
        for (int j = 0; j < n_ec; j++)
        {
            y[j] = y_nueva[j];
            a[j] = a_nueva[j];
        }

        tiempo += h_step; // ti+1.    
    }

    // Libera el espacio de memoria reservado.
    delete[] y, y_nueva, a, a_nueva; 

    return 0;
}


