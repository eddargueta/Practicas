// 100-cuerpos.cpp
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
uniform_real_distribution<long double> random_menos1_a_1(-1.0, 1.0);

// Variables globales.
double *xp, *yp, *xpn, *ypn; // Posicion en x, y.
double *vx, *vy, *vxn, *vyn; // Velocidad en x, y.
double *ax, *ay, *axn, *ayn; // Aceleracion en x, y.
double *masa; // Masa de cada particula.
double E, P, L; 
// double K, U;
int n_cuerpos; // Numero de cuerpos.

// Constantes.
const double G  = 6.6743e-11;
const double ua = 1.495978707e11; //1.50e11;

// Inicializar masas.
void init_masa()
{
    /*masa[0] = 7.349e22; // Luna.
    masa[1] = 5.972e24; // Tierra.
    masa[2] = 1.989e30; // Sol.*/
    
    // Las masas inician con el mismo valor.
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
    yp[2] = 0.0; */
    
    // Valores aleatorios entre [-1, 1] UA en un cuadrado xy.
    for (int i = 0; i < n_cuerpos; i++)
    {
        xp[i] = random_menos1_a_1(generator) * ua;
        yp[i] = random_menos1_a_1(generator) * ua;
    }
}

// Inicializar velocidades.
void init_velocidad()
{
    //K = 0.0;
    
    /*vx[0] = 0.0; // Velocidad Luna.
    vy[0] = 3.0802e4;
    vx[1] = 0.0; // Velocidad Tierra.
    vy[1] = 2.9780e4;
    vx[2] = 0.0; // Velocidad Sol.
    vy[2] = 0.0;*/
       
    for (int i = 0; i < n_cuerpos; i++)
    {
        double r = sqrt(xp[i]*xp[i] + yp[i]*yp[i]);
        double v = sqrt(G * M_PI * n_cuerpos * masa[i] * r)/(2*ua);

        vx[i] = v * (-yp[i]/r + random_menos1_a_1(generator));
        vy[i] = v * (xp[i]/r + random_menos1_a_1(generator));
    }
    
    //for (int i = 0; i < n_cuerpos; i++)
      //  K += 0.5*masa[i]*(pow(vx[i], 2) + pow(vy[i], 2)); // Energia cinetica inicial.
      
}

void energia_Pmomento_Lmomento()
{
    E, P, L = 0.0;
    //P, L = 0.0;
    
    double K    = 0.0;
    double U_ij = 0.0;
    
    for (int i = 0; i < n_cuerpos; i++) 
    {
        double v = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
        double cross_pro = xp[i]*vy[i]-yp[i]*vx[i]; // Producto cruz r_jxp_j.
        K += 0.5*masa[i]*v*v;
        
        P += masa[i]*v;
        L += masa[i]*abs(cross_pro);
        
        for (int j = i+1; j < n_cuerpos; j++)
        {
            double dx  = xp[i] - xp[j];
            double dy  = yp[i] - yp[j];
            double rij = sqrt(dx*dx + dy*dy); 	
            
            U_ij -= G*masa[i]*masa[j]/rij;
        }
    }
    
    E = K + U_ij;
}

void aceleraciones(double *px, double *py, double *xa, double *ya, 
    const double t, double *ac)
{
    //U = 0.0;
   	
    // Aceleracion debido a la gravedad.
    for (int i = 0; i < n_cuerpos; i++)
    {
        xa[i] = 0;
        ya[i] = 0;

        for (int j = i+1; j < n_cuerpos; j++)
        {
            double dx = px[i] - px[j];
            double dy = py[i] - py[j];
            double r2 = dx*dx + dy*dy;
            
            xa[i] -= (G*masa[j]*dx)/pow(r2, 1.5);
            ya[i] -= (G*masa[j]*dy)/pow(r2, 1.5);
            xa[j] -= (G*masa[i]*dx)/pow(r2, 1.5);
            ya[j] -= (G*masa[i]*dy)/pow(r2, 1.5);
            
            // U -= G*masa[i]*masa[j]/sqrt(r2); // Energia potencial.   
        }   
    }
  
    // Llena los arreglos de a(t)/a(t+dt).
    for(int i = 0; i < n_cuerpos; i++)
    {
        ac[i            ] = xa[i];
        ac[i + n_cuerpos] = ya[i];
    }
}

void act_posiciones(double *y_nueva, const double t, const double h)
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

void act_velocidades(double *y_nueva, const double t, const double h)
{
    //K = 0.0;
    
    // Calcula velocidad v(t+dt) y energia cinetica.
    for (int i = 0; i < n_cuerpos; i++)
    {
        vxn[i] = vx[i] + 0.5*(ax[i] + axn[i])*h;
        vyn[i] = vy[i] + 0.5*(ay[i] + ayn[i])*h;
        
        //K += 0.5*masa[i]*(pow(vxn[i], 2) + pow(vyn[i], 2));
    }
    
    // Llena arreglo v(t+dt).
    for (int i = 0; i < n_cuerpos; i++)
    {
    	y_nueva[i + 2*n_cuerpos] = vxn[i];
        y_nueva[i + 3*n_cuerpos] = vyn[i];
    }          
}

/*void energia()
{
    E = 0.0;
    E = K + U;
}*/

void verificar_colisiones(double t)
{
    double dist_lim = 1.0e7;

    for (int i = 0; i < n_cuerpos; i++)
    {
        if (masa[i] != 0.0)
        {
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
    const double h_step  = 6*3600; // Tamanio de paso (6 horas).
    const long int Niter = 7.3e6; // Numero de iteraciones (5000 anios).  
    const int out_cada   = 1000;   // Out_cada iteraciones.
    
    n_cuerpos = 100; // Numero de cuerpos.

    // Otras variables.
    const int n_ec = n_cuerpos * 4; // Numero ecuaciones.
    double tiempo = 0.0;

    // Archivos de salida.
    ofstream of_posicion("posicion.dat", ios::out);
    // ofstream of_velocidad("velocidad.dat", ios::out);
    ofstream of_energia("energia.dat", ios::out);
    ofstream of_momentumLineal("momLineal.dat", ios::out);
    ofstream of_momentumAngular("momAngular.dat", ios::out);
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
    aceleraciones(xp, yp, ax, ay, tiempo, a);
    
    //energia();
    energia_Pmomento_Lmomento();
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
        act_posiciones(y_nueva, tiempo + h_step, h_step); // r(t+dt).
        aceleraciones(xpn, ypn, axn, ayn, tiempo + h_step, a_nueva); // a(t+dt).
        act_velocidades(y_nueva, tiempo + h_step, h_step); // v(t+dt).
            
        //energia();
        energia_Pmomento_Lmomento();
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


