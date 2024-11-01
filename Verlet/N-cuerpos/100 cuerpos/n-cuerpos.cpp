// verlet2.cpp
// Aplicacion del metodo de Verlet de velocidad
// para el movimiento gravitacional en 2D.

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
double E, P, L, Pcm; // Energia y momenta del sistema.
int n_cuerpos; // Numero de cuerpos.

// Constantes.
const double G  = 6.6743e-11;
const double ua = 1.495978707e11; //1.50e11;

// Inicializar masas.
void init_masa()
{
    
    /*masa[0] = 7.349e22; // Luna.
    masa[1] = 5.972e24; // Tierra.
    masa[2] = 1.989e30; // Sol.
    masa[3] = 3.302e23; // Mercurio.
    masa[4] = 4.87e24;  // Venus.*/
    
    Las masas inician con el mismo valor.
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
    yp[2] = 0.0; 
    xp[3] = 0.388*ua; // Mercurio.
    yp[3] = 0.0;
    xp[4] = 0.722*ua; // Venus.
    yp[4] = 0.0;*/
    
    
    // srand(time(0));

    // Valores aleatorios entre [-1, 1] UA en un cuadrado xy.
    for (int i = 0; i < n_cuerpos; i++)
    {
        // xp[i] = ((double)2*rand()/RAND_MAX)*ua - ua;
        // yp[i] = ((double)2*rand()/RAND_MAX)*ua - ua;
        xp[i] = random_menos1_a_1(generator) * ua;
        yp[i] = random_menos1_a_1(generator) * ua;
    }   
}

  // Inicializar velocidades.
void init_velocidad()
{
    
    /*vx[0] = 0.0; // Velocidad Luna.
    vy[0] = 3.0802e4;
    vx[1] = 0.0; // Velocidad Tierra.
    vy[1] = 2.9780e4;
    vx[2] = 0.0; // Velocidad Sol.
    vy[2] = 0.0;
    vx[3] = 0.0; // Velocidad Mercurio.
    vy[3] = 4.7872e4;
    vx[4] = 0.0; // Velocidad Venus.
    vy[4] = 3.5000e4;*/
       
    for (int i = 0; i < n_cuerpos; i++)
    {
        double r = sqrt(pow(xp[i], 2) + pow(yp[i], 2));
        double v = sqrt(G * M_PI * n_cuerpos * masa[i] * r)/(2*ua);

        vx[i] = v * (-yp[i]/r + random_menos1_a_1(generator));
        vy[i] = v * (xp[i]/r + random_menos1_a_1(generator));
    }    
}

void energia_Pmomento_Lmomento()
{
    E, P, L = 0.0;

    double K    = 0.0; // Energia cinetica de las is.
    double U_ij = 0.0; // Energia potencial ij.
   
    for (int j = 0; j < n_cuerpos; j++)
    {        
        double v2 = pow(vx[j], 2) + pow(vy[j], 2);  // v^2.
        double ri = sqrt(pow(xp[j], 2) + pow(yp[j], 2)); // r_i
        
        K += 0.5*masa[j]*v2;
        P += masa[j]*sqrt(v2);
        L += masa[j]*ri*sqrt(v2);
        
        for (int i = 0; i < j; i++)
        {
            double r = sqrt(pow((xp[i] - xp[j]), 2) + pow((yp[i] - yp[j]), 2)); 
            U_ij -= G*masa[i]*masa[j]/r;
        }        
    }

    E = K + U_ij;       
}


void nCuerposGrav(double *y, double *px, double *py, double *xv, 
    double *yv, double *xa, double *ya, const double t, double *dydt)
{
    px = y;
    py = y + n_cuerpos;
    xv = y + 2*n_cuerpos;
    yv = y + 3*n_cuerpos;
    xa = dydt;
    ya = dydt + n_cuerpos;
    
    // Aceleracion debido a la gravedad.
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
  
    for(int j = 0; j < n_cuerpos; j++)
    {
        dydt[j            ] = xa[j];
        dydt[j + n_cuerpos] = ya[j];
    }
}

void verletV(double *y,
                double *y_nueva,
                double *a,
                double *a_nueva,
                const int n_ec,
                const double t,
                const double h,
                void (*derivada)(double *, double *, double *, double *,
                    double *, double *, double *, const double, double *))
{
    // Calcula aceleracion a(t).
    if (t == 0.0)
    	derivada(y, xp, yp, vx, vy, ax, ay, t, a); 
    
    // Calcula posicion x(t+dt)
    for (int i = 0; i < 2*n_cuerpos; i++)
       y_nueva[i] = y[i] + y[i + 2*n_cuerpos]*h + 0.5*a[i]*h*h;
    
    // Calcula aceleracion a(t+dt).
    derivada(y_nueva, xpn, ypn, vxn, vyn, axn, ayn, t+h, a_nueva);

    // Calcula velocidad v(t+dt).
    for (int i = 0; i < 2*n_cuerpos; i++)
        y_nueva[i + 2*n_cuerpos] = y[i + 2*n_cuerpos] + 0.5*(a[i] + a_nueva[i])*h;   
}

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


int main()
{
    // Parametros.
    const double h_step  = 6*3600;// Tamanio de paso 6 horas.
    const long int Niter = 7.3e6; // Numero de iteraciones. (5000 anios)  
    const int out_cada   = 1000; // Out_cada iteraciones.
    
    n_cuerpos = 100; // Numero de cuerpos.

    // Otras variables.
    const int n_ec = n_cuerpos * 4; // Numero ecuaciones.
    double tiempo = 0.0;

    // Archivos de salida.
    ofstream of_posicion("posicion-100-cuerpos.dat", ios::out);
    // ofstream of_velocidad("velocidad-100c.dat", ios::out);
    ofstream of_energia("energia-100-cuerpos.dat", ios::out);
    ofstream of_momentumLineal("momLineal-100-cuerpos.dat", ios::out);
    ofstream of_momentumAngular("momAngular-100-cuerpos.dat", ios::out);
    // ofstream of_masa("masa-n-cuerpos-2.dat", ios::out);
    // ofstream of_centroMasa("centroMasa.dat", ios::out);

    // Reserva espacio para y (posiciones y velocidades).
    // y a (aceleraciones).
    double *y	     = new double[n_ec];
    double *y_nueva  = new double[n_ec];
    double *a        = new double[2*n_cuerpos];
    double *a_nueva  = new double[2*n_cuerpos];
    
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

    // Inicializar masa.
    init_masa();

    // Inicializar posicion.
    init_posicion();

    // Inicializar velocidad.
    init_velocidad();
    
    // Asigna la dir. de memoria de la funcion.
    void (*derivada)(double *y, double *xp, double *yp, double *vx, 
    double *vy, double *ax, double *ay, const double t, double *dydt);
    derivada = nCuerposGrav;

    // salidaSol(tiempo, y, n_ec);
    
    // Ciclo de iteraciones.
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Niter; i++)
    {
        if (i%1000 == 0)
        {
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            //cout << duration.count() << "s " << endl;
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
    
        energia_Pmomento_Lmomento();
        verificar_colisiones(tiempo);

        if (i % out_cada == 0)
        {
            salidaSolucion(tiempo, of_posicion);
            //salidaVelocidad(tiempo, of_velocidad);
            salidaEnergia(tiempo, of_energia);
            salidaMomentumLineal(tiempo, of_momentumLineal);
            salidaMomentumAngular(tiempo, of_momentumAngular);
            //salidaPCM(tiempo, of_centroMasa);
            //salidaMasa(tiempo, of_masa);
        }
            
        verletV(y, y_nueva, a, a_nueva, n_ec, tiempo, h_step, derivada);

        for (int j = 0; j < n_ec; j++)
        {
            y[j] = y_nueva[j];
            a[j] = a_nueva[j];
        }

        tiempo += h_step; // ti+1    
    }

    // Libera el espacio de memoria reservado.
    delete[] y, y_nueva, a, a_nueva; 

    return 0;
}


