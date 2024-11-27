// Aplicacion del metodo RK4 al pendulo doble.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

void RK4(double *y, 
    double *y_nueva, 
    const int n_ec,
    const double t,
    const double h,
    void (*derivada)(const double *, const double, double *))
{
    double *k0 = new double[n_ec];
    double *k1 = new double[n_ec];
    double *k2 = new double[n_ec];
    double *k3 = new double[n_ec];
    double *z  = new double[n_ec];

    (*derivada)(y, t, k0); // k0 = f(ti,yi).

    for (int i = 0; i < n_ec; i++)
        z[i] = y[i] + 0.5*h*k0[i]; // yi + 0.5hk0.

    (*derivada)(z, t+0.5*h, k1); // k1 = f(ti+h/2,yi+0.5hk0).

    for (int i = 0; i < n_ec; i++)
        z[i] = y[i] + 0.5*h*k1[i]; // yi + 0.5hk1.

    (*derivada)(z, t+0.5*h, k2); // k2 = f(ti+h/2,yi+0.5hk1).

    for (int i = 0; i < n_ec; i++)
        z[i] = y[i] + h*k2[i]; // yi + hk2.

    (*derivada)(z, t+h, k3); // k3 = f(ti+h,yi+hk2).

    for (int i = 0; i < n_ec; i++)
        y_nueva[i] = y[i] + (h/6.)*(k0[i] + 2*k1[i] + 2*k2[i] + k3[i]);

    delete[] k0;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] z;
}

void pendulo_doble(const double *y, const double t, double *dydt)
{
    const double m = 0.5; // Masa de la barra.
    const double l = 0.3; // Longitud de la barra.
    const double g = 9.8; // Aceleracion gravedad.
    const double ml2 = m*l*l;

    dydt[0] = (6./ml2)*(2.*y[2] - 3.*y[3]*cos(y[0] - y[1]))/(16. - 9.*pow(cos(y[0] - y[1]), 2)); // Tetha1 punto.
    dydt[1] = (6./ml2)*(8.*y[3] - 3.*y[2]*cos(y[0] - y[1]))/(16 - 9.*pow(cos(y[0] - y[1]), 2)); // Tetha2 punto.
    dydt[2] = (-0.5*ml2)*(dydt[0]*dydt[1]*sin(y[0] - y[1]) + (3.*g/l)*sin(y[0])); // p1 punto.
    dydt[3] = (-0.5*ml2)*(-dydt[0]*dydt[1]*sin(y[0] - y[1]) + (g/l)*sin(y[1]));   // p2 punto.
}

void energia_pendulo_doble(const double *y, const double t, double E)
{
    const double m = 0.5; // Masa de la barra.
    const double l = 0.3; // Longitud de la barra.
    const double g = 9.8; // Aceleracion gravedad.
    const double ml2 = m*l*l;

    double K, U;
    double dydt[4];

    pendulo_doble(y, t, dydt);

    K = (1/6)*ml2*(dydt[1]*dydt[1] + 4.*dydt[0]*dydt[0] + 3.*dydt[0]*dydt[1]*cos(y[0] - y[1]));
    U = -0.5*m*g*l*(3*cos(y[0]) + cos(y[1]));

    E = K + U;
}

void salidaSolucion(const double tiempo, const double *y, const int N, ofstream &of)
{
  of << fixed << setprecision(3) << tiempo;

  for(int i = 0; i < N; i++)
    of << scientific << setprecision(9) << "\t" << y[i];

  of << endl;  
}

void salidaEnergia(const double tiempo, const double E, ofstream &of)
{
    of << fixed << setprecision(3) << tiempo;

    of << scientific << setprecision(16) << "\t" << E;

    of << endl;  
}

int main()
{   
    // Parametros.
    const double h_step  = 0.01; // Tamanio de paso.
    const long int Niter = 12000; // Numero de iteraciones. 
    const int out_cada   = 1;   // Out_cada iteraciones.
    
    // Otras variables.
    const int n_ec = 4; // Numero ecuaciones.
    double tiempo  = 0.0;
    double E       = 0.0; // Energia del sistema.

    // Archivo de salida.
    ofstream of_solucion("solucion.dat", ios::out);
    ofstream of_energia("energia.dat", ios::out); 
    
    // Reserva espacio para
    double *y	     = new double[n_ec]; // Angulos y momentos.
    double *y_nueva  = new double[n_ec]; 
    
    // Condiciones iniciales.
    y[0] = 0.5*M_PI;
    y[1] = 0.5*M_PI;
    y[2] = 0.0;
    y[3] = 0.0;
    
    for (int i = 0; i < n_ec; i++)
        y_nueva[i] = 0.0;

    // Puntero a la funcion "derivada".
    void (*derivada)(const double *, const double, double *);
    derivada = pendulo_doble;

    //energia_pendulo_doble(y, tiempo, E);

    // Ciclo principal.
    for (int i = 0; i < Niter; i++)
    {
        energia_pendulo_doble(y, tiempo, E);
        
        if (i % out_cada == 0)
        {
            salidaSolucion(tiempo, y, n_ec, of_solucion);
            salidaEnergia(tiempo, E, of_energia);
        }
        
        RK4(y, y_nueva, n_ec, tiempo, h_step, derivada);

        for (int j = 0; j < n_ec; j++)
            y[j] = y_nueva[j];

        tiempo += h_step; // ti+1.
    }
    
    // Libera el espacio de memoria reservado.
    delete[] y, y_nueva; 
    
    return 0;
}


