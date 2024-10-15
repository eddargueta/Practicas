// Aplicacion del metodo de diferencias finitas
// a la ecuacion de onda.

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f_cond_ini(double x)
{
    double L = 1.0;
    double k = 1.0;

    if (0.0 <= x && x <= L/2)
        return 2*k*x/L;
    else if (L/2 <= x && x <= L)
        return 2*k*(L-x)/L;
    else
        return 0.0;
}

double g_cond_ini(double x)
{
    return 0.0;
}

void condiciones_iniciales(double *u_vieja, double *u, int N, double *x, double dt)
{
    for (int i = 0; i < N+1; i++) // u(x,t=0) = u_{i0}.
        u_vieja[i] = f_cond_ini(x[i]);
    
    for (int i = 0; i < N+1; i++) // v(x,t=0) = u_{i1}.
        u[i] = u_vieja[i] + g_cond_ini(x[i])*dt;
}

double p_cond_frontera(double t)
{
    return 0.0;
}

double q_cond_frontera(double t)
{
    return 0.0;
}

void condiciones_frontera(double *u_nueva, int N, double tiempo)
{
    u_nueva[0] = p_cond_frontera(tiempo);
    u_nueva[N] = q_cond_frontera(tiempo);
}

void salidaSolucion(ostream &of, double *u, double *x, double t, int N)
{
    for(int i = 0; i < N+1; i++)
            of << t << "\t" << x[i] << "\t" << u[i] << endl;

    of << endl << endl;
}

int main()
{
    int N = 100; // Numero de puntos en x.
    int out_cada = 1; // Output cada no. de iteraciones.
    double L     = 1.0; // Longitud en el dominio x.
    double dx    = 0.01; // L/N;
    double alfa  = 1.25;
    double v     = 1.0; // Velocidad de la onda.
    double dt    = alfa*dx/v;
    int Nt = 400; // Numero de iteraciones en el tiempo.
    double tiempo = 0.0; // LLeva la cuenta del tiempo.

    ofstream of;
    of.open("solucion-1d-alpha5.dat", ios::out);

    // Informacion sobre los datos.
    of << "alpha = 1.25" << endl;
    of << "# columna 1: tiempo" << endl;
    of << "# columna 2: coordenada x" << endl;
    of << "# columna 3: u(x,t)" << endl;

    // Variables para u.
    double *u_nueva = new double[N+1]; // u_{i,j+1}.
    double *u       = new double[N+1]; // u_{i,j}.
    double *u_vieja = new double[N+1]; // u_{i,j-1}.
    double *x       = new double[N+1]; // Coordenada x.
    
    // Coordenada x.
    for (int i = 0; i < N+1; i++)
        x[i] = i*dx;

    // Condiciones iniciales u_{i,j-1}, u_{i,j}.
    condiciones_iniciales(u_vieja, u, N, x, dt);
    
    // Condiciones de frontera.
    condiciones_frontera(u_nueva, N, tiempo);
    
    // Salida en t=0.
    salidaSolucion(of, u, x, tiempo, N);

    tiempo += dt;
    
    salidaSolucion(of, u, x, tiempo, N);
    
    // Ciclo principal.
    for (int j = 0; j <= Nt; j++)
    {
        for (int i = 0; i < N; i++)
            u_nueva[i] = 2.*(1.-alfa*alfa)*u[i] + alfa*alfa*(u[i-1]+u[i+1]) - u_vieja[i];

        // Actualiza las condiciones de frontera.
        condiciones_frontera(u_nueva, N, tiempo + dt);
        
        // Intercambio entre las Us para la siguiente iteracion.
        for (int i = 0; i < N+1; i++)
        {
            u_vieja[i] = u[i];
            u[i]       = u_nueva[i];
        }

        tiempo += dt;

        // Salida.
        if (j % out_cada == 0)
        {
            salidaSolucion(of, u, x, tiempo, N);
            cout << "it = " << j << " / " << Nt << endl;
        }
    }
    
    return 0;
}


