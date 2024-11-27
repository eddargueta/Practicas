// Aplicacion del metodo de diferencias finitas
// a la ecuacion de calor en 2D.
// Para correr con varios nucleos y optimizacion (-O2):
// 1. Compilar con g++ -O2 -fopenmp -o ec-calor-2D.out ec-calor-2D.cpp.
// 2. Escribir en la terminal OMP_NUM_THREADS=X (X: numero de nucleos).
// 3. ./ec-calor-2D.out.

#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

double f_cond_ini(double x, double y)
{
    return 0.0;
}

void condiciones_iniciales(double **u, int Nx, int Ny, double *x, double *y, double dt)
{
    for(int i = 0; i < Nx+1; i++)
    {
        for(int j = 0; j < Ny+1; j++)
            u[i][j] = f_cond_ini(x[i], y[j]);
    }
}

double p1(double x, double t)
{
    return 0.0;
}

double p2(double x, double t)
{
    return sin(M_PI*x);
    // return 0.0;
}

double q1(double y, double t)
{
    return 0.0;
}

double q2(double y, double t)
{
    return 0.0;
}

void condiciones_frontera(double **UU, int Nx, int Ny, double *x, double *y, double tiempo)
{
    // Borde inferior y superior.
    #pragma omp parallel for
    for (int i = 0; i < Nx+1; i++)
    {
        UU[i][0]  = p1(x[i], tiempo);
        UU[i][Ny] = p2(x[i], tiempo); 
    }
    
    // Borde izquierdo y derecho.
    #pragma omp parallel for
    for (int j = 0; j < Ny+1; j++)
    {
        UU[0][j]  = q1(y[j], tiempo);
        UU[Nx][j] = q2(y[j], tiempo); 
    }
}

void salidaSolucion(ostream &of, double **u, double *x, double *y, double t, int Nx, int Ny)
{
    for(int i = 0; i < Nx+1; i++)
    {
        for(int j = 0; j < Ny+1; j++)
            of << t << "\t" << x[i] << "\t" << y[j] << "\t" << u[i][j] << endl;

        of << endl;
    }
    // Separa los bloques de datos.
    of << endl;
}

int main()
{
    int Nx = 100; // Numero de puntos en x.
    int Ny = 100; // Numero de puntos en y.
    int out_cada   = 1; // Output cada no. de iteraciones.
    double Lx      = 1.0; // Longitud en el dominio x.
    double Ly      = 1.0; // Longitud en el dominio x.
    double dx      = Lx/Nx; // dx = dy;
    double alfa    = 0.25;
    double c_calor = 1.0; // Coeficiente de dist. de calor.
    double dt      = (alfa*dx*dx)/(c_calor)*(c_calor);
    int Nt = 2000; // Numero de iteraciones en el tiempo.
    double tiempo  = 0.0; // LLeva la cuenta del tiempo.

    ofstream of;
    of.open("solucion-1.dat", ios::out);

    // Informacion sobre los datos.
    of << "alfa = 0.25" << endl;
    of << "# columna 1: tiempo" << endl;
    of << "# columna 2: coordenada x" << endl;
    of << "# columna 3: coordenada y" << endl;
    of << "# columna 4: u(x,y,t)" << endl;

    // Variables para u.
    double **u_nueva = new double*[Nx+1]; // u_{i,j+1}.
    double **u       = new double*[Nx+1]; // u_{i,j}.
    double *x        = new double[Nx+1]; // Coordenada x.
    double *y        = new double[Ny+1]; // Coordenada y.

    // Asigna espacio a cada puntero.
    for(int i = 0; i < Nx+1; i++)
    {
        u_nueva[i] = new double[Nx+1];
        u[i]       = new double[Nx+1];
    }

    // Coordenada x.
    for (int i = 0; i < Nx+1; i++)
        x[i] = i*dx;
        
    // Coordenada y.
    for(int i = 0; i < Ny+1; i++)
        y[i] = i*dx;

    // Condicion inicial u_{ijk}.
    condiciones_iniciales(u, Nx, Ny, x, y, dt);
    
    // Condiciones de frontera.
    condiciones_frontera(u, Nx, Ny, x, y, tiempo);
    
    // Salida en t=0.
    salidaSolucion(of, u, x, y, tiempo, Nx, Ny);
    			
    // Ciclo principal.
    for (int k = 0; k <= Nt; k++) // Iteraciones en el tiempo.
    {
        // Calcula los puntos interiores en u_{ij,k+1}.
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < Nx; i++) // Iteraciones en x.
        {
            for (int j = 1; j < Ny; j++) // Iteraciones en y.
            	u_nueva[i][j] = u[i][j] + alfa*
            		(u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] - 4.*u[i][j]);
        }

        // Calula los puntos en los bordes en u_{ij,k+1}.
        condiciones_frontera(u_nueva, Nx, Ny, x, y, tiempo + dt);
        
        // Intercambio entre las Us para la siguiente iteracion.
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx+1; i++)
        {
            for (int j = 0; j < Ny+1; j++)
                u[i][j] = u_nueva[i][j];
        }

        tiempo += dt;

        // Salida.
        if (k % out_cada == 0)
        {
            salidaSolucion(of, u, x, y, tiempo, Nx, Ny);
            cout << "it = " << k << " / " << Nt << endl;
        }
    }
   
    return 0;
}


