// Aplicacion del metodo de diferencias finitas
// a la ecuacion de onda en 2D en distintos me
// dios.

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f_cond_ini(double x, double y)
{
    double Lx  = 1.0; // Longitud en x.
    double Ly  = 1.0; // Longitud en y.
    double w   = 800; // Controla el ancho del pulso.
    double amp = 1.0; // Amplitud del pulso.

    // Pulso en forma de una gaussiana.
    // return amp*exp(-w*(pow(x-0.5*Lx,2) + pow(y-0.5*Ly,2)));
    return 0.0;
}

double g_cond_ini(double x, double y)
{
    double L = 1.0; // Longitud de la cuerda.
    return 0.0;
}

void condiciones_iniciales(double **u_vieja, double **u, int Nx, 
    int Ny, double *x, double *y, double dt)
{
    for(int i = 0; i < Nx+1; i++)
    {
        for(int j = 0; j < Ny+1; j++)
        {
            u_vieja[i][j] = f_cond_ini(x[i], y[j]);
            u[i][j]       = u_vieja[i][j] +  g_cond_ini(x[i], y[j])*dt;
        }
    }
}

double q1(double y, double t)
{
    double Ly  = 1.0; // Longitud en y.
    double w   = 800; // Ancho del pulso.
    double amp = 1.0; // Amplitud de la onda.
    double omega = 50.0;
     
    return exp(-w*pow((y - 0.5*Ly), 2))*sin(omega*t);
    // return 0.0;
}

double p1(double x, double t)
{
    return 0.0;
}

double p2(double x, double t)
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
    for (int i = 0; i < Nx+1; i++)
    {
        UU[i][0]  = p1(x[i], tiempo);
        UU[i][Ny] = p2(x[i], tiempo); 
    }
    
    // Borde izquierdo y derecho.
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
    int Nx = 200; // Numero de puntos en x.
    int Ny = 200; // Numero de puntos en y.
    int out_cada  = 10;
    double Lx     = 1.0; // Longitud del dominio en x.
    double Ly     = 1.0; // Longitud del dominio en y.
    double dx     = Lx/Nx; // dx = dy.
    double alfa   = 0.5;
    // double dt  = alfa*dx/v; // Tamanio de paso en el tiempo.
    double dt;
    double tiempo = 0.0; // Lleva la cuenta del tiempo.  
    int Nt = 2000; // Numero de iteraciones en el tiempo.

    ofstream of;
    of.open("solucion.dat", ios::out);

    // Informacion sobre los datos.
    of << "# columna 1: tiempo" << endl;
    of << "# columna 2: coordenada x" << endl;
    of << "# columna 3: coordenada y" << endl;
    of << "# columna 4: u(x,y,t)" << endl;

    // Variables para u.
    double **u_nueva = new double*[Nx+1]; // u_{ij,k+1}.
    double **u       = new double*[Nx+1]; // u_{ij,k}.
    double **u_vieja = new double*[Nx+1]; // u_{ij,k-1}.
    double **v       = new double*[Nx+1];
    double *x        = new double[Nx+1]; // Coordenada x.
    double *y        = new double[Ny+1]; // Coordenada y.

    // Asigna espacio a cada puntero.
    for(int i = 0; i < Nx+1; i++)
    {
        u_nueva[i] = new double[Nx+1];
        u[i]       = new double[Nx+1];
        u_vieja[i] = new double[Nx+1];
        v[i]       = new double[Nx+1];
    }
    
    // Asigna velocidades.
    // Lado izquierdo.
    for(int i = 0; i < (Nx+1)/2; i++)
    {
    	for(int j = 0; j < Ny+1; j++)
            v[i][j] = 1.0;
    }
  
    // Lado derecho.
    for(int i = (Nx+1)/2; i < Nx+1; i++)
    {
    	for(int j = 0; j < Ny+1; j++)
    	    v[i][j] = 0.4;
    }
  
    // Fijar dt.
    dt = alfa * dx/1.0;
    
    // Coordenada x.
    for(int i = 0; i < Nx+1; i++)
        x[i] = i*dx;

    // Coordenada y.
    for(int i = 0; i < Ny+1; i++)
        y[i] = i*dx;

    // Condiciones iniciales u_{ij,k-1}, u_{ijk}.
    condiciones_iniciales(u_vieja, u, Nx, Ny, x, y, dt);

    // Condiciones de frontera.
    condiciones_frontera(u_vieja, Nx, Ny, x, y, tiempo);
    condiciones_frontera(u, Nx, Ny, x, y, tiempo);
    
    // Salida en t=0.
    salidaSolucion(of, u_vieja, x, y, tiempo, Nx, Ny);
    
    tiempo += dt;
    
    salidaSolucion(of, u, x, y, tiempo, Nx, Ny);

    // Ciclo principal.
    for (int k = 0; k <= Nt; k++) // Iteraciones en el tiempo.
    {
        // Calcula los puntos interiores en u_{ij,k+1}.
        for (int i = 1; i < Nx; i++) // Iteraciones en x.
        {
            for (int j = 1; j < Ny; j++) // Iteraciones en y.
            {
            	alfa = v[i][j]*dt/dx;
                u_nueva[i][j] = 2.*u[i][j] - u_vieja[i][j] + alfa*alfa*(
                        u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] - 4.*u[i][j]);
            }
        }

        // Calula los puntos en los bordes en u_{ij,k+1}.
        condiciones_frontera(u_nueva, Nx, Ny, x, y, tiempo + dt);
        
        // Intercambio entre las Us para la siguiente iteracion.
        for (int i = 0; i < Nx+1; i++)
        {
            for (int j = 0; j < Ny+1; j++)
            {
                u_vieja[i][j] = u[i][j];
                u[i][j]       = u_nueva[i][j];
            }
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
