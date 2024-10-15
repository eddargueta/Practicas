// Simulacion del cinturon de asteroides considerando asteroides grandes
// y pequenios. 100 cuerpos en total.

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

// Prototipado funciones.

// Variables globales.
double *xp, *yp, *xpn, *ypn; // Posicion en x, y.
double *vx, *vy, *vxn, *vyn; // Velocidad en x, y.
double *ax, *ay, *axn, *ayn; // Aceleracion en x, y.
double *masa; // Masa de cada particula.
double E, P, L, Pcm; // Energia y momenta del sistema.
int n_cuerpos; // Numero de cuerpos.

// Constantes.
const double G  = 6.6743e-11;
const double ua = 1.495978707e11; 
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
    //masa[6] = 5.6880e26; // Saturno.
    //masa[7] = 8.6860e25; // Urano.
    //masa[8] = 1.0240e26; // Neptuno.
     
    // Asteroides grandes.
    masa[6]  = 9.4300e20; // Ceres. 
    masa[7] = 2.2000e20; // Palas. 
    masa[8] = 2.7100e20; // Vesta. 
    masa[9] = 8.6700e19; // Higia. 
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
    //xp[6]  = 9.5371*ua; // Saturno.
    //yp[6]  = 0.0;
    //xp[7]  = 19.1913*ua; // Urano.
    //yp[7]  = 0.0;
    //xp[8]  = 30.0689*ua; // Neptuno.
   // yp[8]  = 0.0;
    
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
            double xi = random_menos1_a_1(generator) * (3.65*ua);
	    double yi = random_menos1_a_1(generator) * (3.65*ua);

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
    //vx[6]  = 0.0;      // Saturno.
    //vy[6]  = 9.6202e3;
    //vx[7]  = 0.0;      // Urano.
    //vy[7]  = 6.8100e3;
    //vx[8]  = 0.0;      // Neptuno.
    //vy[8]  = 1.3070e4;
    
    
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

void energia_Pmomento_Lmomento()
{
    E, P, L = 0.0;

    double K   = 0.0; // Energia cinetica sistema.
    double U   = 0.0; // Energia potencial sistema.
    double Pcm = 0.0; // Momento lineal total.

    for (int j = 0; j < n_cuerpos; j++)
    {   
        double vx2_mas_vy2      = pow(vx[j], 2) + pow(vy[j], 2); // Magnitud velocidad^2.
        double sqrt_vx2_mas_vy2 = sqrt(vx2_mas_vy2); // Magnitud velocidad.
        double cross_pro = xp[j]*vy[j]-yp[j]*vx[j]; // Producto cruz r_jxp_j.
        double sqrt_cross_pro = sqrt(pow(cross_pro, 2));
        
        K   += 0.5*masa[j]*vx2_mas_vy2;
        Pcm += masa[j] * sqrt_vx2_mas_vy2; 
        //L   += masa[j]*abs(cross_pro);
        L += masa[j]*sqrt_cross_pro;

        for (int i = 0; i < j; i++)
        {
            double x2_menos_y2 = pow((xp[i] - xp[j]), 2) + pow((yp[i] - yp[j]), 2);
            double sqrt_x2_menos_y2 = sqrt(x2_menos_y2);

            U -= G*masa[i]*masa[j]/sqrt_x2_menos_y2;
        }
    }

    E = K + U;
    P = Pcm;     
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
            
            xa[i] -= G*masa[j]*dx/r3;
            ya[i] -= G*masa[j]*dy/r3;
        }   
    }

    /*
    for (int i = 5; i < n_cuerpos; i++)
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
            
            xa[i] -= G*masa[j]*dx/r3;
            ya[i] -= G*masa[j]*dy/r3;
        }   
    }
    */
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

void salidaVelocidad(const double t, ofstream &of)
{
  
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
    const double h_step = 1000; // Tamanio de paso. (1000s)
    const int Niter     = 50*365*24*3600/1000; // Numero de iteraciones. (50 anios)  
    const int out_cada  = 24*3600/1000;    // Out_cada iteraciones. 
    
    n_cuerpos = 100; // Numero de cuerpos.

    // Otras variables.
    const int n_ec = n_cuerpos * 4; // Numero ecuaciones.
    double tiempo = 0.0;

    // Archivos de salida.
    ofstream of_posicion("posicion-asteroides-prueba-5.dat", ios::out);
    // ofstream of_velocidad("velocidad.dat", ios::out);
    ofstream of_energia("energia-asteroides-prueba-5.dat", ios::out);
    ofstream of_momentumLineal("momLineal-prueba-5.dat", ios::out);
    ofstream of_momentumAngular("momAngular-asteroides-prueba-5.dat", ios::out);
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
        // verificar_colisiones(tiempo);

        if (i % out_cada == 0)
        {
            salidaSolucion(tiempo, of_posicion);
            salidaEnergia(tiempo, of_energia);
            salidaMomentumLineal(tiempo, of_momentumLineal);
            salidaMomentumAngular(tiempo, of_momentumAngular);
            // salidaPCM(tiempo, of_centroMasa);
            // salidaMasa(tiempo, of_masa);
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


