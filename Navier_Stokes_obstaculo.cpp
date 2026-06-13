#include <iostream>
#include <fstream>


// Primero creamos una funcion para saber el indice de la memoria.
inline int IDX(int x, int y, int Nx) {
    return x + y * Nx;
}

//Derivada de primer orden en x(Dif. central)
double ddx(double* f, int x, int y, int Nx, double dx) {
    return ( f[IDX(x + 1, y, Nx)] - f[IDX(x - 1, y, Nx)] ) / (2.0 * dx);
}

// Derivada de primer orden en y (Diferencia central)
double ddy(double* f, int x, int y, int Nx, double dy) {
    return ( f[IDX(x, y + 1, Nx)] - f[IDX(x, y - 1, Nx)] ) / (2.0 * dy);
}

// Laplaciano con estrella de 5 picos
double lap(double* f, int x, int y, int Nx, double dy) {
    return ( f[IDX(x + 1, y, Nx)] + f[IDX(x - 1, y, Nx)] + f[IDX(x, y + 1, Nx)] + f[IDX(x, y - 1, Nx)] - 4.0 * f[IDX(x, y, Nx)]) / (dy*dy);
}

// Funcion para aplicar las condiciones de frontera
// Funcion para aplicar las condiciones de frontera
void aplicar_fronteras(double* u, double* v, int Nx, int Ny, double velocidad_tapa) {
    
    // 1. Condicion de velocidad cero en las paredes (pared izquierda y derecha)
    for(int y = 0; y < Ny; y++) {
        u[ IDX(0, y, Nx) ] = 0.0;
        v[ IDX(0, y, Nx) ] = 0.0;
        u[ IDX(Nx-1, y, Nx) ] = 0.0;
        v[ IDX(Nx-1, y, Nx) ] = 0.0;
    }

    // 2. Condicion de velocidades en las tapas superior e inferior
    for(int x = 0; x < Nx; x++) {
        u[ IDX(x, 0, Nx) ] = 0.0;
        v[ IDX(x, 0, Nx) ] = 0.0;
        u[ IDX(x, Ny-1, Nx) ] = velocidad_tapa;
        v[ IDX(x, Ny-1, Nx) ] = 0.0; 
    }

    // 3. Obstáculo central
    // Definimos un bloque cuadrado en el centro de la matriz.
    // Como tu malla es de 300x300, pondremos el bloque entre los píxeles 100 y 150.
    for(int y = 100; y <= 150; y++) {
        for(int x = 100; x <= 150; x++) {
            u[ IDX(x, y, Nx) ] = 0.0; // El fluido choca y se detiene
            v[ IDX(x, y, Nx) ] = 0.0;
        }
    }
}




// Función para guardar el campo de velocidades en CSV
void guardar_csv(double* u, double* v, int Nx, int Ny, double dx, double dy) {
    std::ofstream archivo("resultado_fluido.csv");
    
    // Escribimos la primera línea con los nombres de las columnas
    archivo << "x,y,u,v\n";
    
    // Recorremos TODA la matriz incluyendo bordes para guardar cada punto
    for(int y = 0; y < Ny; y++) {
        for(int x = 0; x < Nx; x++) {
            // Calculamos la posición física real en metros
            double x_fisico = x * dx;
            double y_fisico = y * dy;
            
            // Guardamos: x, y, u, v
            archivo << x_fisico << "," << y_fisico << "," 
                    << u[IDX(x, y, Nx)] << "," << v[IDX(x, y, Nx)] << "\n";
        }
    }
    
    archivo.close();
    std::cout << "Datos guardados en 'resultado_fluido.csv'" << std::endl;
}






int main(){
    // 1. Definimos el tamaño de la caja y constantes físicas
    int Nx = 300;
    int Ny = 300;
    double rho = 1.0;
    double dt = 0.0001;
    double dx = 0.0033;
    double dy = 0.0033;
    int pasos_totales = 20000;     // Para simular 2 segundo entero
    double nu = 0.001;  // Viscosidad cinemática aumentada para estabilidad
    double velocidad_tapa = 1.0;  // Velocidad de la pared superior

    // 2. Reservamos memoria RAM para todos los archiveros
    double* u = new double[Nx * Ny]();
    double* v = new double[Nx * Ny]();
    double* u_star = new double[Nx * Ny]();
    double* v_star = new double[Nx * Ny]();
    double* p = new double[Nx * Ny]();
    double* p_new = new double[Nx * Ny]();
    double* b = new double[Nx * Ny](); 



    // Bucle temporal
    for (int paso = 0; paso < pasos_totales; paso++) {
        
        // PASO A: Calcular velocidad intermedia (u_star, v_star)
        // Usamos ddx, ddy y lap
        // Ciclos para el interior del fluido

        // Ciclos para el interior del fluido
        for(int y = 1; y < Ny - 1; y++) {
            for(int x = 1; x < Nx - 1; x++) {

                // Velocidad actual en esta coordenada
                double u_actual = u[IDX(x, y, Nx)];
                double v_actual = v[IDX(x, y, Nx)];

                // Evolucion temporal ignorando la presion
                u_star[IDX(x, y, Nx)] = u_actual + dt * (-u_actual * ddx(u, x, y, Nx, dx) - v_actual * ddy(u, x, y, Nx, dy) + nu * lap(u, x, y, Nx, dx));
                
                v_star[IDX(x, y, Nx)] = v_actual + dt * (-u_actual * ddx(v, x, y, Nx, dx) - v_actual * ddy(v, x, y, Nx, dy) + nu * lap(v, x, y, Nx, dx));
            }
        }


        // PASO B: Resolver la Ecuacion de Poisson para la presion
        // AQUI ENTRA OPENMP

        // B.1 Primero vamos a ver el lado derecho de la expresión de la ecuación de Poisson
        for(int y = 1; y < Ny - 1; y++){
            for(int x = 1; x < Nx - 1; x++){
                b[IDX(x,y,Nx)] = (rho/dt)*(ddx(u_star, x, y, Nx, dx)+ddy(v_star, x, y, Nx, dy));
            }
        }

        // B.2 Usando el despeje de p_nueva pensando que la p_nueva corresponde al bloque central
        // Y todos los puntos de alrededor son en una fracción infinitesimal de tiempo una p_antigua

        // Solucionador iterativo de Poisson (Jacobi)
        for(int iter = 0; iter < 200; iter++){
            
            // 1. Calcular la nueva matriz de presión
            #pragma omp parallel for // Paralelizamos con Open MP
            for(int y = 1; y < Ny - 1; y++){
                for(int x = 1; x < Nx - 1; x++){
                    p_new[IDX(x, y, Nx)] = 0.25 * (p[IDX(x+1, y, Nx)] + p[IDX(x-1, y, Nx)] + p[IDX(x, y+1, Nx)] + p[IDX(x, y-1, Nx)] - b[IDX(x, y, Nx)] * (dx * dx)); 
                }
            }

            // 2. Actualizar la memoria (Copiar p_new a p)
            #pragma omp parallel for //Paralelizamos con OpenMP
            for(int y = 1; y < Ny - 1; y++){
                for(int x = 1; x < Nx - 1; x++){
                    p[IDX(x, y, Nx)] = p_new[IDX(x, y, Nx)];
                }
            }
        }

        // PASO C: Corregir la velocidad final (u, v) restando el gradiente de presion
        // Usamos ddx y ddy sobre la presion

    for(int y = 1; y < Ny - 1; y++){
        for(int x = 1; x < Nx - 1; x++){

            u[IDX(x, y, Nx)] = u_star[IDX(x, y, Nx)] - (dt/rho) * ddx(p, x, y, Nx, dx);
            v[IDX(x, y, Nx)] = v_star[IDX(x, y, Nx)] - (dt/rho) * ddy(p, x, y, Nx, dy);  
        }
    }


        // Aplicamos las condiciones de frontera para forzar los bordes
        aplicar_fronteras(u, v, Nx, Ny, velocidad_tapa);
    }

    guardar_csv(u, v, Nx, Ny, dx, dy); //guardamos


    delete[] u;
    delete[] v;
    delete[] u_star;
    delete[] v_star;
    delete[] p;
    delete[] p_new;
    delete[] b;

    return 0;


}

