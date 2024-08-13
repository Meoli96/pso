#include "../src/pso.cpp"

// Objective function
double sphere_function(Real2D x) {
    // X^2 + Y^2
    return x[0] * x[0] + x[1] * x[1];
}

int main() {
    // Initialize PSO
    PSO pso(100, 
        0.5, 1.0, 1.0,
        sphere_function,
        {4, 4}, 4, true);

    // Optimize
    pso.optimize(100);

}