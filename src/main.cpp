#include <chrono>
#include <cmath>
#include <iostream>

#ifndef DOD_TEST
#include "../src/pso.cpp"
#endif

#ifdef DOD_TEST
#include "../src/pso_dod.cpp"
#endif

// Objective function
double sphere_function(Real2D x) {
    // X^2 + Y^2
    return x[0] * x[0] + x[1] * x[1];
}

double rastrigin(Real2D x) {
    // A = 10, n = 2
    return 20 + x[0] * x[0] - 10 * cos(2 * M_PI * x[0]) + x[1] * x[1] -
           10 * cos(2 * M_PI * x[1]);
}

double ackley(Real2D x) {
    double A = -20 * exp(-0.2 * sqrt(0.5 * (x[0] * x[0] + x[1] * x[1])));
    double B = -exp(0.5 * (cos(2 * M_PI * x[0]) + cos(2 * M_PI * x[1])));
    double C = 20 + M_E;
    return A + B + C;
}
int main() {
// // Initialize PSO
#ifndef DOD_TEST
    PSO pso(1000000, 0.5, 1.0, 1.0, rastrigin, {250, 250}, 100, false);

    // // Optimize
    // // Start clock to meausre time
    // auto start = std::chrono::high_resolution_clock::now();
    // pso.optimize(100);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = end - start;
    // std::cout << "Elapsed time serial: " << elapsed.count() << "s\n";

    auto start1 = std::chrono::high_resolution_clock::now();
    pso.optimize_parallel(100);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end1 - start1;
    std::cout << "Elapsed time parallel: " << elapsed1.count() << "s\n";

    // // //Save results
    // pso.saveToFile("rastrigin_parallel.csv");

    // run_python("rastrigin_parallel.csv");
#endif

#ifdef DOD_TEST

    // auto start = std::chrono::high_resolution_clock::now();
    // pso_serial(100000, 100, 0.5, 1.0, 1.0, {50, 50}, 100, rastrigin);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = end - start;
    // std::cout << "Elapsed time serial: " << elapsed.count() << "s\n";

    auto start1 = std::chrono::high_resolution_clock::now();
    pso_parallel(1000000, 100, 0.5, 1.0, 1.0, {50, 50}, 100, rastrigin);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end1 - start1;
    std::cout << "Elapsed time parallel: " << elapsed1.count() << "s\n";
#endif
}
