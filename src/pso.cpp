#include "../include/pso.h"




PSO::PSO(int n_particles, double w, double c1, double c2, double (*fun)(Real2D),
         Real2D center, double radius, bool save)
    : n_particles(n_particles),
      w(w),
      c1(c1),
      c2(c2),
      objective_function(fun),
      save(save),
      engine(std::random_device{}()),
      distribution(0, 1) {
    
    // Initialize vector
    this->particles.resize(n_particles);

 
    // Initialize particles
    for (int i = 0; i < n_particles; ++i) {
           // Generate random numbers 
        double rx = this->distribution(engine);
        double ry = this->distribution(engine);
        // !TODO This could use a dedicated linear algebra library
        particles[i].position[0] = center[0] + radius * (2.0 * rx - 1.0);
        particles[i].position[1] = center[1] + radius * (2.0 * ry - 1.0);

        particles[i].best_position = particles[i].position;
        particles[i].best_value = fun(particles[i].position);

        if (particles[i].best_value < global_best_value) {
            global_best_value = particles[i].best_value;
            global_best_position = particles[i].best_position;
        }
    }
    std::cout << "Global best value: " << global_best_value << std::endl;
}
void PSO::update(Particle& particle) {
    // Generate random numbers
    double r1 = this->distribution(engine); 
    double r2 = this->distribution(engine);
    
    // Update velocity
    particle.velocity[0] = w * particle.velocity[0] + 
                            c1 * r1 * (particle.best_position[0] - particle.position[0]) +
                            c2 * r2 * (global_best_position[0] - particle.position[0]);
    particle.velocity[1] = w * particle.velocity[1] +
                            c1 * r1 * (particle.best_position[1] - particle.position[1]) +
                            c2 * r2 * (global_best_position[1] - particle.position[1]);             

    // Update position
    particle.position[0] += particle.velocity[0];
    particle.position[1] += particle.velocity[1];
}
void PSO::update_global_best() {
    for (int i = 1; i < n_particles; ++i) {
        if (particles[i].best_value < global_best_value) {
            global_best_value = particles[i].best_value;
            global_best_position = particles[i].best_position;
        }
    }
}

void PSO::optimize(int n_iterations) {
    // If csv is true, open file, write header (n_particles, n_iterations, w,
    // c1, c2) Each successive line will be
    // Particle[i].position, f(Particle[i].position)
    if (save) {
        csv_file.open("pso.csv");
        csv_file << n_particles << "," << n_iterations << "," << w << "," << c1
                 << "," << c2 << "\n";
    }
    // Iterations
    for (int i = 0; i < n_iterations; ++i) {
        // Particles
        for (int j = 0; j < n_particles; ++j) {
            // Update particle
            update(particles[j]);

            // Check if the new position is better than local and global best
            double f = objective_function(particles[j].position);
            if (f < particles[j].best_value) {
                particles[j].best_value = f;
                particles[j].best_position = particles[j].position;
            }
            if (save && csv_file.is_open()) {
                // Write current particle
                csv_file << particles[j].position[0] << ","
                         << particles[j].position[1] << "," << f << ",";
            }
        }

        // Update global best for next iteration
        update_global_best();

        if (save && csv_file.is_open()) {
            // Go to newline in csv file
            csv_file << "\n";
        }
    }
}



void run_python() {
    // Run python script
    std::system("python3 plot.py");
}