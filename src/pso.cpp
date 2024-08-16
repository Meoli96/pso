#include "../include/pso.h"

#include <algorithm>
#include <thread>

#include "../include/ThreadPool.hpp"

#define N_THREADS std::thread::hardware_concurrency()

// double thread_fun(Particle* particles, int particles_per_thread,
//                   int particles_last_thread, int t, int i,
//                   double global_best_value, ParticleSave* particle_buf,
//                   double (*objective_function)(Real2D), double w, double c1,
//                   double c2, Real2D global_best_position,
//                   double* thread_best_global_f, int n_particles, bool save

// ) {
//     int start = t * particles_per_thread;
//     int end = (t + 1) * particles_per_thread;

//     double best_global_f = global_best_value;
//     if (t == N_THREADS - 1) {
//         end += particles_last_thread;
//     }
//     for (int j = start; j < end; ++j) {
//         // Update particle
//         update(particles[j]);
//         if (save) {
//             // Save particle
//             this->particle_buf[i * n_particles + j].position =
//                 particles[j].position;
//             this->particle_buf[i * n_particles + j].f =
//                 objective_function(particles[j].position);
//         }
//     }
//     // Ideally perform reduction to update global best
//     for (int j = start; j < end; ++j) {
//         if (particles[j].best_value < best_global_f) {
//             best_global_f = particles[j].best_value;
//         }
//     }
//     thread_best_global_f[t] = best_global_f;
// }

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
    initialize(center, radius);
}

void PSO::initialize(Real2D center, double radius) {
    for (int i = 0; i < n_particles; ++i) {
        // Generate random numbers
        double r_theta = this->distribution(engine);
        double r_r = this->distribution(engine);
        // Spread them around a circle of radius r
        // !TODO This could use a dedicated linear algebra library
        particles[i].position[0] =
            center[0] + radius * r_r * cos(2 * M_PI * r_theta);
        particles[i].position[1] =
            center[1] + radius * r_r * sin(2 * M_PI * r_theta);

        particles[i].best_position = particles[i].position;
        particles[i].best_value = objective_function(particles[i].position);

        if (particles[i].best_value < global_best_value) {
            global_best_value = particles[i].best_value;
            global_best_position = particles[i].best_position;
        }
    }
}
double PSO::update_parallel(std::span<Particle> particles,
                            double& thread_best_f) {
    int n_particles = particles.size();
    double best_f = global_best_value;

    for (int i = 0; i < n_particles; ++i) {
        // Update particle
        double f = update(particles[i]);
        if (f < best_f) {
            best_f = f;
        }
        if (save) {
            // Save particle
            this->particle_buf[i].position = particles[i].position;
            this->particle_buf[i].f = f;
        }
    }
    thread_best_f = best_f;
    return best_f;
}

double PSO::update(Particle& particle) {
    // Return value of the objective function
    // Generate random numbers
    double r1 = this->distribution(engine);
    double r2 = this->distribution(engine);

    // Update velocity
    particle.velocity[0] =
        w * particle.velocity[0] +
        c1 * r1 * (particle.best_position[0] - particle.position[0]) +
        c2 * r2 * (global_best_position[0] - particle.position[0]);
    particle.velocity[1] =
        w * particle.velocity[1] +
        c1 * r1 * (particle.best_position[1] - particle.position[1]) +
        c2 * r2 * (global_best_position[1] - particle.position[1]);

    // Update position
    particle.position[0] += particle.velocity[0];
    particle.position[1] += particle.velocity[1];

    double f = objective_function(particle.position);
    if (f < particle.best_value) {
        particle.best_value = f;
        particle.best_position = particle.position;
    }
    return f;
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
        // Initialize ParticleSave[n_iterations][n_particles]
        this->particle_buf =
            std::make_unique<ParticleSave[]>(n_iterations * n_particles);
    }
    // Iterations
    for (int i = 0; i < n_iterations; ++i) {
        // Particles
        for (int j = 0; j < n_particles; ++j) {
            // Update particle
            update(particles[j]);
            if (save) {
                // Save particle
                this->particle_buf[i * n_particles + j].position =
                    particles[j].position;
                this->particle_buf[i * n_particles + j].f =
                    objective_function(particles[j].position);
            }
        }

        // Update global best for next iteration
        update_global_best();
    }
}

void PSO::optimize_parallel(int n_iterations) {
    // !TODO Implement parallel version of the optimization (this one sucks)

    // If csv is true, open file, write header (n_particles, n_iterations, w,
    // c1, c2) Each successive line will be
    // Particle[i].position, f(Particle[i].position)
    if (save) {
        // Initialize ParticleSave[n_iterations][n_particles]
        this->particle_buf =
            std::make_unique<ParticleSave[]>(n_iterations * n_particles);
    }
    // Get spans for threads
    int particles_per_thread = n_particles / N_THREADS;
    int particles_last_thread = n_particles % N_THREADS;

    // Initialize thread best global f
    double thread_best_global_f[N_THREADS];

    // Initialize threadpool
    thread::ThreadPool<double(std::span<Particle>, double&)> pool(
        [this](std::span<Particle> particles, double& thread_best_f) {
            return this->update_parallel(particles, thread_best_f);
        });
    // Iterations
    for (int i = 0; i < n_iterations; ++i) {
        // Particles
        for (int t = 0; t < N_THREADS; t++) {
            int start = t * particles_per_thread;
            int end = (t + 1) * particles_per_thread;
            if (t == N_THREADS - 1) {
                end += particles_last_thread;
            }
            // Make span out of particles
            std::span t_particles(particles.begin() + start,
                                  particles.begin() + end);

            pool.addjob(t_particles, thread_best_global_f[t]);
        }
        // Wait for threads to finish
        for (int t = 0; t < N_THREADS; t++) {
            pool.wait();
        }
        // Update global best for next iteration
        this->global_best_value = *std::min_element(
            thread_best_global_f, thread_best_global_f + N_THREADS);
    }
}

void PSO::saveToFile(std::string filename, int n_iterations) {
    std::ofstream file(outFolder + filename);

    // Write header
    file << this->n_particles << "," << this->w << "," << this->c1 << ","
         << this->c2 << std::endl;

    // Write data
    for (int i = 0; i < n_iterations; ++i) {
        for (int j = 0; j < n_particles; ++j) {
            file << this->particle_buf[i * n_particles + j].position[0] << ","
                 << this->particle_buf[i * n_particles + j].position[1] << ","
                 << this->particle_buf[i * n_particles + j].f << ",";
        }
        file << "\n";
    }
    // Flush and close
    file.flush();
    file.close();
}

void run_python() {
    // Run python script
    std::string plotFile = "../out/pso.csv";
    std::string path = "../scripts/plot.py";
    std::string command = "python3 " + path + " " + plotFile;
    std::system(command.c_str());
}