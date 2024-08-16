#include <array>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <random>
#include <memory>
#include <span>

typedef std::array<double, 2> Real2D;

struct ParticleSave {
    Real2D position;
    double f;
};

std::string outFolder = "../out/";

class Particle {
   public:
    Real2D position;
    Real2D velocity;
    
    Real2D best_position;
    double best_value;

    Particle() : position({0,0}), velocity({0,0}), 
    best_position({0,0}), best_value(std::numeric_limits<double>::max()) {}
};


class PSO { 
private:
    int n_particles;
    double w, c1, c2;
    double (*objective_function)(Real2D);

    Real2D global_best_position = {0,0};
    double global_best_value = std::numeric_limits<double>::max();
    
    // Particles
    std::vector<Particle> particles;
    std::unique_ptr<ParticleSave[]> particle_buf = nullptr;

    bool save;

    // Random number generator members
    std::default_random_engine engine;
    std::uniform_real_distribution<double> distribution;
    
    void update_global_best();
    double update(Particle& particle);
    void initialize(Real2D center, double radius);
    double update_parallel(std::span<Particle> particles, double& thread_best_f);

public:
    // Constructor
    PSO(int n_particles,
        double w, double c1, double c2,
        double (*fun)(Real2D),
        Real2D center = {0,0}, double radius = 1, 
        bool save = false);
    // Destructor
    ~PSO() = default;

    // Setters
    void setSave(bool save) { this->save = save; }
    void setNParticles(int n_particles) { this->n_particles = n_particles; }
    void setParameters(double w, double c1, double c2) { this->w = w; this->c1 = c1; this->c2 = c2; }
    void setObjective(double (*fun)(Real2D)) { this->objective_function = fun; }
    
    // Public methods
    void optimize(int n_iterations);
    void optimize_parallel(int n_iterations);
    void saveToFile(std::string filename, int n_iterations);
    // Reset method should take same arguments as constructor 
    // defaulted to the current values
    


};
