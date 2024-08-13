#include <array>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <random>

typedef std::array<double, 2> Real2D;

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
    
    std::vector<Particle> particles;

    bool save;
    std::ofstream csv_file;


    Real2D global_best_position = {0,0};
    double global_best_value = std::numeric_limits<double>::max();
    void update_global_best();

    // Random number generator members
    std::default_random_engine engine;
    std::uniform_real_distribution<double> distribution;

    void update(Particle& particle);

public:
    // Constructor
    PSO(int n_particles,
        double w, double c1, double c2,
        double (*fun)(Real2D),
        Real2D center = {0,0}, double radius = 1, 
        bool save = false);
    // Destructor
    ~PSO() = default;

    void optimize(int n_iterations);
    void setSave(bool save) { this->save = save; }
};
