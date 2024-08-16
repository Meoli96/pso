#include <array>
#include <memory>
#include <random>
#include <span>

#include "../include/ThreadPool.hpp"

#define NUM_THREADS std::thread::hardware_concurrency()
//  Lets try this Data Oriented Design all the young kids are talking about


typedef std::array<double, 2> Real2D;

struct Input {
    double w, c1, c2;
    double r1, r2;
    Real2D position;
    Real2D velocity;
    Real2D best_position;
    Real2D best_global_position;

    Input()
        : w(0),
          c1(0),
          c2(0),
          r1(0),
          r2(0),
          position({0, 0}),
          velocity({0, 0}),
          best_position({0, 0}),
          best_global_position({0, 0}) {}
};

struct Out1 {
    Real2D position;
    Real2D velocity;

    Out1() : position({0, 0}), velocity({0, 0}) {}
};

std::span<Out1> update(std::span<Input> particles,
                       std::span<Out1> out = std::span<Out1>()) {
    if (out.empty()) {
        out = std::span<Out1>(new Out1[particles.size()], particles.size());
    }
    for (int i = 0; i < particles.size(); ++i) {
        // Update velocity
        out[i].velocity[0] =
            particles[i].w * particles[i].velocity[0] +
            particles[i].c1 * particles[i].r1 *
                (particles[i].best_position[0] - particles[i].position[0]) +
            particles[i].c2 * particles[i].r2 *
                (particles[i].best_global_position[0] -
                 particles[i].position[0]);
        out[i].velocity[1] =
            particles[i].w * particles[i].velocity[1] +
            particles[i].c1 * particles[i].r1 *
                (particles[i].best_position[1] - particles[i].position[1]) +
            particles[i].c2 * particles[i].r2 *
                (particles[i].best_global_position[1] -
                 particles[i].position[1]);

        // Update position
        out[i].position[0] += out[i].velocity[0];
        out[i].position[1] += out[i].velocity[1];

        // Update particles
        particles[i].position = out[i].position;
        particles[i].velocity = out[i].velocity;
    }
    return out;
}

double fitness(Real2D x) {
    return 20 + x[0] * x[0] - 10 * cos(2 * M_PI * x[0]) + x[1] * x[1] -
           10 * cos(2 * M_PI * x[1]);
}
struct BestValue {
    Real2D position;
    double value;
};

std::span<BestValue> update_bestValue(std::span<Out1> out,
                                      std::span<BestValue> best) {
    for (int i = 0; i < out.size(); ++i) {
        double f = fitness(out[i].position);
        if (f < best[i].value) {
            best[i].value = f;
            best[i].position = out[i].position;
        }
    }
    return best;
}

double update_global_best(double global_best, Real2D& global_best_position,
                          std::span<BestValue> bests) {
    // Iterate over bests and update global_best
    for (auto best : bests) {
        if (best.value < global_best) {
            global_best = best.value;
            global_best_position = best.position;
        }
    }
    return global_best;
}

void update_random_and_out(std::span<Input> particles, std::span<Out1> out,
                           Real2D global_best_position = {0, 0}) {
    // Random number generator members
    std::default_random_engine engine(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0, 1);

    for (int i = 0; i < particles.size(); ++i) {
        // Generate random numbers
        particles[i].r1 = distribution(engine);
        particles[i].r2 = distribution(engine);
        particles[i].position = out[i].position;
        particles[i].velocity = out[i].velocity;
        particles[i].best_global_position = global_best_position;
    }
}

void generate_particles(std::span<Input> particles, std::span<BestValue> bests,
                        double w, double c1, double c2, Real2D center,
                        double radius) {
    BestValue global_best = {center, std::numeric_limits<double>::max()};

    // Random number generator members
    std::default_random_engine engine(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0, 1);

    for (int i = 0; i < particles.size(); ++i) {
        double r_r = distribution(engine);
        double r_theta = distribution(engine);

        Real2D rand_pos = {center[0] + radius * r_r * cos(2 * M_PI * r_theta),
                           center[1] + radius * r_r * sin(2 * M_PI * r_theta)};

        particles[i].w = w;
        particles[i].c1 = c1;
        particles[i].c2 = c2;
        particles[i].position = rand_pos;
        particles[i].velocity = {0, 0};
        particles[i].best_position = rand_pos;
        particles[i].best_global_position = center;
        bests[i] = {rand_pos, fitness(rand_pos)};
    }
}

int run_dod() {
    int n_particles = 100000;  // 1 hundred thousand
    double global_best = std::numeric_limits<double>::max();
    Real2D global_best_position = {0, 0};
    std::unique_ptr<Input[]> particles = std::make_unique<Input[]>(n_particles);
    std::unique_ptr<Out1[]> out = std::make_unique<Out1[]>(n_particles);
    std::unique_ptr<BestValue[]> bests =
        std::make_unique<BestValue[]>(n_particles);

    // Make span out of unique_ptrs
    std::span<Input> particles_span(particles.get(), n_particles);
    std::span<Out1> out_span(out.get(), n_particles);
    std::span<BestValue> bests_span(bests.get(), n_particles);

    Real2D center = {50, 50};
    double radius = 100;

    // Initialize particles
    generate_particles(particles_span, bests_span, 0.5, 1.0, 1.0, center,
                       radius);
    for (int i = 0; i < 100; ++i) {
        update(particles_span, out_span);
        update_bestValue(out_span, bests_span);
        global_best =
            update_global_best(global_best, global_best_position, bests_span);
        update_random_and_out(particles_span, out_span);
    }
    return 0;
}

void thread_run(std::span<Input> particle, std::span<Out1> out,
                std::span<BestValue> bests, double global_best) {
    update(particle, out);
    update_bestValue(out, bests);
}

int run_parallel() {
    const int n_particles = 100000;  // 1 million
    double global_best = std::numeric_limits<double>::max();
    Real2D global_best_position = {0, 0};
    std::unique_ptr<Input[]> particles = std::make_unique<Input[]>(n_particles);
    std::unique_ptr<Out1[]> out = std::make_unique<Out1[]>(n_particles);
    std::unique_ptr<BestValue[]> bests =
        std::make_unique<BestValue[]>(n_particles);

    // Make span out of unique_ptrs
    std::span<Input> particles_span(particles.get(), n_particles);
    std::span<Out1> out_span(out.get(), n_particles);
    std::span<BestValue> bests_span(bests.get(), n_particles);

    Real2D center = {50, 50};
    double radius = 100;
    // Generate particles
    generate_particles(particles_span, bests_span, 0.5, 1.0, 1.0, center,
                       radius);


    // Initialize thread pool
    thread::ThreadPool<void(std::span<Input>, std::span<Out1>,
                            std::span<BestValue>, double)>
        pool_update(thread_run);
    thread::ThreadPool<void(std::span<Input>, std::span<Out1>, Real2D)> pool(
        update_random_and_out);
    // Split spans into chunks
    for (int it = 0; it < 100; ++it) {
        for (int t = 0; t < NUM_THREADS; t++) {
            int start = t * n_particles / NUM_THREADS;
            int end = (t + 1) * n_particles / NUM_THREADS;
            if (t == NUM_THREADS - 1) {
                end = n_particles;
            }
            std::span<Input> particles_chunk =
                particles_span.subspan(start, end - start);
            std::span<Out1> out_chunk = out_span.subspan(start, end - start);
            std::span<BestValue> bests_chunk =
                bests_span.subspan(start, end - start);
            pool_update.addjob(particles_chunk, out_chunk, bests_chunk,
                               global_best);
        }
        pool_update.wait();
        global_best =
            update_global_best(global_best, global_best_position, bests_span);
        for (int t = 0; t < NUM_THREADS; t++) {
            int start = t * n_particles / NUM_THREADS;
            int end = (t + 1) * n_particles / NUM_THREADS;
            if (t == NUM_THREADS - 1) {
                end = n_particles;
            }
            std::span<Input> particles_chunk =
                particles_span.subspan(start, end - start);
            std::span<Out1> out_chunk = out_span.subspan(start, end - start);
            pool.addjob(particles_chunk, out_chunk, global_best_position);
        }
        pool.wait();
    }
    return 0;
}