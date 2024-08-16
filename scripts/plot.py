#!/usr/bin/env python3
class Particle:
    def __init__(self, x, y, f):
        self.x = x
        self.y = y
        self.f = f

import csv

def read_csv(filename):
    # First line is the header, structure is:
    #   0: n_particles, n_iterations, w, c1, c2
    # Each iteration after that is:
    #   : Particle position, f
    # Read firt line
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        n_particles, w, c1, c2 = map(float, header)
        # n_particles and n_iterations are integers
        n_particles = int(n_particles)
        
        
        data = []
        for row in reader:
            # strip last element of row
            row = row[:-1]
            data.append(list(map(float, row)))

    return n_particles, w, c1, c2, data

# Function from data to list of (lists, iterations) particles
def data_to_particles(data, n_particles, n_iterations):
    particles = []
    for it in range(n_iterations):
        cur_it = []
        # Extract it row from data
        row = data[it]
        for i in range(n_particles):
            # Extract particle from row
            particle = Particle(row[i*3], row[i*3+1], row[i*3+2])
            cur_it.append(particle)
        particles.append(cur_it)

    return particles


import matplotlib.pyplot as plt
import matplotlib.animation as animation


def animate(particles, fps):
    fig, ax = plt.subplots()
    line, = ax.plot([], [], 'o', markersize=3)
    ax.set_xlim(-500, 500)
    ax.set_ylim(-500, 500)

    def frame(i):
        x = [p.x for p in particles[i]]
        y = [p.y for p in particles[i]]
        # color plot by f
        c = [p.f for p in particles[i]]
        line.set_data(x, y)
        return line,
    ani = animation.FuncAnimation(fig, frame, frames=len(particles), interval=1000/fps)
    return ani

if __name__ == '__main__':
    import sys
    print("Starting plot.py")
    # Read filename from command line
    if (len(sys.argv) < 2):
        filename = "../out/pso.csv"
    else:
        filename = sys.argv[1]
    
    # Read name of file without extension
    outFile = filename.split("/")[-1].split(".")[0]
    

    outFile = "../out/" + outFile + ".gif"
    n_particles, w, c1, c2, data = read_csv(filename)
    n_iterations = len(data)
    particles = data_to_particles(data, n_particles, n_iterations)
    ani = animate(particles, 5)
    ani.save(outFile, writer='imagemagick', fps=5)
    