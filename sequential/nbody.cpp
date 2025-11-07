#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include "omp_loop.hpp"
#include <omp.h> // for timing

double G = 6.674*std::pow(10,-11);

struct simulation {
    size_t nbpart;

    std::vector<double> mass;

    // position
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    // velocity
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;

    // force
    std::vector<double> fx;
    std::vector<double> fy;
    std::vector<double> fz;

    simulation(size_t nb)
        : nbpart(nb), mass(nb),
          x(nb), y(nb), z(nb),
          vx(nb), vy(nb), vz(nb),
          fx(nb), fy(nb), fz(nb)
    {}
};

void random_init(simulation& s) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution dismass(0.9, 1.);
    std::normal_distribution dispos(0., 1.);

    for (size_t i = 0; i < s.nbpart; ++i) {
        s.mass[i] = dismass(gen);

        s.x[i] = dispos(gen);
        s.y[i] = dispos(gen);
        s.z[i] = 0.;

        s.vx[i] = s.y[i]*1.5;
        s.vy[i] = -s.x[i]*1.5;
        s.vz[i] = 0.;
    }
}

void init_solar(simulation& s) {
    enum Planets {SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON};
    s = simulation(10);

    s.mass[SUN] = 1.9891e30;
    s.mass[MERCURY] = 3.285e23;
    s.mass[VENUS] = 4.867e24;
    s.mass[EARTH] = 5.972e24;
    s.mass[MARS] = 6.39e23;
    s.mass[JUPITER] = 1.898e27;
    s.mass[SATURN] = 5.683e26;
    s.mass[URANUS] = 8.681e25;
    s.mass[NEPTUNE] = 1.024e26;
    s.mass[MOON] = 7.342e22;

    double AU = 1.496e11;

    s.x = {0, 0.39*AU, 0.72*AU, 1.0*AU, 1.52*AU, 5.20*AU, 9.58*AU, 19.22*AU, 30.05*AU, 1.0*AU + 3.844e8};
    s.y = {0,0,0,0,0,0,0,0,0,0};
    s.z = {0,0,0,0,0,0,0,0,0,0};

    s.vx = {0,0,0,0,0,0,0,0,0,0};
    s.vy = {0,47870,35020,29780,24130,13070,9680,6800,5430,29780+1022};
    s.vz = {0,0,0,0,0,0,0,0,0,0};
}

void reset_force(simulation& s) {
    for (size_t i = 0; i < s.nbpart; ++i) {
        s.fx[i] = 0.;
        s.fy[i] = 0.;
        s.fz[i] = 0.;
    }
}

// Compute the force on particle i due to particle j (thread-local)
inline void compute_force_on_i(const simulation& s, size_t i, size_t j, double &fx, double &fy, double &fz) {
    double softening = 0.1;
    double dx = s.x[j] - s.x[i];
    double dy = s.y[j] - s.y[i];
    double dz = s.z[j] - s.z[i];
    double dist_sq = dx*dx + dy*dy + dz*dz + softening*softening;
    double dist = std::sqrt(dist_sq);
    double F = G * s.mass[i] * s.mass[j] / dist_sq;

    fx += F * dx / dist;
    fy += F * dy / dist;
    fz += F * dz / dist;
}

void apply_force(simulation& s, size_t i, double dt) {
    s.vx[i] += s.fx[i]/s.mass[i]*dt;
    s.vy[i] += s.fy[i]/s.mass[i]*dt;
    s.vz[i] += s.fz[i]/s.mass[i]*dt;
}

void update_position(simulation& s, size_t i, double dt) {
    s.x[i] += s.vx[i]*dt;
    s.y[i] += s.vy[i]*dt;
    s.z[i] += s.vz[i]*dt;
}

void dump_state(simulation& s) {
    std::cout << s.nbpart << '\t';
    for (size_t i = 0; i < s.nbpart; ++i) {
        std::cout << s.mass[i] << '\t';
        std::cout << s.x[i] << '\t' << s.y[i] << '\t' << s.z[i] << '\t';
        std::cout << s.vx[i] << '\t' << s.vy[i] << '\t' << s.vz[i] << '\t';
        std::cout << s.fx[i] << '\t' << s.fy[i] << '\t' << s.fz[i] << '\t';
    }
    std::cout << '\n';
}

void load_from_file(simulation& s, std::string filename) {
    std::ifstream in(filename);
    size_t nbpart;
    in >> nbpart;
    s = simulation(nbpart);
    for (size_t i = 0; i < s.nbpart; ++i) {
        in >> s.mass[i];
        in >> s.x[i] >> s.y[i] >> s.z[i];
        in >> s.vx[i] >> s.vy[i] >> s.vz[i];
        in >> s.fx[i] >> s.fy[i] >> s.fz[i];
    }
    if (!in.good())
        throw "kaboom";
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "usage: " << argv[0] << " <input> <dt> <nbstep> <printevery> <nbthreads>\n"
                  << "input can be:\n"
                  << "a number (random initialization)\n"
                  << "planet (initialize with solar system)\n"
                  << "a filename (load from file in singleline tsv)\n";
        return -1;
    }

    double dt = std::atof(argv[2]);
    size_t nbstep = std::atol(argv[3]);
    size_t printevery = std::atol(argv[4]);
    int nbthreads = std::atoi(argv[5]);

    simulation s(1);

    // parse input
    size_t nbpart = std::atol(argv[1]);
    if (nbpart > 0) {
        s = simulation(nbpart);
        random_init(s);
    } else {
        std::string inputparam = argv[1];
        if (inputparam == "planet")
            init_solar(s);
        else
            load_from_file(s, inputparam);
    }

    OmpLoop loop;
    loop.setNbThread(nbthreads);
    loop.setGranularity(1);

    double start = omp_get_wtime();

    for (size_t step = 0; step < nbstep; ++step) {
        if (step % printevery == 0)
            dump_state(s);

        reset_force(s);

        // Parallelized force computation (no race conditions)
        loop.parfor(0, s.nbpart, [&](size_t i) {
            double fx_i = 0., fy_i = 0., fz_i = 0.;
            for (size_t j = 0; j < s.nbpart; ++j)
                if (i != j)
                    compute_force_on_i(s, i, j, fx_i, fy_i, fz_i);
            s.fx[i] = fx_i;
            s.fy[i] = fy_i;
            s.fz[i] = fz_i;
        });

        // Parallel integrate velocities & positions
        loop.parfor(0, s.nbpart, [&](size_t i) {
            apply_force(s, i, dt);
            update_position(s, i, dt);
        });
    }

    double end = omp_get_wtime();
    std::cout << "Elapsed time: " << (end - start) << " seconds\n";

    return 0;
}