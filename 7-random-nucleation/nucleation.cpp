// nucleation.cpp
// Algorithms for PFHub Nucleation Benchmark
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#ifndef NUCLEATION_UPDATE
#define NUCLEATION_UPDATE
#include <cmath>
#include <map>
#include <random>
#include <set>

#include "CImg.h"
#include "MMSP.hpp"
#include "nucleation.hpp"

double g(double x)
{
	return x*x * (1.0 - x)*(1.0 - x);
}

double p(double x)
{
	return x*x*x * (6.0 * x*x - 15.0 * x + 10.0);
}

double g_prime(double x)
{
	return 2.0 * x * (2.0 * x*x - 3.0 * x + 1.0);
}

double p_prime(double x)
{
	return 30.0 * g(x);
}

double pf_tanh(double x, double r)
{
	return 0.5 * (1.0 - std::tanh((x - r) / std::sqrt(2)));
}

namespace MMSP
{

template <int dim, typename T>
void embed_at(grid<dim,T>& Grid, vector<int> X0, double r0)
{
    const double R = 50 * r0;
    vector<int> x(2, 0);
    for (x[1] = X0[1] - R;  x[1] < X0[1] + R; x[1]++) {
        bool outside_y = x[1] < y0(Grid) || x[1] >= y1(Grid);
        for (x[0] = X0[0] - R; x[0] < X0[0] + R; x[0]++) {
            bool outside_x = x[0] < x0(Grid) || x[0] >= x1(Grid);
            if (outside_x ||outside_y)
                continue;
            const vector<int> z = x - X0;
			const double r = meshres * std::sqrt(z * z);
            const double phi = Grid(x) + pf_tanh(r, r0);
            Grid(x) = std::min(phi, 1.0);
        }
    }
}

template <int dim, typename T>
double free_energy(grid<dim,T>& Grid)
{
    grid<dim,T> nrgGrid(Grid);

	double dV = 1.0;
	for (int d = 0; d < dim; d++)
		dV *= dx(Grid, d);

	double energy = 0.0;

	#pragma omp parallel for reduction(+:energy)
	for (int n = 0; n < nodes(Grid); n++) {
		vector<int> x = position(Grid, n);
		double phi = Grid(n);
		vector<T> grad_phi = gradient(Grid, x);
		double grad_phi_sq = grad_phi * grad_phi;
        double dE = 0.5 * grad_phi_sq + g(phi) - df * p(phi);
		energy += dE;
		nrgGrid(n) = dV * dE;
	}

	energy *= dV;

	#ifdef MPI_VERSION
	double local(energy);
	MPI_Allreduce(&local, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif

    output(nrgGrid, "energy.dat");

	return energy;
}

template <int dim, typename T>
double solid_frac(grid<dim,T>& Grid)
{
    double f = 0.0;
    double N = nodes(Grid);

	for (int n = 0; n < nodes(Grid); n++)
		f += Grid(n);

	#ifdef MPI_VERSION
	double local(f);
	MPI_Allreduce(&local, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    local = N;
	MPI_Allreduce(&local, &N, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif

	return f / N;
}

template <int dim, typename T>
int count_particles(grid<dim,T>& Grid)
{
    int x, y;
    std::set<short> shades;

    cimg_library::CImg<short> img(xlength(Grid), ylength(Grid));

    #pragma omp parallel for
    for (int n = 0; n < nodes(Grid); n++) {
        vector<int> X = position(Grid, n);
        X[0] -= g0(Grid, 0);
        X[1] -= g0(Grid, 1);
        const short phi = (Grid(n) < 0.5) ? 0 : 1;
        img.atXY(X[0], X[1]) = phi;
    }

    cimg_library::CImg<short> lbl(img);
    lbl.label(0, 0);

    cimg_forXY(img, x, y) {
        if (img(x, y) > 0)
            shades.insert(lbl(x, y));
    }

    return shades.size();
}

void generate(int dim, const char* filename)
{
	int rank = 0;
	#ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	FILE* fh;

	if (dim==2) {
        const double dt = stability * meshres*meshres / 2;
		if (rank == 0) {
			std::cout << "r* = " << r_star << " and r0 = " << r0 << std::endl;
			std::cout << "dt = " << dt << std::endl;
			std::cout << "Run " << std::ceil(1.00 / dt) << " steps to hit unit time, "
			          << std::ceil(200. / dt) << " steps to 200, "
                      << std::ceil(600. / dt) << " steps to 600." << std::endl;
		}

		const int half_domain = std::ceil(L / (2.0 * meshres));

		GRID2D initGrid(0, -half_domain, half_domain, -half_domain, half_domain);

		for (int d = 0; d < dim; d++) {
            dx(initGrid, d) = meshres;
			if (x0(initGrid, d) == g0(initGrid, d))
				b0(initGrid, d) = Neumann;
			if (x1(initGrid, d) == g1(initGrid, d))
				b1(initGrid, d) = Neumann;
		}

        // flatten the field
        #pragma omp parallel for
		for (int n = 0; n < nodes(initGrid); n++)
            initGrid(n) = 0.0;

        // Generate random positions and times for 100 seeds

        std::random_device                     rand_dev;
        std::mt19937                           generator(rand_dev());
        std::uniform_real_distribution<double> t_distr(0.0, 600.0);
        std::uniform_int_distribution<int>     x_distr(-half_domain, half_domain);

		if (rank == 0) {
			fh = fopen("seeds.csv", "w+");

            for (int i = 0; i < 100; i++) {
                const double t = t_distr(generator);
                const int x = x_distr(generator);
                const int y = x_distr(generator);

                fprintf(fh, "%lf,%d,%d\n", t, x, y);
            }

            fclose(fh);
        }

        // Embed one random seed

        vector<double> x(dim, 0.0);

        x[0] = x_distr(generator);
        x[1] = x_distr(generator);

		#ifdef MPI_VERSION
        MPI_Bcast(&x[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&x[1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif

        embed_at(initGrid, x, r0);

        ghostswap(initGrid);

		output(initGrid, filename);

        int count = count_particles(initGrid);

		double F = free_energy(initGrid);
        double f = solid_frac(initGrid);

		if (rank == 0) {
			fh = fopen("free_energy.csv", "w+");
			fprintf(fh, "time,energy,fraction,particles\n");
			fprintf(fh, "%lf,%lf,%lf,%d\n", 0.0, F, f, count);
			fclose(fh);
		}
	}

	else {
		if (rank == 0)
			std::cout << "Error: The Benchmark is only defined for 2D." << std::endl;
		Abort(EXIT_FAILURE);
	}
}

template <int dim, typename T>
void update(grid<dim,T>& oldGrid, int steps)
{
	FILE* fh;
	int rank = 0;
	#ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif
	static double elapsed = 0.0;
    const double io_dt = 1.0; // only write stats every unit of time
	static double io_elapsed = io_dt - 1.0e-10;
    const double dt = stability * meshres*meshres / 2;

    // Read in set of seeds
    std::map<double, vector<int> > seeds;

    fh = fopen("seeds.csv", "r");

    for (int i = 0; i < 100; i++) {
        vector<int> x(2);
        double t;
        fscanf(fh, "%lf,%d,%d\n", &t, &x[0], &x[1]);
        seeds[t] = x;
    }

    fclose(fh);

    std::map<double, vector<int> >::const_iterator seed = seeds.begin();
    while (seed != seeds.end() && elapsed >= seed->first)
        seed++;

	grid<dim,T> newGrid(oldGrid);

	if (rank == 0)
		fh = fopen("free_energy.csv", "a");

	for (int step = 0; step < steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		#pragma omp parallel for
		for (int n = 0; n < nodes(oldGrid); n++) {
            vector<int> x = position(oldGrid, n);
			T phi = oldGrid(n);
			newGrid(n) = phi + dt * (laplacian(oldGrid, x)
			                         - g_prime(phi)
			                         + df * p_prime(phi));
		}

		swap(oldGrid, newGrid);
		ghostswap(oldGrid);

		elapsed += dt;

        while (seed != seeds.end() && elapsed >= seed->first) {
            embed_at(oldGrid, seed->second, r0);
            ghostswap(oldGrid);
            seed++;
        }

        if (elapsed >= io_elapsed || step == steps-1) {
        	io_elapsed += io_dt;
            double F = free_energy(oldGrid);
            double f = solid_frac(oldGrid);
            int count = count_particles(oldGrid);

            if (rank == 0) {
                fprintf(fh, "%lf,%lf,%lf,%d\n", elapsed, F, f, count);
                fflush(fh);
            }
        }
	}

	if (rank == 0)
		fclose(fh);

	double dV = 1.0;
	for (int d = 0; d < dim; d++)
		dV *= dx(oldGrid, d);

    #ifdef DEBUG
    grid<dim,T> nrgGrid(oldGrid);
    ghostswap(nrgGrid);

	#pragma omp parallel for
    for (int n = 0; n < nodes(nrgGrid); n++) {
        vector<int> x = position(nrgGrid, n);
        double phi = oldGrid(x);
		vector<T> grad_phi = gradient(oldGrid, x);
		double grad_phi_sq = grad_phi * grad_phi;
		nrgGrid(n) = dV * (0.5 * grad_phi_sq + g(phi) - df * p(phi));
    }

    output(nrgGrid, "energy.dat");
    #endif
}

} // namespace MMSP

#endif

#include "MMSP.main.hpp"
