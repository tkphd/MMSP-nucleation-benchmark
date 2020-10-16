// nucleation.cpp
// Algorithms for PFHub Nucleation Benchmark
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#ifndef NUCLEATION_UPDATE
#define NUCLEATION_UPDATE
#include "MMSP.hpp"
#include <cmath>
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
double free_energy(grid<dim,T>& Grid)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    grid<dim,T> nrgGrid(Grid);

	double dV = 1.0;
	for (int d = 0; d < dim; d++)
		dV *= dx(Grid, d);

	double energy = 0.0;

	for (int n = 0; n < nodes(Grid); n++) {
		vector<int> x = position(Grid, n);
		double phi = Grid(n);
		vector<T> grad_phi = gradient(Grid, x);
		double grad_phi_sq = grad_phi * grad_phi;
        double dE = 0.5 * grad_phi_sq + g(phi) - df * p(phi);
		energy += dE;;
		nrgGrid(n) = dV * dE;
	}

	energy *= dV;

	#ifdef MPI_VERSION
	double local(energy);
	MPI::COMM_WORLD.Allreduce(&local, &energy, 1, MPI_DOUBLE, MPI_SUM);
	#endif

    output(nrgGrid, "energy.dat");

	return energy;
}

template <int dim, typename T>
double solid_frac(grid<dim,T>& Grid)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    double f = 0.0;
    double N = nodes(Grid);

	for (int n = 0; n < nodes(Grid); n++)
		f += Grid(n);

	#ifdef MPI_VERSION
	double local(f);
	MPI::COMM_WORLD.Allreduce(&local, &f, 1, MPI_DOUBLE, MPI_SUM);

    local = N;
	MPI::COMM_WORLD.Allreduce(&local, &N, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	return f / N;
}

void generate(int dim, const char* filename)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	FILE* fh;

	if (dim==2) {
        const double dt = stability * meshres*meshres / 2;
		if (rank == 0) {
			std::cout << "dt = " << dt << std::endl;
			std::cout << "Run " << std::ceil(1.00 / dt) << " steps to hit unit time, "
			          << std::ceil(100. / dt) << " steps to 100, "
                      << std::ceil(200. / dt) << " steps to 200." << std::endl;
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

		// Embed a single seed in the middle of the domain
		for (int n = 0; n < nodes(initGrid); n++) {
			const vector<int> x = position(initGrid, n);
			const double r = meshres * std::sqrt(x * x);
			initGrid(n) = pf_tanh(r, r0);
		}

        ghostswap(initGrid);

		output(initGrid, filename);

		double F = free_energy(initGrid);
        double f = solid_frac(initGrid);

		if (rank == 0) {
			fh = fopen("free_energy.csv", "w+");
			fprintf(fh, "time,energy,fraction\n");
			fprintf(fh, "%f,%f,%f\n", 0.0, F, f);
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
	static double elapsed = 0.0;
    const double dt = stability * meshres*meshres / 2;

	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	ghostswap(oldGrid);

	grid<dim,T> newGrid(oldGrid);

	if (rank == 0)
		fh = fopen("free_energy.csv", "a");

    const int io_steps = 0.5 / dt; // only write stats every half-unit of time

	for (int step = 0; step < steps; step++) {
		if (rank==0)
			print_progress(step, steps);

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

        if (step % io_steps ==0 || step == steps-1) {
            double F = free_energy(oldGrid);
            double f = solid_frac(oldGrid);

            if (rank == 0) {
                fprintf(fh, "%f,%f,%f\n", elapsed, F, f);
                fflush(fh);
            }
        }
	}

	if (rank == 0)
		fclose(fh);

	double dV = 1.0;
	for (int d = 0; d < dim; d++)
		dV *= dx(oldGrid, d);

    grid<dim,T> nrgGrid(oldGrid);
    ghostswap(nrgGrid);
    for (int n = 0; n < nodes(nrgGrid); n++) {
        vector<int> x = position(nrgGrid, n);
        double phi = oldGrid(x);
		vector<T> grad_phi = gradient(oldGrid, x);
		double grad_phi_sq = grad_phi * grad_phi;
		nrgGrid(n) = dV * (0.5 * grad_phi_sq + g(phi) - df * p(phi));
    }
    output(nrgGrid, "energy.dat");
}

} // namespace MMSP

#endif

#include "MMSP.main.hpp"
