// nucleation.cpp
// Algorithms for PFHub Nucleation Benchmark
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#ifndef NUCLEATION_UPDATE
#define NUCLEATION_UPDATE
#include <cmath>
#include <map>
#include <random>
#include <set>

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
T boundary_value(const MMSP::grid<dim, T>& GRID, const vector<int>& x)
{
    const bool on_boundary = (x[1] == g0(GRID, 1) - 1);
    if (!on_boundary)
        return GRID(x);

    const bool on_prewet = (x[0] >= g0(GRID, 0) / 2) && (x[0] < g1(GRID, 0) / 2);
    if (on_prewet)
        return wetted;

    return 0.0;
}

template <int dim, typename T>
T neu_laplacian(const MMSP::grid<dim, T>& GRID, const vector<int>& x)
{
    T laplacian = 0.0;
    MMSP::vector<int> s = x;
    const T& y = GRID(x);

    for (int i=0; i<dim; i++) {
        s[i] += 1;
        const T& yh = GRID(s);
        s[i] -= 2;
        const T yl = boundary_value(GRID, s);
        s[i] += 1;

        double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
        laplacian += weight * (yh - 2.0 * y + yl);
    }
    return laplacian;
}

template <int dim, typename T>
MMSP::vector<T> neu_gradient(const grid<dim, T>& GRID, const vector<int>& x)
{
    vector<T> gradient(dim);
    vector<int> s = x;

    for (int i=0; i<dim; i++) {
        s[i] += 1;
        const T& yh = GRID(s);
        s[i] -= 2;
        const T yl = boundary_value(GRID, s);
        s[i] += 1;

        double weight = 1.0 / (2.0 * dx(GRID, i));
        gradient[i] = weight * (yh - yl);
    }
    return gradient;
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
		vector<T> grad_phi = neu_gradient(Grid, x);
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
			std::cout << "dt = " << dt << std::endl;
			std::cout << "Run " << unsigned(std::ceil(6500. / dt))
                      << " steps to hit 6500 units of time. " << std::endl;
		}

		const int half_domain = std::ceil(L / (2.0 * meshres));

		GRID2D initGrid(0, -half_domain, half_domain, 0, half_domain);

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

        ghostswap(initGrid);

		output(initGrid, filename);

		double F = free_energy(initGrid);
        double f = solid_frac(initGrid);

		if (rank == 0) {
			fh = fopen("free_energy.csv", "w+");
			fprintf(fh, "time,energy,fraction\n");
			fprintf(fh, "%lf,%lf,%lf\n", 0.0, F, f);
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
    const double io_dt = 10.0; // only write stats every 10 units of time
	static double io_elapsed = io_dt - 1.0e-10;
    const double dt = stability * meshres*meshres / 2;

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
			newGrid(n) = phi + dt * (neu_laplacian(oldGrid, x)
			                         - g_prime(phi)
			                         + df * p_prime(phi));
		}

		swap(oldGrid, newGrid);
		ghostswap(oldGrid);

		elapsed += dt;

        if (elapsed >= io_elapsed || step == steps-1) {
        	io_elapsed += io_dt;
            double F = free_energy(oldGrid);
            double f = solid_frac(oldGrid);

            if (rank == 0) {
                fprintf(fh, "%lf,%lf,%lf\n", elapsed, F, f);
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

    double nrg_tot = 0.0;

	#pragma omp parallel for reduction(+:nrg_tot)
    for (int n = 0; n < nodes(nrgGrid); n++) {
        vector<int> x = position(nrgGrid, n);
        double phi = oldGrid(x);
		vector<T> grad_phi = neu_gradient(oldGrid, x);
		double grad_phi_sq = grad_phi * grad_phi;
        double dE = dV * (0.5 * grad_phi_sq + g(phi) - df * p(phi));
		nrgGrid(n) = dE;
        nrg_tot += dE;
    }

    for (int n = 0; n < nodes(nrgGrid); n++)
        nrgGrid(n) /= nrg_tot;

    output(nrgGrid, "energy.dat");
}

} // namespace MMSP

#endif

#include "MMSP.main.hpp"
