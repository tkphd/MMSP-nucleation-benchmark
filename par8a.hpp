// Parameters for PFHub Nucleation Benchmark
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#include <cmath>
const double meshres = 0.4;
const double L = 100.0;
const double df = std::sqrt(2) / 30;
const double r_star = std::sqrt(2) / (6 * df);
const double r0 = 0.99 * r_star;
const double stability = 0.125;
