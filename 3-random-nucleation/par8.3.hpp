// Parameters for PFHub Nucleation Benchmark
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)
#include <cmath>

const double meshres = 0.4;
const double L = 1000.0;
const double df = std::sqrt(2) / 12;
const double r_star = std::sqrt(2) / (6 * df);
const double r0 = 2.2 * r_star;
const double stability = 0.125;

