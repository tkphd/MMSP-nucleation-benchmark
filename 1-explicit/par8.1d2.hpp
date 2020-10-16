// Parameters for PFHub Nucleation Benchmark
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#include <cmath>

const double meshres = 0.2;
const double L = 100.0;
const double df = std::sqrt(2) / 30;
const double r_star = std::sqrt(2) / (6 * df);
const double r0 = 1.01 * r_star;
const double stability = 0.125;
