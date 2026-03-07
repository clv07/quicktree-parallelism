#include "distancemat.hpp"
#include <tuple>
#include <cstdint>

std::tuple<uint32_t, uint32_t, double, double> compute_min_ij(const DistanceMatrix* mat);
DistanceMatrix update_matrix(const DistanceMatrix* mat);