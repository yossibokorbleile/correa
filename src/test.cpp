#include "hera/wasserstein.h"
#include "hera/wasserstein/auction_params.h"


using Diagram = std::vector<std::pair<double, double>>;

std::vector<std::pair<double, double>> pd1;
std::vector<std::pair<double, double>> pd2;
hera::AuctionParams<double> params;
params.wasserstein_power = 2;  // which W_q you want

double w_q_dist = hera::wasserstein_dist<std::vector<std::pair<double, double>>>(pd1, pd2, params);
