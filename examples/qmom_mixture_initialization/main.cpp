#include "qmomkit/distributions.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <iostream>

int main()
{
    const std::size_t node_count = 3;
    const std::size_t maximum_order = 2 * node_count - 1;
    const double gamma_density = 1.5;
    const double lognormal_density = 0.8;

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::GammaDistribution{2.0, 0.4}, gamma_density}, {qmomkit::LogNormalDistribution{-0.2, 0.5}, lognormal_density}}};

    const qmomkit::RawMomentResult moments = qmomkit::compute_raw_moments(mixture, maximum_order);
    if (!moments.ok())
    {
        std::cerr << "Mixture moment generation failed: " << moments.message << '\n';
        return 1;
    }

    std::cout << "Absolute-density mixture -> raw moments -> QMOM\n";
    std::cout << "component densities: " << gamma_density << " and " << lognormal_density << '\n';
    std::cout << "m0 = " << moments.moments[0] << " = sum(component number densities)\n";

    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(moments.moments, node_count);
    if (!qmom.ok())
    {
        std::cerr << "QMOM failed: " << qmom.message << '\n';
        return 1;
    }

    for (std::size_t node = 0; node < qmom.weights.size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << qmom.weights[node] << " abscissa=" << qmom.abscissas[node] << '\n';
    }

    return 0;
}
