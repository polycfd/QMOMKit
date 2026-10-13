#include "qmomkit/distributions.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <iostream>

int main()
{
    // An n-node QMOM inversion needs the moments m0 .. m(2n-1), so the highest
    // order to generate follows from the node count.
    const std::size_t node_count = 2;
    const std::size_t maximum_order = 2 * node_count - 1;
    const double number_density = 2.5;

    // Step 1: turn an analytical distribution into raw moments. QMOMKit scales
    // them by the number density rather than normalizing it away, so m0 comes
    // back as the total number density itself.
    const qmomkit::DistributionSpecification distribution = qmomkit::LogNormalDistribution{-0.2, 0.4};
    const qmomkit::RawMomentResult moments = qmomkit::compute_raw_moments(distribution, maximum_order, number_density);
    if (!moments.ok())
    {
        std::cerr << "Moment generation failed: " << moments.message << '\n';
        return 1;
    }

    std::cout << "Log-normal distribution -> raw moments -> QMOM\n";
    std::cout << "number_density = " << number_density << '\n';
    std::cout << "m0 = " << moments.moments[0] << " (QMOMKit does not normalize it away)\n";
    for (std::size_t order = 0; order < moments.moments.size(); ++order)
    {
        std::cout << "m" << order << " = " << moments.moments[order] << '\n';
    }

    // Step 2: invert those moments into a quadrature. Both steps report success
    // through ok(); a valid distribution can still produce moments that the
    // requested node count cannot represent, so check each one.
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
