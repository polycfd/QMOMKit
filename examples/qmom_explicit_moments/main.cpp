#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <iostream>

int main()
{
    const std::size_t node_count = 2;

    // Raw moments m0..m3 for weights {0.4, 0.6} at abscissas {-1, 2}.
    const qmomkit::MomentSet moments{1.0, 0.8, 2.8, 4.4};

    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, node_count);
    if (!result.ok())
    {
        std::cerr << "QMOM failed: " << result.message << '\n';
        return 1;
    }

    std::cout << "Explicit raw moments -> QMOM\n";
    std::cout << "m0 is the total weight: " << moments[0] << '\n';
    for (std::size_t node = 0; node < result.weights.size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << result.weights[node] << " abscissa=" << result.abscissas[node] << '\n';
    }

    return 0;
}
