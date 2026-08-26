#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>

int main()
{
    const auto configured = qmomkit::config::load_qmom(std::filesystem::path{"configuration.qmom"});
    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    std::cout << "Configuration-defined absolute-density mixture -> raw moments -> QMOM\n";
    std::cout << "m0 = " << qmom_case.moments[0] << " = sum(component number densities)\n";

    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmom_case);
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
