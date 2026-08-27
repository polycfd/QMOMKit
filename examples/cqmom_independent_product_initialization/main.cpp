#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <utility>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    auto configured = qmomkit::config::load_cqmom(configuration_file);
    if (!configured)
    {
        std::cerr << "CQMOM configuration failed: " << configured.error() << '\n';
        return 1;
    }

    // CqmomCase is the immutable configuration/problem definition.
    const qmomkit::CqmomCase cqmom_case = std::move(configured).value();

    // compute_cqmom(...) computes the initial quadrature for that case.
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
    if (!cqmom.ok())
    {
        std::cerr << "CQMOM failed: " << cqmom.message << '\n';
        return 1;
    }

    // CqmomState is the runtime state passed to phase-space evolution workflows.
    const qmomkit::CqmomState state = cqmom_case.initial_state();

    const qmomkit::CqmomQuadrature& quadrature = cqmom.quadrature;
    const qmomkit::PhaseSpaceHierarchy& hierarchy = cqmom_case.hierarchy();
    std::cout << "configuration -> CqmomCase -> compute_cqmom -> initial_state -> CqmomState\n";
    std::cout << "CqmomCase: immutable configuration/problem definition\n";
    std::cout << "CqmomState: evolving phase-space evolution state created from the case\n";
    std::cout << "node_counts = {" << hierarchy.node_counts()[0] << ", " << hierarchy.node_counts()[1] << "}\n";
    std::cout << "conditioning_order = {" << hierarchy.conditioning_order()[0] << ", " << hierarchy.conditioning_order()[1] << "}\n";
    (void)state;

    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            std::cout << " x" << dimension << '=' << quadrature.abscissas()[dimension][node];
        }
        std::cout << '\n';
    }

    return 0;
}
