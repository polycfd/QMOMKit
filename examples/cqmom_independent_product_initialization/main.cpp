#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <utility>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    auto configured = qmomkit::config::load_cqmom_prepared(configuration_file);
    if (!configured)
    {
        std::cerr << "CQMOM configuration failed: " << configured.error() << '\n';
        return 1;
    }

    // CqmomPreparedCase is the immutable configuration/problem definition.
    const qmomkit::CqmomPreparedCase cqmom_case = std::move(configured).value();

    // CqmomState is the runtime state passed to phase-space evolution workflows.
    const qmomkit::CqmomState state = cqmom_case.initial_state();

    const qmomkit::PhaseSpaceHierarchy& hierarchy = cqmom_case.hierarchy();
    std::cout << "configuration -> CqmomPreparedCase -> initial_state -> CqmomState\n";
    std::cout << "CqmomPreparedCase: immutable configuration/problem definition\n";
    std::cout << "CqmomState: evolving phase-space evolution state created from the case\n";
    std::cout << "node_counts = {" << hierarchy.node_counts()[0] << ", " << hierarchy.node_counts()[1] << "}\n";
    std::cout << "conditioning_order = {" << hierarchy.conditioning_order()[0] << ", " << hierarchy.conditioning_order()[1] << "}\n";

    for (std::size_t node = 0; node < state.quadrature().weights().size(); ++node)
    {
        std::cout << "node " << node << ": weight=" << state.quadrature().weights()[node];
        for (std::size_t dimension = 0; dimension < state.quadrature().number_of_dimensions(); ++dimension)
        {
            std::cout << " x" << dimension << '=' << state.quadrature().abscissas()[dimension][node];
        }
        std::cout << '\n';
    }

    return 0;
}
