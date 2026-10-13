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
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = std::move(configured).value();

    // CqmomState is the runtime state passed to phase-space evolution workflows.
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();

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
