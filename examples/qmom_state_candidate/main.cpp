#include "qmomkit/qmom.hpp"
#include "qmomkit/realizability_policy.hpp"

#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

int main()
{
    // Weights 0.5 at abscissas 1 and 3; two nodes need the moments m0..m3.
    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet{1.0, 2.0, 5.0, 14.0};
    qmom_case.node_count = 2;
    qmomkit::QmomState state = qmom_case.initial_state();

    // Accept only moments of a distribution on [0, inf), such as particle sizes.
    qmomkit::RealizabilityPolicy nonnegative;
    nonnegative.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    nonnegative.support = qmomkit::MomentSupport::nonnegative;

    // Moments as the application might have changed them.
    const std::vector<qmomkit::MomentSet> changed_moments{
        {2.0, 5.0, 14.0, 41.0}, // weights 0.5 and 1.5 at abscissas 1 and 3
        {0.0, 0.0, 0.0, 0.0},   // nothing left
        {1.0, -1.0, 2.0, -4.0}, // a negative mean on [0, inf)
    };

    for (const qmomkit::MomentSet& moments : changed_moments)
    {
        qmomkit::QmomStateCandidate candidate = qmom_case.state_candidate(moments, nonnegative);
        if (candidate.ok())
        {
            state = std::move(*candidate.state);
            std::cout << "accepted, new state:\n";
            const qmomkit::QmomResult& quadrature = state.quadrature();
            for (std::size_t node = 0; node < state.active_node_count(); ++node)
            {
                std::cout << "  weight " << quadrature.weights[node] << " at abscissa " << quadrature.abscissas[node] << '\n';
            }
        }
        else if (candidate.status == qmomkit::StateCandidateStatus::empty)
        {
            std::cout << "empty: nothing to evolve, the state is kept\n";
        }
        else
        {
            std::cout << "not accepted, the state is kept: " << candidate.message << '\n';
        }
    }

    return 0;
}
