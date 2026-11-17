#include <cstddef>
#include <monte_carlo/moments.hpp>

// A shared library of a host that calls into the installed static companion. On
// Linux it links only if the companion was compiled as position-independent
// code.
std::size_t monte_carlo_package_consumer_shared_library_moment_count()
{
    const qmomkit::monte_carlo::Ensemble ensemble{{-1.0, 2.0}, {3.0, 4.0}};
    const qmomkit::RawMomentRequest request{{{1, 0}, {0, 1}}};
    const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 2.0, request);
    return estimate.ok() ? estimate.values.size() : 0;
}
