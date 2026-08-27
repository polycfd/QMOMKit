#include <cstddef>
#include <monte_carlo/birth.hpp>
#include <monte_carlo/moments.hpp>

int main()
{
    const qmomkit::monte_carlo::Ensemble ensemble{{-1.0, 2.0}, {3.0, 4.0}};
    const qmomkit::RawMomentRequest request{{{1, 0}, {0, 1}, {1, 0}}};
    const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 2.0, request);

    if (!estimate.ok() || estimate.indices != request.indices || estimate.values.size() != request.indices.size() ||
        estimate.standard_errors.size() != request.indices.size())
    {
        return 1;
    }

    const qmomkit::monte_carlo::PopulationBirthResult no_birth = qmomkit::monte_carlo::inject_population(ensemble, 2.0, 0.0, {}, 1701);
    if (!no_birth.ok() || no_birth.ensemble != ensemble || no_birth.number_density != 2.0)
    {
        return 1;
    }

    return 0;
}
