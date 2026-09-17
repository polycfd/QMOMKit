#include "monte_carlo/birth.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_case_dispatch.hpp"
#include "test_checks.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace
{

using qmomkit::monte_carlo::BirthCoordinateSampler;
using qmomkit::monte_carlo::Ensemble;
using qmomkit::monte_carlo::PopulationBirthResult;
using qmomkit::monte_carlo::Sample;
using qmomkit::monte_carlo::SamplingStatus;

void test_zero_birth_is_a_no_op()
{
    const Ensemble ensemble{{1.0}, {3.0}};
    std::size_t sampler_calls = 0;
    const BirthCoordinateSampler sampler = [&sampler_calls](std::mt19937_64&)
    {
        ++sampler_calls;
        return Sample{5.0};
    };

    const PopulationBirthResult result = qmomkit::monte_carlo::inject_population(ensemble, 4.0, 0.0, sampler, 1201);

    CHECK(result.ok());
    CHECK(result.ensemble == ensemble);
    CHECK(result.number_density == 4.0);
    CHECK(result.seed == 1201);
    CHECK(sampler_calls == 0);
}

void test_deterministic_birth_preserves_weight_and_moments()
{
    const Ensemble ensemble{{1.0}, {3.0}};
    constexpr double number_density = 4.0;
    constexpr double representative_weight = 2.0;
    const BirthCoordinateSampler sampler = [](std::mt19937_64&) { return Sample{5.0}; };

    const PopulationBirthResult result = qmomkit::monte_carlo::inject_population(ensemble, number_density, 40.0, sampler, 1202);

    CHECK(result.ok());
    CHECK(result.ensemble.size() > ensemble.size());
    const std::size_t birth_count = result.ensemble.size() - ensemble.size();
    CHECK(result.ensemble[0] == ensemble[0]);
    CHECK(result.ensemble[1] == ensemble[1]);
    for (std::size_t birth = 0; birth < birth_count; ++birth)
    {
        CHECK(result.ensemble[ensemble.size() + birth] == Sample{5.0});
    }
    CHECK_CLOSE(result.number_density, representative_weight * static_cast<double>(result.ensemble.size()), 0.0);

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate moments =
        qmomkit::monte_carlo::compute_ensemble_moments(result.ensemble, result.number_density, 2);
    CHECK(moments.ok());
    CHECK_CLOSE(moments.moments[0], representative_weight * static_cast<double>(ensemble.size() + birth_count), 0.0);
    CHECK_CLOSE(moments.moments[1], representative_weight * (1.0 + 3.0 + 5.0 * static_cast<double>(birth_count)), 0.0);
    CHECK_CLOSE(moments.moments[2], representative_weight * (1.0 + 9.0 + 25.0 * static_cast<double>(birth_count)), 0.0);
}

void test_multidimensional_birth_preserves_phase_space_order()
{
    const Ensemble ensemble{{1.0, 2.0}, {3.0, 4.0}};
    const BirthCoordinateSampler sampler = [](std::mt19937_64&) { return Sample{5.0, 6.0}; };

    const PopulationBirthResult result = qmomkit::monte_carlo::inject_population(ensemble, 2.0, 24.0, sampler, 1203);

    CHECK(result.ok());
    CHECK(result.ensemble.size() > ensemble.size());
    const std::size_t birth_count = result.ensemble.size() - ensemble.size();
    const Sample born_coordinate{5.0, 6.0};
    for (std::size_t birth = 0; birth < birth_count; ++birth)
    {
        CHECK(result.ensemble[ensemble.size() + birth] == born_coordinate);
    }

    const qmomkit::monte_carlo::TensorMomentEstimate moments =
        qmomkit::monte_carlo::compute_ensemble_moments(result.ensemble, result.number_density, std::vector<std::size_t>{1, 1});
    CHECK(moments.ok());
    const std::size_t mixed_index = qmomkit::flatten_multi_index(qmomkit::MultiIndex{1, 1}, moments.moments.extents);
    CHECK_CLOSE(moments.moments.values[mixed_index], 14.0 + 30.0 * static_cast<double>(birth_count), 0.0);
}

void test_seed_reproducibility()
{
    const Ensemble ensemble{{0.5, -1.0}, {1.5, 2.0}};
    const BirthCoordinateSampler sampler = [](std::mt19937_64& engine)
    {
        std::uniform_real_distribution<double> distribution(-2.0, 3.0);
        return Sample{distribution(engine), distribution(engine)};
    };

    const PopulationBirthResult first = qmomkit::monte_carlo::inject_population(ensemble, 3.0, 18.0, sampler, 1204);
    const PopulationBirthResult second = qmomkit::monte_carlo::inject_population(ensemble, 3.0, 18.0, sampler, 1204);

    CHECK(first.ok());
    CHECK(second.ok());
    CHECK(first.ensemble == second.ensemble);
    CHECK(first.number_density == second.number_density);
    CHECK(first.seed == second.seed);
}

void test_invalid_inputs()
{
    const Ensemble ensemble{{1.0}};
    const BirthCoordinateSampler valid_sampler = [](std::mt19937_64&) { return Sample{2.0}; };

    const PopulationBirthResult negative = qmomkit::monte_carlo::inject_population(ensemble, 1.0, -1.0, valid_sampler, 1205);
    CHECK(!negative.ok());
    CHECK(negative.status == SamplingStatus::invalid_input);
    CHECK(negative.ensemble.empty());

    const PopulationBirthResult non_finite =
        qmomkit::monte_carlo::inject_population(ensemble, 1.0, std::numeric_limits<double>::infinity(), valid_sampler, 1205);
    CHECK(!non_finite.ok());
    CHECK(non_finite.status == SamplingStatus::invalid_input);
    CHECK(non_finite.ensemble.empty());

    const PopulationBirthResult nan_increment =
        qmomkit::monte_carlo::inject_population(ensemble, 1.0, std::numeric_limits<double>::quiet_NaN(), valid_sampler, 1205);
    CHECK(!nan_increment.ok());
    CHECK(nan_increment.status == SamplingStatus::invalid_input);
    CHECK(nan_increment.ensemble.empty());

    const PopulationBirthResult empty = qmomkit::monte_carlo::inject_population(Ensemble{}, 1.0, 10.0, valid_sampler, 1205);
    CHECK(!empty.ok());
    CHECK(empty.status == SamplingStatus::invalid_input);
    CHECK(empty.ensemble.empty());

    const PopulationBirthResult zero_density = qmomkit::monte_carlo::inject_population(ensemble, 0.0, 10.0, valid_sampler, 1205);
    CHECK(!zero_density.ok());
    CHECK(zero_density.status == SamplingStatus::invalid_input);
    CHECK(zero_density.ensemble.empty());

    const Ensemble inconsistent{{1.0}, {2.0, 3.0}};
    const PopulationBirthResult inconsistent_result = qmomkit::monte_carlo::inject_population(inconsistent, 2.0, 10.0, valid_sampler, 1205);
    CHECK(!inconsistent_result.ok());
    CHECK(inconsistent_result.status == SamplingStatus::invalid_input);
    CHECK(inconsistent_result.ensemble.empty());

    const PopulationBirthResult reference = qmomkit::monte_carlo::inject_population(ensemble, 1.0, 100.0, valid_sampler, 1206);
    CHECK(reference.ok());
    CHECK(reference.ensemble.size() > ensemble.size());

    const BirthCoordinateSampler wrong_rank = [](std::mt19937_64&) { return Sample{2.0, 3.0}; };
    const PopulationBirthResult wrong_rank_result = qmomkit::monte_carlo::inject_population(ensemble, 1.0, 100.0, wrong_rank, 1206);
    CHECK(!wrong_rank_result.ok());
    CHECK(wrong_rank_result.status == SamplingStatus::invalid_input);
    CHECK(wrong_rank_result.ensemble.empty());

    const BirthCoordinateSampler non_finite_coordinate = [](std::mt19937_64&) { return Sample{std::numeric_limits<double>::quiet_NaN()}; };
    const PopulationBirthResult non_finite_coordinate_result =
        qmomkit::monte_carlo::inject_population(ensemble, 1.0, 100.0, non_finite_coordinate, 1206);
    CHECK(!non_finite_coordinate_result.ok());
    CHECK(non_finite_coordinate_result.status == SamplingStatus::non_finite_sample);
    CHECK(non_finite_coordinate_result.ensemble.empty());
}

void test_injected_samples_are_evolved()
{
    const Ensemble ensemble{{0.0}};
    const BirthCoordinateSampler sampler = [](std::mt19937_64&) { return Sample{2.0}; };
    const PopulationBirthResult birth = qmomkit::monte_carlo::inject_population(ensemble, 1.0, 20.0, sampler, 1207);
    CHECK(birth.ok());
    CHECK(birth.ensemble.size() > ensemble.size());

    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{1.0}; });
    const qmomkit::monte_carlo::EnsembleEvolutionResult evolved = qmomkit::monte_carlo::evolve_ensemble(birth.ensemble, velocity, 0.0, 0.5);
    CHECK(evolved.ok());
    CHECK(evolved.ensemble.size() == birth.ensemble.size());
    CHECK_CLOSE(evolved.ensemble.front()[0], 0.5, 1.0e-12);
    for (std::size_t sample = 1; sample < evolved.ensemble.size(); ++sample)
    {
        CHECK_CLOSE(evolved.ensemble[sample][0], 2.5, 1.0e-12);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::monte_carlo::test::dispatch_test_case(
        argc, argv,
        {
            {"zero_birth_is_a_no_op", test_zero_birth_is_a_no_op},
            {"deterministic_birth_preserves_weight_and_moments", test_deterministic_birth_preserves_weight_and_moments},
            {"multidimensional_birth_preserves_phase_space_order", test_multidimensional_birth_preserves_phase_space_order},
            {"seed_reproducibility", test_seed_reproducibility},
            {"invalid_inputs", test_invalid_inputs},
            {"injected_samples_are_evolved", test_injected_samples_are_evolved},
        });
}
