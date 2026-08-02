#include "monte_carlo/moments.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace
{

using qmomkit::RawMomentRequest;
using qmomkit::monte_carlo::Ensemble;
using qmomkit::monte_carlo::EnsembleMomentStatus;
using qmomkit::monte_carlo::SelectedEnsembleMomentEstimate;

void check_empty_failure(const SelectedEnsembleMomentEstimate& result)
{
    CHECK(!result.ok());
    CHECK(result.indices.empty());
    CHECK(result.values.empty());
    CHECK(result.standard_errors.empty());
    CHECK(!result.message.empty());
}

void check_selected_matches_dense(const Ensemble& ensemble, const double number_density, const RawMomentRequest& request,
                                  const std::vector<std::size_t>& maximum_orders)
{
    const SelectedEnsembleMomentEstimate selected = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, request);
    const qmomkit::monte_carlo::TensorMomentEstimate dense = qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_orders);

    CHECK(selected.ok());
    CHECK(dense.ok());
    CHECK(selected.indices == request.indices);
    CHECK(selected.values.size() == request.indices.size());
    CHECK(selected.standard_errors.size() == request.indices.size());

    for (std::size_t requested = 0; requested < request.indices.size(); ++requested)
    {
        const std::size_t flat = qmomkit::flatten_multi_index(request.indices[requested], dense.moments.extents);
        const double value_tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                       std::max({1.0, std::abs(selected.values[requested]), std::abs(dense.moments.values[flat])});
        const double error_tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                       std::max({1.0, std::abs(selected.standard_errors[requested]), std::abs(dense.standard_errors[flat])});
        CHECK_CLOSE(selected.values[requested], dense.moments.values[flat], value_tolerance);
        CHECK_CLOSE(selected.standard_errors[requested], dense.standard_errors[flat], error_tolerance);
    }
}

void test_selected_one_dimensional_moments()
{
    const Ensemble ensemble{{-2.0}, {3.0}, {0.0}};
    const RawMomentRequest request{{{2}, {0}, {1}, {2}}};
    const SelectedEnsembleMomentEstimate result = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 2.5, request);

    CHECK(result.ok());
    CHECK(result.indices == request.indices);
    CHECK(result.values.size() == request.indices.size());
    CHECK(result.standard_errors.size() == request.indices.size());
    CHECK_CLOSE(result.values[0], 65.0 / 6.0, 1.0e-14);
    CHECK_CLOSE(result.values[1], 2.5, 0.0);
    CHECK_CLOSE(result.values[2], 5.0 / 6.0, 1.0e-14);
    CHECK_CLOSE(result.values[3], result.values[0], 0.0);
    CHECK_CLOSE(result.standard_errors[0], result.standard_errors[3], 0.0);
}

void test_selected_multidimensional_moments_match_dense()
{
    const Ensemble ensemble{{-2.0, 1.0}, {3.0, -4.0}};
    const RawMomentRequest request{{{2, 1}, {0, 0}, {1, 2}, {2, 1}, {1, 0}}};

    check_selected_matches_dense(ensemble, 4.0, request, {2, 2});

    const SelectedEnsembleMomentEstimate result = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 4.0, request);
    CHECK_CLOSE(result.values[1], 4.0, 0.0);
    CHECK_CLOSE(result.values[3], result.values[0], 0.0);
}

void test_selected_zero_density_and_integer_power_edges()
{
    const Ensemble ensemble{{0.0, -2.0}, {0.0, 3.0}};
    const RawMomentRequest request{{{0, 0}, {0, 1}, {0, 2}}};
    const SelectedEnsembleMomentEstimate result = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 2.0, request);

    CHECK(result.ok());
    CHECK_CLOSE(result.values[0], 2.0, 0.0);
    CHECK_CLOSE(result.values[1], 1.0, 0.0);
    CHECK_CLOSE(result.values[2], 13.0, 0.0);

    const RawMomentRequest high_order{{{10}}};
    const SelectedEnsembleMomentEstimate high_order_result =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{2.0}}, 1.0, high_order);
    CHECK(high_order_result.ok());
    CHECK(high_order_result.standard_errors.size() == 1);
    CHECK_CLOSE(high_order_result.values.front(), 1024.0, 0.0);
    CHECK_CLOSE(high_order_result.standard_errors.front(), 0.0, 0.0);

    const SelectedEnsembleMomentEstimate zero_density = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 0.0, request);
    CHECK(zero_density.ok());
    for (const double value : zero_density.values)
    {
        CHECK_CLOSE(value, 0.0, 0.0);
    }
    for (const double error : zero_density.standard_errors)
    {
        CHECK_CLOSE(error, 0.0, 0.0);
    }
}

void test_empty_request_requires_valid_ensemble_and_density()
{
    const RawMomentRequest request;
    const SelectedEnsembleMomentEstimate empty_request = qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, 3.0, request);
    CHECK(empty_request.ok());
    CHECK(empty_request.indices.empty());
    CHECK(empty_request.values.empty());
    CHECK(empty_request.standard_errors.empty());

    const SelectedEnsembleMomentEstimate empty_ensemble = qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{}, 3.0, request);
    CHECK(empty_ensemble.status == EnsembleMomentStatus::invalid_input);
    CHECK(empty_ensemble.message == "Ensemble must contain at least one sample.");

    const SelectedEnsembleMomentEstimate invalid_density =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, -1.0, request);
    check_empty_failure(invalid_density);
    CHECK(invalid_density.status == EnsembleMomentStatus::invalid_input);
    CHECK(invalid_density.message == "Number density must be finite and non-negative.");

    const SelectedEnsembleMomentEstimate nan_density =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, std::numeric_limits<double>::quiet_NaN(), request);
    check_empty_failure(nan_density);
    CHECK(nan_density.status == EnsembleMomentStatus::invalid_input);
    CHECK(nan_density.message == "Number density must be finite and non-negative.");

    const SelectedEnsembleMomentEstimate infinite_density =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, std::numeric_limits<double>::infinity(), request);
    check_empty_failure(infinite_density);
    CHECK(infinite_density.status == EnsembleMomentStatus::invalid_input);
    CHECK(infinite_density.message == "Number density must be finite and non-negative.");
}

void test_selected_invalid_inputs_have_no_partial_payload()
{
    const RawMomentRequest valid_request{{{0, 0}}};

    const SelectedEnsembleMomentEstimate zero_dimension = qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{}}, 1.0, valid_request);
    check_empty_failure(zero_dimension);
    CHECK(zero_dimension.status == EnsembleMomentStatus::invalid_input);
    CHECK(zero_dimension.message == "Samples must contain at least one phase-space coordinate.");

    const SelectedEnsembleMomentEstimate inconsistent =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0}, {2.0, 3.0}}, 1.0, RawMomentRequest{{{0}}});
    check_empty_failure(inconsistent);
    CHECK(inconsistent.status == EnsembleMomentStatus::invalid_input);
    CHECK(inconsistent.message == "All samples must have the same phase-space dimension.");

    const SelectedEnsembleMomentEstimate non_finite_before_inconsistent = qmomkit::monte_carlo::compute_requested_ensemble_moments(
        Ensemble{{std::numeric_limits<double>::quiet_NaN()}, {1.0, 2.0}}, 1.0, RawMomentRequest{{{0}}});
    check_empty_failure(non_finite_before_inconsistent);
    CHECK(non_finite_before_inconsistent.status == EnsembleMomentStatus::invalid_input);
    CHECK(non_finite_before_inconsistent.message == "All samples must have the same phase-space dimension.");

    const SelectedEnsembleMomentEstimate rank_mismatch =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, 1.0, RawMomentRequest{{{0}}});
    check_empty_failure(rank_mismatch);
    CHECK(rank_mismatch.status == EnsembleMomentStatus::invalid_input);
    CHECK(rank_mismatch.message == "Requested raw moment rank must match the sample phase-space dimension.");

    const SelectedEnsembleMomentEstimate non_finite_sample =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, std::numeric_limits<double>::quiet_NaN()}}, 1.0, valid_request);
    check_empty_failure(non_finite_sample);
    CHECK(non_finite_sample.status == EnsembleMomentStatus::non_finite_sample);
    CHECK(non_finite_sample.message == "Sample coordinate is non-finite.");

    const SelectedEnsembleMomentEstimate nan_density =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, std::numeric_limits<double>::quiet_NaN(), valid_request);
    check_empty_failure(nan_density);
    CHECK(nan_density.status == EnsembleMomentStatus::invalid_input);

    const SelectedEnsembleMomentEstimate infinite_density =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, std::numeric_limits<double>::infinity(), valid_request);
    check_empty_failure(infinite_density);
    CHECK(infinite_density.status == EnsembleMomentStatus::invalid_input);

    const RawMomentRequest oversized_request{{{0, std::numeric_limits<std::size_t>::max()}}};
    const SelectedEnsembleMomentEstimate oversized =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0, 2.0}}, 1.0, oversized_request);
    check_empty_failure(oversized);
    CHECK(oversized.status == EnsembleMomentStatus::invalid_input);
    CHECK(oversized.message == "Requested raw moment exponent is too large.");
}

void test_selected_numerical_failures()
{
    const RawMomentRequest monomial_overflow{{{2}}};
    const SelectedEnsembleMomentEstimate monomial_result =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{std::numeric_limits<double>::max()}}, 1.0, monomial_overflow);
    check_empty_failure(monomial_result);
    CHECK(monomial_result.status == EnsembleMomentStatus::non_finite_moment);
    CHECK(monomial_result.message == "Tensor moment accumulation became non-finite.");

    const RawMomentRequest sum_squares_overflow{{{1}}};
    const SelectedEnsembleMomentEstimate sum_squares_result =
        qmomkit::monte_carlo::compute_requested_ensemble_moments(Ensemble{{1.0e200}, {1.0e200}}, 1.0, sum_squares_overflow);
    check_empty_failure(sum_squares_result);
    CHECK(sum_squares_result.status == EnsembleMomentStatus::non_finite_moment);
}

void test_selected_inputs_are_unchanged()
{
    const Ensemble ensemble{{-1.0, 2.0}, {3.0, 4.0}};
    const RawMomentRequest request{{{1, 0}, {0, 1}, {1, 0}}};
    const Ensemble original_ensemble = ensemble;
    const RawMomentRequest original_request = request;

    const SelectedEnsembleMomentEstimate result = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, 2.0, request);
    CHECK(result.ok());
    CHECK(ensemble == original_ensemble);
    CHECK(request.indices == original_request.indices);
}

void test_failed_selection_does_not_mutate_inputs()
{
    const Ensemble ensemble{{-1.0, 2.0}, {3.0, 4.0}};
    const RawMomentRequest request{{{1, 0}, {0, 1}, {1, 0}}};
    const Ensemble original_ensemble = ensemble;
    const RawMomentRequest original_request = request;

    const SelectedEnsembleMomentEstimate result = qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, -1.0, request);
    check_empty_failure(result);
    CHECK(ensemble == original_ensemble);
    CHECK(request.indices == original_request.indices);
}

} // namespace

int main()
{
    test_selected_one_dimensional_moments();
    test_selected_multidimensional_moments_match_dense();
    test_selected_zero_density_and_integer_power_edges();
    test_empty_request_requires_valid_ensemble_and_density();
    test_selected_invalid_inputs_have_no_partial_payload();
    test_selected_numerical_failures();
    test_selected_inputs_are_unchanged();
    test_failed_selection_does_not_mutate_inputs();
    return 0;
}
