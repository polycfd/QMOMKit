#include "internal/cqmom_case.hpp"
#include "internal/performance.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_relative_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)));
}

bool contains_index(const qmomkit::MultiIndexTable& indices, const qmomkit::MultiIndexView index)
{
    return std::find(indices.begin(), indices.end(), index) != indices.end();
}

void append_unique(std::vector<qmomkit::MultiIndex>& indices, qmomkit::MultiIndex index)
{
    if (std::find(indices.begin(), indices.end(), index) == indices.end())
    {
        indices.push_back(std::move(index));
    }
}

void assert_indices_equal(const qmomkit::MultiIndexTable& actual, const std::vector<qmomkit::MultiIndex>& expected)
{
    assert(actual.to_multi_indices() == expected);
}

void assert_unique_indices(const qmomkit::MultiIndexTable& indices)
{
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        for (std::size_t j = i + 1; j < indices.size(); ++j)
        {
            assert(indices[i] != indices[j]);
        }
    }
}

template <typename Function> void assert_throws_out_of_range(Function function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::out_of_range&)
    {
        threw = true;
    }
    assert(threw);
}

template <typename Function> void assert_throws_invalid_argument(Function function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
}

std::vector<double> raw_moments_1d(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < abscissas.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], static_cast<int>(order));
        }
    }
    return moments;
}

void test_moment_set_public_container_behavior()
{
    const qmomkit::MomentSet empty;
    assert(empty.empty());
    assert(empty.size() == 0);
    assert(empty.values().empty());

    const qmomkit::MomentSet moments{1.0, -2.5, 4.25};
    assert(!moments.empty());
    assert(moments.size() == 3);
    assert(moments.values() == std::vector<double>({1.0, -2.5, 4.25}));
    assert_near(moments[0], 1.0, 0.0);
    assert_near(moments[2], 4.25, 0.0);
    assert_throws_out_of_range([&moments]() { (void)moments[3]; });
}

std::vector<double> indexed_values(const std::vector<std::size_t>& extents)
{
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        values[flat] = 10.0 + static_cast<double>(flat);
    }
    return values;
}

std::vector<double> indexed_count_values(const std::size_t count)
{
    std::vector<double> values(count, 0.0);
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        values[index] = 10.0 + static_cast<double>(index);
    }
    return values;
}

std::vector<double> raw_moments_from_grid(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<std::size_t>& extents,
                                          const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
        {
            const std::size_t flat_node = hierarchy.flatten_node_index(node_index);
            double contribution = weights[flat_node];
            for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][flat_node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

std::vector<double> raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                               const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            moments[flat_moment] += contribution;
        }
    }
    return moments;
}

// Build a case from a raw moment tensor by selecting the moments its hierarchy
// reconstructs from. The tests below state their moments as tensors because
// that is how their analytical references are computed; the parameters mirror
// the case fields so that each call reads as the case it builds. Unlike the
// case, the helper is fixed unless the QMOM options ask for adaptive
// reconstruction.
qmomkit::CqmomCase cqmom_case_from_tensor(const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::RawMomentTensor& moments,
                                          const qmomkit::ConditionalMomentOptions& conditional_options = {},
                                          const qmomkit::QmomOptions& qmom_options = qmomkit::QmomOptions{false},
                                          const qmomkit::CqmomDiagnosticLevel diagnostic_level = qmomkit::CqmomDiagnosticLevel::none,
                                          const qmomkit::VandermondeSolver vandermonde_solver = qmomkit::VandermondeSolver::automatic,
                                          qmomkit::CqmomAbscissaBounds coordinate_support = {})
{
    return {qmomkit::extract_cqmom_moments(moments, qmomkit::create_cqmom_moment_layout(hierarchy)),
            conditional_options,
            qmom_options,
            diagnostic_level,
            vandermonde_solver,
            std::move(coordinate_support)};
}

void assert_same_quadrature(const qmomkit::CqmomQuadrature& lhs, const qmomkit::CqmomQuadrature& rhs, const double tolerance)
{
    assert(lhs.number_of_dimensions() == rhs.number_of_dimensions());
    assert(lhs.total_nodes() == rhs.total_nodes());
    assert(lhs.adaptive_branch_counts_enabled() == rhs.adaptive_branch_counts_enabled());
    assert(lhs.has_rectangular_storage() == rhs.has_rectangular_storage());
    assert(lhs.terminal_paths() == rhs.terminal_paths());
    for (std::size_t node = 0; node < lhs.weights().size(); ++node)
    {
        assert_near(lhs.weights()[node], rhs.weights()[node], tolerance);
    }
    for (std::size_t dimension = 0; dimension < lhs.number_of_dimensions(); ++dimension)
    {
        for (std::size_t node = 0; node < lhs.total_nodes(); ++node)
        {
            assert_near(lhs.abscissas()[dimension][node], rhs.abscissas()[dimension][node], tolerance);
        }
    }
}

qmomkit::CqmomMomentData moment_data_from_tensor(const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::RawMomentTensor& moment_tensor)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    return qmomkit::extract_cqmom_moments(moment_tensor, layout);
}

qmomkit::RawMomentTensor public_cqmom_promotion_moments()
{
    return {{4, 4},
            {
                1.0,
                0.44,
                0.424,
                0.5552,
                0.55,
                0.242,
                0.2332,
                0.30536,
                1.2125,
                0.5335,
                0.5141,
                0.67318,
                1.121875,
                0.493625,
                0.475675,
                0.622865,
            }};
}

// Assigning through a reference keeps the self-assignment cases out of reach of
// -Wself-assign while still calling the same operator the caller would.
template <class T> void assign_from(T& target, const T& source)
{
    target = source;
}

qmomkit::RawMomentTensor scaled_moments(qmomkit::RawMomentTensor moments, const double factor)
{
    // Scaling every raw moment by the same factor scales the quadrature weights
    // and leaves the abscissas alone, so the result stays realizable.
    for (double& value : moments.values)
    {
        value *= factor;
    }
    return moments;
}

// CqmomPreparedCase and CqmomState are pimpl types whose copy operations are written by
// hand: the payload variant has to be deep-copied while the prepared metadata is
// deliberately shared. Copy assignment is also how every CQMOM example advances
// its loop (`state = result.state`), and none of the four assignment operators
// had a caller anywhere in the library, the tests or the examples.
void test_cqmom_prepared_case_and_state_assignment_preserve_independent_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor first_moments = public_cqmom_promotion_moments();
    const qmomkit::RawMomentTensor second_moments = scaled_moments(first_moments, 2.0);
    assert(first_moments.values != second_moments.values);

    const qmomkit::CqmomPreparedCase first_cqmom_prepared_case{cqmom_case_from_tensor(hierarchy, first_moments)};
    const qmomkit::CqmomPreparedCase second_cqmom_prepared_case{cqmom_case_from_tensor(hierarchy, second_moments)};

    // `initial_state()` stores the canonical moments derived from the accepted
    // quadrature rather than echoing the input tensor, so compare states to states.
    const std::vector<double> first_values = first_cqmom_prepared_case.initial_state().moment_values();
    const std::vector<double> second_values = second_cqmom_prepared_case.initial_state().moment_values();
    assert(first_values != second_values);

    // Copy assignment replaces the target's payload, leaves the source alone, and
    // keeps sharing the source's immutable prepared metadata.
    qmomkit::CqmomPreparedCase assigned_cqmom_prepared_case = first_cqmom_prepared_case;
    assert(assigned_cqmom_prepared_case.initial_state().moment_values() == first_values);
    assigned_cqmom_prepared_case = second_cqmom_prepared_case;
    assert(assigned_cqmom_prepared_case.initial_state().moment_values() == second_values);
    assert(second_cqmom_prepared_case.initial_state().moment_values() == second_values);
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(assigned_cqmom_prepared_case).get() ==
           qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(second_cqmom_prepared_case).get());

    // The self-assignment guard: without it the payload is released before it is
    // read back, so the case has to survive this and stay usable afterwards.
    assign_from(assigned_cqmom_prepared_case, assigned_cqmom_prepared_case);
    assert(assigned_cqmom_prepared_case.initial_state().moment_values() == second_values);
    const qmomkit::CqmomResult after_self_assignment = qmomkit::compute_cqmom(assigned_cqmom_prepared_case);
    assert(after_self_assignment.ok());

    qmomkit::CqmomPreparedCase movable_cqmom_prepared_case = first_cqmom_prepared_case;
    qmomkit::CqmomPreparedCase move_assigned_cqmom_prepared_case = second_cqmom_prepared_case;
    move_assigned_cqmom_prepared_case = std::move(movable_cqmom_prepared_case);
    assert(move_assigned_cqmom_prepared_case.initial_state().moment_values() == first_values);

    // A copy is independent of later edits to the original, which is what
    // `config::with_initial_moments` returning a new case relies on.
    const qmomkit::CqmomPreparedCase copied_cqmom_prepared_case = first_cqmom_prepared_case;
    const auto replaced = qmomkit::config::with_initial_moments(first_cqmom_prepared_case, second_moments);
    assert(replaced.ok());
    assert(replaced.value().initial_state().moment_values() == second_values);
    assert(copied_cqmom_prepared_case.initial_state().moment_values() == first_values);
    assert(first_cqmom_prepared_case.initial_state().moment_values() == first_values);

    // The evolution-loop pattern every CQMOM example uses: snapshot the state,
    // advance, then assign the result back over the live state.
    // The buffer is packed over every terminal node, so fill all of it: a
    // partially filled buffer is a contract violation, not a small drift.
    const qmomkit::PhaseSpaceVelocity velocity = [](const double, const std::vector<double>&, std::vector<double>& velocities, const std::size_t)
    { std::fill(velocities.begin(), velocities.end(), 0.5); };

    qmomkit::CqmomState state = first_cqmom_prepared_case.initial_state();
    const std::vector<double> initial_values = state.moment_values();
    const qmomkit::CqmomState snapshot = state;

    const qmomkit::CqmomPhaseSpaceEvolutionResult advanced = qmomkit::advance_cqmom_advection(state, 0.0, 0.1, velocity);
    assert(advanced.ok());
    state = advanced.state;
    assert(state.moment_values() != initial_values);
    assert(advanced.state.moment_values() == state.moment_values());
    // The snapshot is the rollback target, so it must not follow the live state.
    assert(snapshot.moment_values() == initial_values);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(snapshot).get() ==
           qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(state).get());

    const qmomkit::CqmomPhaseSpaceEvolutionResult second_step = qmomkit::advance_cqmom_advection(state, 0.1, 0.1, velocity);
    assert(second_step.ok());
    const std::vector<double> after_first_step = state.moment_values();
    state = second_step.state;
    assert(state.moment_values() != after_first_step);
    assert(snapshot.moment_values() == initial_values);

    const std::vector<double> before_self_assignment = state.moment_values();
    assign_from(state, state);
    assert(state.moment_values() == before_self_assignment);
    const qmomkit::CqmomPhaseSpaceEvolutionResult after_self_step = qmomkit::advance_cqmom_advection(state, 0.2, 0.1, velocity);
    assert(after_self_step.ok());

    qmomkit::CqmomState movable_state = first_cqmom_prepared_case.initial_state();
    qmomkit::CqmomState move_assigned_state = state;
    move_assigned_state = std::move(movable_state);
    assert(move_assigned_state.moment_values() == initial_values);
    assert(state.moment_values() == before_self_assignment);
}

void test_cqmom_case_promotes_to_state_and_exposes_quadrature()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moment_tensor = public_cqmom_promotion_moments();
    qmomkit::CqmomCase cqmom_case{moment_data_from_tensor(hierarchy, moment_tensor)};
    const qmomkit::CqmomResult expected = qmomkit::compute_cqmom(cqmom_case);
    assert(expected.ok());

    const qmomkit::CqmomPreparedCase cqmom_prepared_case{std::move(cqmom_case)};
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();

    assert(state.hierarchy().node_counts() == hierarchy.node_counts());
    assert(state.moment_indices() == qmomkit::create_cqmom_moment_layout(hierarchy).required_indices.to_multi_indices());
    assert(state.moment_values().size() == state.moment_indices().size());
    assert_same_quadrature(state.quadrature(), expected.quadrature, 1.0e-10);
}

// Moment data of four weighted points on a two-by-two CQMOM hierarchy.
qmomkit::CqmomMomentData grid_moment_data(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<double>& weights,
                                          const std::vector<std::vector<double>>& abscissas)
{
    const std::vector<std::size_t> extents{4, 4};
    return moment_data_from_tensor(hierarchy, {extents, raw_moments_from_grid(hierarchy, extents, weights, abscissas)});
}

// Counts the CQMOM reconstructions of a call in instrumented builds; elsewhere
// the count is not checked.
template <typename Function> void assert_cqmom_reconstructions(const std::size_t expected, Function function)
{
    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        function();
    }
    if (qmomkit::internal::performance::enabled())
    {
        assert(session.snapshot().reconstruction.cqmom_reconstructions == expected);
    }
}

// The state of moment data on a prepared case: the quadrature compute_cqmom
// reconstructs from the data, and the moments of that quadrature.
void assert_cqmom_state_of(const qmomkit::CqmomState& state, const qmomkit::CqmomPreparedCase& prepared_case, const qmomkit::CqmomMomentData& data)
{
    const qmomkit::CqmomResult expected = qmomkit::compute_cqmom(qmomkit::config::with_initial_moments(prepared_case, data).value());
    assert(expected.ok());
    assert_same_quadrature(state.quadrature(), expected.quadrature, 0.0);
    const std::vector<qmomkit::MultiIndex> indices = data.layout.required_indices.to_multi_indices();
    const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(expected.quadrature, qmomkit::RawMomentRequest{indices});
    assert(moments.ok());
    assert(state.moment_indices() == indices);
    assert(state.moment_values() == moments.values);
}

void test_cqmom_state_candidate_matches_initial_state()
{
    // A state from other moments is the quadrature of those moments with the
    // settings of the case, and the moments of that quadrature; it shares the
    // branch structure of the prepared case. The moments of the case give the
    // state of initial_state().
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    for (const bool adaptive : {false, true})
    {
        qmomkit::CqmomCase cqmom_case{grid_moment_data(hierarchy, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 1.5, 1.5}, {0.25, 2.0, 0.25, 2.0}})};
        cqmom_case.qmom_options.adaptive = adaptive;
        const qmomkit::CqmomPreparedCase prepared_case{cqmom_case};
        const auto& prepared_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(prepared_case);

        const qmomkit::CqmomMomentData other = grid_moment_data(hierarchy, {0.3, 0.2, 0.4, 0.1}, {{0.5, 0.5, 1.5, 1.5}, {-0.5, 2.0, 0.5, 1.0}});
        qmomkit::CqmomStateCandidate from_other;
        assert_cqmom_reconstructions(1, [&] { from_other = prepared_case.state_candidate(other); });
        assert(from_other.ok());
        assert(!from_other.realizability.has_value());
        assert(!from_other.rejected_coordinate.has_value());
        assert(from_other.reconstruction_status == qmomkit::CqmomStatus::success);
        assert(from_other.message.empty());
        assert_cqmom_state_of(*from_other.state, prepared_case, other);
        assert(from_other.state->qmom_options().adaptive == adaptive);
        assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(*from_other.state).get() == prepared_metadata.get());

        const qmomkit::CqmomStateCandidate own = prepared_case.state_candidate(cqmom_case.moments);
        const qmomkit::CqmomState initial = prepared_case.initial_state();
        assert(own.state->moment_values() == initial.moment_values());
        assert_same_quadrature(own.state->quadrature(), initial.quadrature(), 0.0);
    }
}

void test_cqmom_state_candidate_applies_coordinate_policies()
{
    // The second coordinate takes the value -0.5 at one point, so only a
    // nonnegative check of that coordinate rejects the moments, before any
    // reconstruction; a policy that is not applied, or one for the first
    // coordinate, accepts them, and no verdict is kept on success.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomPreparedCase prepared_case{
        qmomkit::CqmomCase{grid_moment_data(hierarchy, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 1.5, 1.5}, {0.25, 2.0, 0.25, 2.0}})}};
    const qmomkit::CqmomMomentData moments = grid_moment_data(hierarchy, {0.3, 0.2, 0.4, 0.1}, {{0.5, 0.5, 1.5, 1.5}, {-0.5, 2.0, 0.5, 1.0}});
    qmomkit::RealizabilityPolicy nonnegative;
    nonnegative.support = qmomkit::MomentSupport::nonnegative;
    assert(prepared_case.state_candidate(moments, {std::nullopt, nonnegative}).ok());
    nonnegative.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    const qmomkit::CqmomStateCandidate first_only = prepared_case.state_candidate(moments, {nonnegative});
    assert(first_only.ok());
    assert(!first_only.realizability.has_value());
    assert(!first_only.rejected_coordinate.has_value());

    for (const qmomkit::RealizabilityCheckMode mode : {qmomkit::RealizabilityCheckMode::FinalOnly, qmomkit::RealizabilityCheckMode::EveryStage})
    {
        nonnegative.check_mode = mode;
        qmomkit::CqmomStateCandidate rejected;
        assert_cqmom_reconstructions(0, [&] { rejected = prepared_case.state_candidate(moments, {std::nullopt, nonnegative}); });
        assert(rejected.status == qmomkit::StateCandidateStatus::realizability_rejected);
        assert(!rejected.state.has_value());
        assert(rejected.rejected_coordinate == std::size_t{1});
        assert(rejected.realizability->status == qmomkit::RealizabilityStatus::not_realizable);
        assert(!rejected.reconstruction_status.has_value());
        assert(rejected.message == "Coordinate 1: " + rejected.realizability->message);
    }

    // With the first coordinate at 0 for two points, its marginal is on the
    // boundary of the nonnegative support: accepted unless boundary moments
    // are excluded. An interval is read only for a bounded-interval support.
    const qmomkit::CqmomMomentData at_origin = grid_moment_data(hierarchy, {0.3, 0.2, 0.4, 0.1}, {{0.0, 0.0, 1.5, 1.5}, {0.5, 2.0, 0.5, 1.0}});
    nonnegative.interval = qmomkit::BoundedInterval{0.0, 1.0};
    assert(prepared_case.state_candidate(at_origin, {nonnegative, nonnegative}).ok());
    nonnegative.accept_boundary = false;
    const qmomkit::CqmomStateCandidate strict_only = prepared_case.state_candidate(at_origin, {nonnegative});
    assert(strict_only.status == qmomkit::StateCandidateStatus::realizability_rejected);
    assert(strict_only.rejected_coordinate == std::size_t{0});
    assert(strict_only.realizability->status == qmomkit::RealizabilityStatus::boundary);

    // A bounded interval is checked against its bounds, and is invalid input
    // without them; the tolerance of a policy applies.
    qmomkit::RealizabilityPolicy interval;
    interval.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    interval.support = qmomkit::MomentSupport::bounded_interval;
    const qmomkit::CqmomStateCandidate no_bounds = prepared_case.state_candidate(at_origin, {std::nullopt, interval});
    assert(no_bounds.status == qmomkit::StateCandidateStatus::invalid_input);
    assert(no_bounds.message == "CQMOM coordinate 1 realizability policy: A bounded-interval realizability policy requires an interval.");
    assert(!no_bounds.rejected_coordinate.has_value());
    assert(!no_bounds.realizability.has_value());
    interval.interval = qmomkit::BoundedInterval{0.0, 1.5};
    assert(prepared_case.state_candidate(at_origin, {std::nullopt, interval}).status == qmomkit::StateCandidateStatus::realizability_rejected);
    interval.interval = qmomkit::BoundedInterval{0.0, 2.5};
    assert(prepared_case.state_candidate(at_origin, {std::nullopt, interval}).ok());
    qmomkit::RealizabilityPolicy coarse;
    coarse.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    assert(prepared_case.state_candidate(moments, {coarse}).ok());
    coarse.options.tolerance = 0.9;
    const qmomkit::CqmomStateCandidate indeterminate = prepared_case.state_candidate(moments, {coarse});
    assert(indeterminate.status == qmomkit::StateCandidateStatus::realizability_rejected);
    assert(indeterminate.realizability->status == qmomkit::RealizabilityStatus::numerically_indeterminate);
    coarse.accept_indeterminate = true;
    assert(prepared_case.state_candidate(moments, {coarse}).ok());

    // A policy beyond the hierarchy is invalid input, whatever its mode and
    // before any coordinate is checked, as in the moment-source update.
    const std::string beyond_message = "CQMOM coordinate realizability policy index is outside the phase-space dimension.";
    for (const qmomkit::RealizabilityCheckMode mode : {qmomkit::RealizabilityCheckMode::None, qmomkit::RealizabilityCheckMode::FinalOnly})
    {
        qmomkit::RealizabilityPolicy beyond;
        beyond.check_mode = mode;
        const qmomkit::CqmomStateCandidate candidate = prepared_case.state_candidate(moments, {std::nullopt, nonnegative, beyond});
        assert(candidate.status == qmomkit::StateCandidateStatus::invalid_input);
        assert(candidate.message == beyond_message);
    }

    // An invalid policy is invalid input before any coordinate is checked,
    // even when an earlier coordinate would reject its moments.
    qmomkit::RealizabilityPolicy unknown_mode;
    unknown_mode.check_mode = static_cast<qmomkit::RealizabilityCheckMode>(7);
    nonnegative.accept_boundary = false;
    assert(prepared_case.state_candidate(at_origin, {nonnegative}).status == qmomkit::StateCandidateStatus::realizability_rejected);
    const qmomkit::CqmomStateCandidate unknown = prepared_case.state_candidate(at_origin, {nonnegative, unknown_mode});
    assert(unknown.status == qmomkit::StateCandidateStatus::invalid_input);
    assert(unknown.message == "CQMOM coordinate 1 realizability policy: Unknown realizability check mode.");
}

void test_cqmom_state_candidate_reports_failures_without_throwing()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    qmomkit::CqmomCase cqmom_case{grid_moment_data(hierarchy, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 1.5, 1.5}, {0.25, 2.0, 0.25, 2.0}})};
    cqmom_case.qmom_options.adaptive = false;
    const qmomkit::CqmomPreparedCase prepared_case{cqmom_case};
    const qmomkit::CqmomMomentData moments = grid_moment_data(hierarchy, {0.3, 0.2, 0.4, 0.1}, {{0.5, 0.5, 1.5, 1.5}, {-0.5, 2.0, 0.5, 1.0}});
    const auto assert_failed = [](const qmomkit::CqmomStateCandidate& candidate, const qmomkit::StateCandidateStatus status,
                                  const std::optional<qmomkit::CqmomStatus> reconstruction_status)
    {
        assert(candidate.status == status);
        assert(!candidate.ok());
        assert(!candidate.state.has_value());
        assert(candidate.reconstruction_status == reconstruction_status);
        assert(!candidate.message.empty());
    };

    // Moment data that does not fit the hierarchy is reported before the
    // reconstruction, which trusts its input.
    const auto assert_invalid = [&](const qmomkit::CqmomMomentData& candidate_moments, const std::string& message)
    {
        qmomkit::CqmomStateCandidate candidate;
        assert_cqmom_reconstructions(0, [&] { candidate = prepared_case.state_candidate(candidate_moments); });
        assert_failed(candidate, qmomkit::StateCandidateStatus::invalid_input, std::nullopt);
        assert(candidate.message == message);
    };
    qmomkit::CqmomMomentData short_data = moments;
    short_data.values.pop_back();
    assert_invalid(short_data, "CQMOM moment data value count does not match its required CQMOM moment layout.");
    qmomkit::CqmomMomentData nonfinite = moments;
    nonfinite.values[3] = std::numeric_limits<double>::quiet_NaN();
    assert_invalid(nonfinite, "CQMOM moment data values must be finite.");
    const qmomkit::PhaseSpaceHierarchy reordered({2, 2}, {1, 0});
    assert_invalid(grid_moment_data(reordered, {0.3, 0.2, 0.4, 0.1}, {{0.5, 0.5, 1.5, 1.5}, {-0.5, 2.0, 0.5, 1.0}}),
                   "CQMOM moment data layout must match the required CQMOM moment layout of the base hierarchy.");

    // Moment data without number density, as in an empty cell, is empty, also
    // when a check rejects boundary moments: its marginals are the zero measure.
    const qmomkit::CqmomMomentData zero_data{moments.layout, std::vector<double>(moments.values.size(), 0.0)};
    qmomkit::RealizabilityPolicy strict;
    strict.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    strict.accept_boundary = false;
    assert_failed(prepared_case.state_candidate(zero_data), qmomkit::StateCandidateStatus::empty, qmomkit::CqmomStatus::qmom_failed);
    assert_failed(prepared_case.state_candidate(zero_data, {strict, strict}), qmomkit::StateCandidateStatus::empty,
                  qmomkit::CqmomStatus::qmom_failed);

    // It is empty whatever the verdict on its marginals, also when they are
    // not realizable.
    const std::vector<qmomkit::MultiIndex> zero_indices = zero_data.layout.required_indices.to_multi_indices();
    qmomkit::CqmomMomentData unrealizable_zero_data = zero_data;
    unrealizable_zero_data
        .values[static_cast<std::size_t>(std::find(zero_indices.begin(), zero_indices.end(), qmomkit::MultiIndex{1, 0}) - zero_indices.begin())] =
        1.0;
    const qmomkit::CqmomStateCandidate emptied = prepared_case.state_candidate(unrealizable_zero_data, {strict, strict});
    assert_failed(emptied, qmomkit::StateCandidateStatus::empty, qmomkit::CqmomStatus::qmom_failed);
    assert(!emptied.realizability.has_value());

    // A negative variance of the first coordinate fails the fixed
    // reconstruction.
    const std::vector<qmomkit::MultiIndex> indices = moments.layout.required_indices.to_multi_indices();
    const auto position = [&indices](const qmomkit::MultiIndex& index)
    { return static_cast<std::size_t>(std::find(indices.begin(), indices.end(), index) - indices.begin()); };
    qmomkit::CqmomMomentData negative_variance = moments;
    negative_variance.values[position({2, 0})] = 0.1;
    const qmomkit::CqmomStateCandidate failed = prepared_case.state_candidate(negative_variance);
    assert_failed(failed, qmomkit::StateCandidateStatus::reconstruction_failed, qmomkit::CqmomStatus::qmom_failed);

    // The moments of the accepted quadrature can overflow although the given
    // ones do not: the first coordinate reduces to one node at 1e200, whose
    // second moment is 1e400.
    qmomkit::CqmomMomentData overflowing{moments.layout, std::vector<double>(moments.values.size(), 0.0)};
    overflowing.values[position({0, 0})] = 1.0;
    overflowing.values[position({1, 0})] = 1.0e200;
    const qmomkit::CqmomPreparedCase adaptive_case{qmomkit::CqmomCase{overflowing}};
    assert_failed(adaptive_case.state_candidate(overflowing), qmomkit::StateCandidateStatus::reconstruction_failed, qmomkit::CqmomStatus::success);

    // initial_state() is the throwing form of the same operation.
    bool threw = false;
    try
    {
        static_cast<void>(qmomkit::config::with_initial_moments(prepared_case, negative_variance).value().initial_state());
    }
    catch (const std::runtime_error& error)
    {
        threw = std::string(error.what()) == "CQMOM initial state quadrature reconstruction failed: " + failed.message;
    }
    assert(threw);
}

void test_cqmom_moment_layout_two_dimensional_reduced_set()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    assert(layout.hierarchy.number_of_dimensions() == 2);
    assert(layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(layout.hierarchy.conditioning_order() == std::vector<std::size_t>({0, 1}));
    assert_indices_equal(layout.required_indices,
                         std::vector<qmomkit::MultiIndex>({{0, 0}, {1, 0}, {2, 0}, {3, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}}));
    assert_unique_indices(layout.required_indices);
}

void test_cqmom_moment_layout_three_dimensional_recursive_indices()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    std::vector<qmomkit::MultiIndex> expected;
    for (std::size_t i = 0; i < 4; ++i)
    {
        append_unique(expected, {i, 0, 0});
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 4; ++j)
        {
            append_unique(expected, {i, j, 0});
        }
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            for (std::size_t k = 0; k < 4; ++k)
            {
                append_unique(expected, {i, j, k});
            }
        }
    }

    assert(expected.size() == 22);
    assert_indices_equal(layout.required_indices, expected);
    assert_unique_indices(layout.required_indices);
}

void test_cqmom_moment_layout_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    std::vector<qmomkit::MultiIndex> expected;
    for (std::size_t order_dim2 = 0; order_dim2 < 4; ++order_dim2)
    {
        append_unique(expected, {0, 0, order_dim2});
    }
    for (std::size_t order_dim2 = 0; order_dim2 < 2; ++order_dim2)
    {
        for (std::size_t order_dim0 = 0; order_dim0 < 4; ++order_dim0)
        {
            append_unique(expected, {order_dim0, 0, order_dim2});
        }
    }
    for (std::size_t order_dim2 = 0; order_dim2 < 2; ++order_dim2)
    {
        for (std::size_t order_dim0 = 0; order_dim0 < 2; ++order_dim0)
        {
            for (std::size_t order_dim1 = 0; order_dim1 < 4; ++order_dim1)
            {
                append_unique(expected, {order_dim0, order_dim1, order_dim2});
            }
        }
    }

    assert(layout.hierarchy.conditioning_order() == std::vector<std::size_t>({2, 0, 1}));
    assert(expected.size() == 22);
    assert_indices_equal(layout.required_indices, expected);
    assert_unique_indices(layout.required_indices);
}

void test_cqmom_moment_layout_unequal_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 4});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    assert(layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 3, 4}));
    assert(layout.required_indices.size() == 56);
    assert(contains_index(layout.required_indices, {3, 0, 0}));
    assert(!contains_index(layout.required_indices, {4, 0, 0}));
    assert(contains_index(layout.required_indices, {1, 5, 0}));
    assert(!contains_index(layout.required_indices, {2, 5, 0}));
    assert(contains_index(layout.required_indices, {1, 2, 7}));
    assert(!contains_index(layout.required_indices, {1, 3, 7}));
    assert_unique_indices(layout.required_indices);
}

void test_extract_cqmom_moments_takes_required_tensor_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor moment_tensor{extents, indexed_values(extents)};

    const qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments(moment_tensor, layout);

    assert(selected.layout.required_indices == layout.required_indices);
    assert(selected.values.size() == layout.required_indices.size());
    for (std::size_t offset = 0; offset < layout.required_indices.size(); ++offset)
    {
        const qmomkit::MultiIndexView index = layout.required_indices[offset];
        const double expected = moment_tensor.values[qmomkit::flatten_multi_index(index, moment_tensor.extents)];
        assert(selected.values[offset] == expected);
        assert(selected.contains_required_index(index));
        assert(selected.moment(index) == expected);
    }
    assert(!selected.contains_required_index({2, 2}));
    assert_throws_out_of_range([&selected]() { (void)selected.moment({2, 2}); });
}

void test_cqmom_moment_data_lookup_reflects_public_value_updates()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 4};
    qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments({extents, indexed_values(extents)}, layout);

    assert(selected.moment({1, 0}) == 14.0);
    selected.values[1] = 1234.0;
    assert(selected.moment({1, 0}) == 1234.0);
    assert(!selected.contains_required_index({0}));
}

void test_extract_cqmom_moments_respects_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 4, 4};
    const qmomkit::RawMomentTensor moment_tensor{extents, indexed_values(extents)};

    const qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments(moment_tensor, layout);

    assert(selected.values.size() == 22);
    for (std::size_t offset = 0; offset < layout.required_indices.size(); ++offset)
    {
        const qmomkit::MultiIndexView index = layout.required_indices[offset];
        assert(selected.values[offset] == moment_tensor.values[qmomkit::flatten_multi_index(index, moment_tensor.extents)]);
    }
    assert(selected.moment({0, 0, 3}) == moment_tensor.values[qmomkit::flatten_multi_index({0, 0, 3}, moment_tensor.extents)]);
    assert(selected.moment({3, 0, 1}) == moment_tensor.values[qmomkit::flatten_multi_index({3, 0, 1}, moment_tensor.extents)]);
    assert(selected.moment({1, 3, 1}) == moment_tensor.values[qmomkit::flatten_multi_index({1, 3, 1}, moment_tensor.extents)]);
}

void test_extract_cqmom_moments_handles_unequal_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 4});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{4, 6, 8};
    const qmomkit::RawMomentTensor moment_tensor{extents, indexed_values(extents)};

    const qmomkit::CqmomMomentData selected = qmomkit::extract_cqmom_moments(moment_tensor, layout);

    assert(selected.values.size() == 56);
    assert(selected.moment({3, 0, 0}) == moment_tensor.values[qmomkit::flatten_multi_index({3, 0, 0}, moment_tensor.extents)]);
    assert(selected.moment({1, 5, 0}) == moment_tensor.values[qmomkit::flatten_multi_index({1, 5, 0}, moment_tensor.extents)]);
    assert(selected.moment({1, 2, 7}) == moment_tensor.values[qmomkit::flatten_multi_index({1, 2, 7}, moment_tensor.extents)]);
}

void test_extract_cqmom_moments_reports_missing_tensor_moment()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<std::size_t> extents{3, 4};
    const qmomkit::RawMomentTensor moment_tensor{extents, indexed_values(extents)};

    assert_throws_out_of_range([&moment_tensor, &layout]() { (void)qmomkit::extract_cqmom_moments(moment_tensor, layout); });
}

void test_create_cqmom_moment_data_converts_valid_selected_result()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentRequest request{layout.required_indices.to_multi_indices()};
    const qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result(request, indexed_count_values(layout.required_indices.size()));

    const qmomkit::CqmomMomentData data = qmomkit::create_cqmom_moment_data(layout, selected);

    assert(data.layout.required_indices == layout.required_indices);
    assert(data.values == selected.values);
    assert(data.moment(layout.required_indices[0]) == selected.values.front());
}

void test_create_cqmom_moment_data_rejects_failed_selected_result()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices.to_multi_indices()}, indexed_count_values(layout.required_indices.size()));
    selected.status = qmomkit::RawMomentStatus::invalid_input;
    selected.message = "failed upstream";

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_index_order_mismatch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::RawMomentRequest request{layout.required_indices.to_multi_indices()};
    std::swap(request.indices.front(), request.indices.back());
    const qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result(request, indexed_count_values(layout.required_indices.size()));

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_missing_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices.to_multi_indices()}, indexed_count_values(layout.required_indices.size()));
    selected.values.pop_back();

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_extra_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices.to_multi_indices()}, indexed_count_values(layout.required_indices.size()));
    selected.values.push_back(123.0);

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_rejects_non_finite_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::SelectedRawMomentResult selected =
        qmomkit::make_selected_raw_moment_result({layout.required_indices.to_multi_indices()}, indexed_count_values(layout.required_indices.size()));
    selected.values.front() = std::numeric_limits<double>::infinity();

    assert_throws_invalid_argument([&layout, &selected]() { (void)qmomkit::create_cqmom_moment_data(layout, selected); });
}

void test_create_cqmom_moment_data_from_independent_product_selected_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::IndependentProductDistribution product{{qmomkit::GammaDistribution{1.5, 0.5}, qmomkit::BetaDistribution{2.0, 3.0}}, 2.0};

    const qmomkit::SelectedRawMomentResult selected =
        qmomkit::compute_raw_moments(product, qmomkit::RawMomentRequest{layout.required_indices.to_multi_indices()});
    const qmomkit::CqmomMomentData data = qmomkit::create_cqmom_moment_data(layout, selected);
    const qmomkit::TensorMomentResult tensor_result = qmomkit::compute_raw_moments(product, {3, 3});
    const qmomkit::CqmomMomentData required_moment_data = qmomkit::extract_cqmom_moments(tensor_result.moments, layout);

    assert(selected.ok());
    assert(tensor_result.ok());
    assert(data.layout.required_indices == layout.required_indices);
    assert(data.values.size() == required_moment_data.values.size());
    for (std::size_t index = 0; index < data.values.size(); ++index)
    {
        assert_near(data.values[index], required_moment_data.values[index], 1.0e-12);
    }
}

void test_cqmom_case_public_api_rejects_missing_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    qmomkit::CqmomMomentData moment_data = moment_data_from_tensor(hierarchy, {extents, indexed_values(extents)});
    moment_data.values.pop_back();

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomCase{moment_data});

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(result.message == "CQMOM moment data values must match the layout size.");
}

void test_cqmom_case_public_api_rejects_non_finite_values()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    qmomkit::CqmomMomentData moment_data = moment_data_from_tensor(hierarchy, {extents, indexed_values(extents)});
    moment_data.values.front() = std::numeric_limits<double>::infinity();

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomCase{moment_data});

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(result.message == "CQMOM moment data values must be finite.");
}

void test_cqmom_case_public_api_rejects_invalid_options()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::CqmomMomentData moment_data = moment_data_from_tensor(hierarchy, {extents, indexed_values(extents)});
    qmomkit::CqmomCase cqmom_case{moment_data};
    cqmom_case.conditional_options.singular_tolerance = 0.0;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(!result.message.empty());
}

void test_cqmom_case_public_api_rejects_incompatible_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const qmomkit::PhaseSpaceHierarchy incompatible_hierarchy({2, 2, 2}, {2, 0, 1});
    const qmomkit::CqmomMomentData moment_data = moment_data_from_tensor(incompatible_hierarchy, {{4, 4, 4}, indexed_values({4, 4, 4})});
    qmomkit::CqmomCase cqmom_case{moment_data};
    cqmom_case.moments.layout.hierarchy = hierarchy;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::invalid_input);
    assert(result.message == "CQMOM moment data layout must match the recursive hierarchy.");
}

template <typename Case> void assert_cqmom_rejects_unrepresentable_required_moment_extent(const Case& cqmom_case)
{
    bool leaked_length_error = false;
    try
    {
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);
        assert(!result.ok());
        assert(result.status == qmomkit::CqmomStatus::invalid_input);
        assert(result.message == "CQMOM node count is too large to derive the required moment extent.");
    }
    catch (const std::length_error&)
    {
        leaked_length_error = true;
    }
    assert(!leaked_length_error);
}

void test_cqmom_public_api_rejects_unrepresentable_required_moment_extent_before_storage()
{
    const std::size_t node_count = std::numeric_limits<std::size_t>::max() / 2 + 1;
    const qmomkit::PhaseSpaceHierarchy hierarchy({node_count});

    const qmomkit::CqmomMomentData moments{{hierarchy, {}}, {}};
    const qmomkit::CqmomCase cqmom_case{moments};
    assert_cqmom_rejects_unrepresentable_required_moment_extent(cqmom_case);
}

void test_qmom_public_api_matches_wheeler()
{
    // Public-wrapper consistency: the public result must preserve the underlying Wheeler result; analytical Wheeler correctness is verified
    // separately.
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.75}, {0.25, 0.75}, 4);

    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2);
    const qmomkit::Quadrature1D wheeler = qmomkit::compute_quadrature(qmomkit::MomentSet(moments), 2);

    assert(qmom.ok());
    assert(wheeler.ok());
    assert(qmom.abscissas.size() == 2);
    for (std::size_t node = 0; node < qmom.abscissas.size(); ++node)
    {
        assert_near(qmom.abscissas[node], wheeler.abscissas[node], 1.0e-13);
        assert_near(qmom.weights[node], wheeler.weights[node], 1.0e-13);
    }
}

// Counts the QMOM reconstructions of a call in instrumented builds; elsewhere
// the count is not checked.
template <typename Function> void assert_qmom_reconstructions(const std::size_t expected, Function function)
{
    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        function();
    }
    if (qmomkit::internal::performance::enabled())
    {
        assert(session.snapshot().reconstruction.qmom_reconstructions == expected);
    }
}

void test_qmom_state_candidate_matches_initial_state()
{
    // A state from other moments is the quadrature of those moments, with the
    // node count, options and bounds of the case, and the moments of that
    // quadrature, as many as were given. The moments of the case give the
    // state of initial_state().
    qmomkit::QmomOptions options;
    options.rmin = 1.0e-8;
    const qmomkit::AbscissaBounds bounds{-1.0, 4.0};
    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet(raw_moments_1d({-0.5, 1.75}, {0.25, 0.75}, 4)), 2, options, bounds};
    const qmomkit::MomentSet other(raw_moments_1d({0.25, 3.0}, {0.6, 0.4}, 6));
    qmomkit::QmomStateCandidate from_other;
    assert_qmom_reconstructions(1, [&] { from_other = qmom_case.state_candidate(other); });
    assert(from_other.ok());
    assert(!from_other.realizability.has_value());
    assert(from_other.reconstruction_status == qmomkit::QmomStatus::success);
    assert(from_other.message.empty());

    const qmomkit::QmomResult quadrature = qmomkit::compute_qmom(other, 2, options, bounds);
    const qmomkit::QmomState& state = *from_other.state;
    assert(state.quadrature().abscissas == quadrature.abscissas);
    assert(state.quadrature().weights == quadrature.weights);
    assert(state.moment_values() ==
           qmomkit::reconstruct_raw_moments(qmomkit::Quadrature1D{quadrature.abscissas, quadrature.weights, qmomkit::QuadratureStatus::success, ""},
                                            other.size()));
    assert(state.moment_values().size() == 6);
    assert(state.maximum_node_count() == 2);
    assert(state.qmom_options().rmin == options.rmin);
    assert(state.abscissa_bounds()->lower == bounds.lower);
    assert(state.abscissa_bounds()->upper == bounds.upper);

    const qmomkit::QmomStateCandidate own = qmom_case.state_candidate(qmom_case.moments);
    const qmomkit::QmomState initial = qmom_case.initial_state();
    assert(own.state->moment_values() == initial.moment_values());
    assert(own.state->quadrature().abscissas == initial.quadrature().abscissas);
    assert(own.state->quadrature().weights == initial.quadrature().weights);

    // The bounds of the case apply to the reconstruction: two points at 5 and
    // 6 lie outside them, as does any rule for their moments.
    assert(qmom_case.state_candidate(qmomkit::MomentSet(raw_moments_1d({5.0, 6.0}, {0.5, 0.5}, 4))).status ==
           qmomkit::StateCandidateStatus::reconstruction_failed);
}

void test_qmom_state_candidate_applies_the_realizability_policy()
{
    // A point at -1 is realizable on the real line, so the reconstruction
    // accepts it; only a requested check on a nonnegative support rejects it,
    // before any reconstruction.
    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet{1.0, 0.0}, 1};
    const qmomkit::MomentSet negative_point{2.0, -2.0};
    qmomkit::RealizabilityPolicy nonnegative;
    nonnegative.support = qmomkit::MomentSupport::nonnegative;
    assert(qmom_case.state_candidate(negative_point).ok());
    const qmomkit::QmomStateCandidate unchecked = qmom_case.state_candidate(negative_point, nonnegative);
    assert(unchecked.ok());
    assert(!unchecked.realizability.has_value());

    for (const qmomkit::RealizabilityCheckMode mode : {qmomkit::RealizabilityCheckMode::FinalOnly, qmomkit::RealizabilityCheckMode::EveryStage})
    {
        nonnegative.check_mode = mode;
        qmomkit::QmomStateCandidate rejected;
        assert_qmom_reconstructions(0, [&] { rejected = qmom_case.state_candidate(negative_point, nonnegative); });
        assert(rejected.status == qmomkit::StateCandidateStatus::realizability_rejected);
        assert(!rejected.state.has_value());
        assert(rejected.realizability->status == qmomkit::RealizabilityStatus::not_realizable);
        assert(!rejected.reconstruction_status.has_value());
        assert(rejected.message == rejected.realizability->message);
    }

    // An interval is read only for a bounded-interval support: a point at 2
    // is nonnegative whatever interval the policy also holds.
    nonnegative.interval = qmomkit::BoundedInterval{0.0, 1.0};
    const qmomkit::QmomStateCandidate accepted = qmom_case.state_candidate(qmomkit::MomentSet{1.0, 2.0}, nonnegative);
    assert(accepted.ok());
    assert(accepted.realizability->status == qmomkit::RealizabilityStatus::strictly_realizable);

    // A point at 0 is on the boundary of that support: accepted unless boundary
    // moments are excluded, and the verdict is kept either way.
    const qmomkit::MomentSet origin_point{2.0, 0.0};
    const qmomkit::QmomStateCandidate boundary = qmom_case.state_candidate(origin_point, nonnegative);
    assert(boundary.ok());
    assert(boundary.realizability->status == qmomkit::RealizabilityStatus::boundary);
    nonnegative.accept_boundary = false;
    const qmomkit::QmomStateCandidate strict_only = qmom_case.state_candidate(origin_point, nonnegative);
    assert(strict_only.status == qmomkit::StateCandidateStatus::realizability_rejected);
    assert(strict_only.realizability->status == qmomkit::RealizabilityStatus::boundary);

    // A bounded interval is checked against its bounds, and a policy without
    // them is invalid input.
    qmomkit::RealizabilityPolicy unit_interval;
    unit_interval.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    unit_interval.support = qmomkit::MomentSupport::bounded_interval;
    const qmomkit::QmomStateCandidate no_bounds = qmom_case.state_candidate(qmomkit::MomentSet{1.0, 2.0}, unit_interval);
    assert(no_bounds.status == qmomkit::StateCandidateStatus::invalid_input);
    assert(no_bounds.message == "QMOM realizability policy: A bounded-interval realizability policy requires an interval.");
    assert(!no_bounds.realizability.has_value());

    // An invalid policy is invalid input whatever the moments, also when they
    // are too few for the node count, and a check mode outside the enumeration
    // is invalid rather than applied.
    assert(qmom_case.state_candidate(qmomkit::MomentSet{1.0}, unit_interval).message == no_bounds.message);
    qmomkit::RealizabilityPolicy unknown_mode;
    unknown_mode.check_mode = static_cast<qmomkit::RealizabilityCheckMode>(7);
    const qmomkit::QmomStateCandidate unknown = qmom_case.state_candidate(qmomkit::MomentSet{1.0, 0.5}, unknown_mode);
    assert(unknown.status == qmomkit::StateCandidateStatus::invalid_input);
    assert(unknown.message == "QMOM realizability policy: Unknown realizability check mode.");
    unit_interval.interval = qmomkit::BoundedInterval{0.0, 1.0};
    assert(qmom_case.state_candidate(qmomkit::MomentSet{1.0, 2.0}, unit_interval).status == qmomkit::StateCandidateStatus::realizability_rejected);
    const qmomkit::QmomStateCandidate inside = qmom_case.state_candidate(qmomkit::MomentSet{1.0, 0.5}, unit_interval);
    assert(inside.ok());
    assert(inside.realizability->status == qmomkit::RealizabilityStatus::strictly_realizable);

    // The tolerance of the policy applies: with 0.4, the middle eigenvalue of
    // the standardized H_2 of two points at -1 and 1 is indeterminate, which
    // is not accepted.
    const qmomkit::QmomCase two_node_case{qmomkit::MomentSet{1.0, 0.0, 1.0, 0.0}, 2};
    qmomkit::RealizabilityPolicy coarse;
    coarse.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    const qmomkit::MomentSet two_points{1.0, 0.0, 1.0, 0.0, 1.0};
    assert(two_node_case.state_candidate(two_points, coarse).ok());
    coarse.options.tolerance = 0.4;
    const qmomkit::QmomStateCandidate indeterminate = two_node_case.state_candidate(two_points, coarse);
    assert(indeterminate.status == qmomkit::StateCandidateStatus::realizability_rejected);
    assert(indeterminate.realizability->status == qmomkit::RealizabilityStatus::numerically_indeterminate);
    coarse.accept_indeterminate = true;
    assert(two_node_case.state_candidate(two_points, coarse).ok());

    // An adaptive reconstruction keeps the nodes the moments support, so a
    // negative variance gives a single node at the mean without a check; the
    // check rejects it.
    const qmomkit::MomentSet negative_variance{1.0, 0.0, -1.0, 0.0};
    assert(two_node_case.state_candidate(negative_variance).state->active_node_count() == 1);
    qmomkit::RealizabilityPolicy real_line;
    real_line.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    assert(two_node_case.state_candidate(negative_variance, real_line).status == qmomkit::StateCandidateStatus::realizability_rejected);
}

void test_qmom_state_candidate_accepts_indeterminate_moments_when_allowed()
{
    // A log-normal with sigma = 0.05 has 12 moments that six-node Wheeler
    // inverts to rounding, but too few digits for the realizability check to
    // decide at order 11: the verdict is indeterminate.
    const qmomkit::RawMomentResult narrow = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{0.0, 0.05}, 11);
    assert(narrow.ok());
    const qmomkit::QmomCase six_node_case{narrow.moments, 6};

    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    policy.support = qmomkit::MomentSupport::nonnegative;
    const qmomkit::QmomStateCandidate rejected = six_node_case.state_candidate(narrow.moments, policy);
    assert(rejected.status == qmomkit::StateCandidateStatus::realizability_rejected);
    assert(rejected.realizability->status == qmomkit::RealizabilityStatus::numerically_indeterminate);

    policy.accept_indeterminate = true;
    const qmomkit::QmomStateCandidate accepted = six_node_case.state_candidate(narrow.moments, policy);
    assert(accepted.ok());
    assert(accepted.state->active_node_count() == 6);
    for (std::size_t k = 0; k < narrow.moments.size(); ++k)
    {
        assert(std::abs(accepted.state->moment_values()[k] / narrow.moments[k] - 1.0) < 1.0e-13);
    }
}

void test_qmom_state_candidate_reports_failures_without_throwing()
{
    qmomkit::QmomOptions fixed;
    fixed.adaptive = false;
    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet{1.0, 0.0, 1.0, 0.0}, 2, fixed};
    qmomkit::RealizabilityPolicy real_line;
    real_line.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    const auto assert_failed = [](const qmomkit::QmomStateCandidate& candidate, const qmomkit::StateCandidateStatus status,
                                  const std::optional<qmomkit::QmomStatus> reconstruction_status)
    {
        assert(candidate.status == status);
        assert(!candidate.ok());
        assert(!candidate.state.has_value());
        assert(candidate.reconstruction_status == reconstruction_status);
        assert(!candidate.message.empty());
    };

    // Too few moments are invalid input, checked or not; non-finite moments
    // are, even above order 2n - 1, which the reconstruction does not read.
    assert_failed(qmom_case.state_candidate(qmomkit::MomentSet{1.0, 0.0, -1.0}), qmomkit::StateCandidateStatus::invalid_input,
                  qmomkit::QmomStatus::invalid_input);
    assert_failed(qmom_case.state_candidate(qmomkit::MomentSet{1.0, 0.0, -1.0}, real_line), qmomkit::StateCandidateStatus::invalid_input,
                  qmomkit::QmomStatus::invalid_input);
    assert_failed(qmom_case.state_candidate(qmomkit::MomentSet{1.0, 0.0, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()}),
                  qmomkit::StateCandidateStatus::invalid_input, std::nullopt);

    // Moments without number density, as in an empty cell, are empty rather
    // than failed, with or without a check, also when the check rejects
    // boundary moments: such moments are boundary-realizable.
    qmomkit::RealizabilityPolicy strict = real_line;
    strict.accept_boundary = false;
    for (const qmomkit::MomentSet& empty_moments : {qmomkit::MomentSet{0.0, 0.0, 0.0, 0.0}, qmomkit::MomentSet{1.0e-31, 0.0, 0.0, 0.0}})
    {
        assert_failed(qmom_case.state_candidate(empty_moments), qmomkit::StateCandidateStatus::empty, qmomkit::QmomStatus::not_realizable);
        assert_failed(qmom_case.state_candidate(empty_moments, real_line), qmomkit::StateCandidateStatus::empty, qmomkit::QmomStatus::not_realizable);
        const qmomkit::QmomStateCandidate strictly_checked = qmom_case.state_candidate(empty_moments, strict);
        assert_failed(strictly_checked, qmomkit::StateCandidateStatus::empty, qmomkit::QmomStatus::not_realizable);
        assert(strictly_checked.realizability->status == qmomkit::RealizabilityStatus::boundary);
    }
    // Without number density, moments are empty whatever the verdict, also
    // when they are not realizable, as rounding-level moments of a cell
    // emptied by transport can be. A negative m0 is no empty cell, a
    // non-finite moment stays invalid input, and with number density a strict
    // check still rejects boundary moments.
    for (const qmomkit::MomentSet& not_realizable : {qmomkit::MomentSet{0.0, 1.0, 0.0, 0.0}, qmomkit::MomentSet{1.0e-31, 0.0, -1.0e-31, 0.0}})
    {
        const qmomkit::QmomStateCandidate emptied = qmom_case.state_candidate(not_realizable, real_line);
        assert_failed(emptied, qmomkit::StateCandidateStatus::empty, qmomkit::QmomStatus::not_realizable);
        assert(emptied.realizability->status == qmomkit::RealizabilityStatus::not_realizable);
    }
    assert_failed(qmom_case.state_candidate(qmomkit::MomentSet{-1.0e-31, 0.0, 1.0e-31, 0.0}, real_line),
                  qmomkit::StateCandidateStatus::realizability_rejected, std::nullopt);
    assert_failed(qmom_case.state_candidate(qmomkit::MomentSet{1.0e-31, std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}, real_line),
                  qmomkit::StateCandidateStatus::invalid_input, std::nullopt);
    assert_failed(qmom_case.state_candidate(qmomkit::MomentSet{1.0e-20, 0.0, 0.0, 0.0}, strict),
                  qmomkit::StateCandidateStatus::realizability_rejected, std::nullopt);

    const qmomkit::QmomStateCandidate negative_variance = qmom_case.state_candidate(qmomkit::MomentSet{1.0, 0.0, -1.0, 0.0});
    assert_failed(negative_variance, qmomkit::StateCandidateStatus::reconstruction_failed, qmomkit::QmomStatus::not_realizable);

    // The moments of the accepted quadrature can overflow although the given
    // ones do not: a single node at 1e200 has a second moment of 1e400.
    const qmomkit::QmomCase adaptive_case{qmomkit::MomentSet{1.0, 0.0, 1.0, 0.0}, 2};
    const qmomkit::QmomStateCandidate overflowing = adaptive_case.state_candidate(qmomkit::MomentSet{1.0, 1.0e200, 0.0, 0.0});
    assert_failed(overflowing, qmomkit::StateCandidateStatus::reconstruction_failed, qmomkit::QmomStatus::success);

    // initial_state() is the throwing form of the same operation.
    bool threw = false;
    try
    {
        static_cast<void>(qmomkit::QmomCase{qmomkit::MomentSet{1.0, 0.0, -1.0, 0.0}, 2, fixed}.initial_state());
    }
    catch (const std::runtime_error& error)
    {
        threw = std::string(error.what()) == "QMOM initial state reconstruction failed: " + negative_variance.message;
    }
    assert(threw);
}

void test_qmom_rejects_node_count_with_unrepresentable_moment_requirement()
{
    const std::size_t impossible_node_count = std::numeric_limits<std::size_t>::max() / 2 + 1;
    const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet{1.0, 0.0}, impossible_node_count);

    assert(!result.ok());
    assert(result.status == qmomkit::QmomStatus::invalid_input);
    assert(result.message.find("too large") != std::string::npos);
}

void test_qmom_rejects_zero_node_count()
{
    // A zero node count is invalid input for every combination of inversion
    // options. The adaptive bounded-support combination is the one that used to
    // report success with an empty quadrature, which a caller that checks ok()
    // before reading weights and abscissas cannot defend against.
    const qmomkit::MomentSet moments{1.0, 0.8, 2.8, 4.4};
    const qmomkit::AbscissaBounds bounds{-10.0, 10.0};

    const auto assert_rejected = [](const qmomkit::QmomResult& result)
    {
        assert(!result.ok());
        assert(result.status == qmomkit::QmomStatus::invalid_input);
        assert(result.message.find("positive") != std::string::npos);
        assert(result.weights.empty());
        assert(result.abscissas.empty());
    };

    for (const bool adaptive : {false, true})
    {
        qmomkit::QmomOptions options;
        options.adaptive = adaptive;

        assert_rejected(qmomkit::compute_qmom(moments, 0, options));
        assert_rejected(qmomkit::compute_qmom(moments, 0, options, bounds));
        assert_rejected(qmomkit::compute_qmom(qmomkit::QmomCase{moments, 0, options, bounds}));
    }
}

void test_qmom_rejects_fewer_than_2n_moments_with_or_without_bounds()
{
    // An n-node rule needs 2n moments, with or without bounds and in fixed or
    // adaptive mode: four moments support two nodes and nothing more. Bounds do
    // not turn a request the moments cannot support into a smaller rule, and do
    // not change how it is rejected. The largest counts would overflow 2n; they
    // are rejected without any attempt.
    const qmomkit::MomentSet moments{1.0, 0.8, 2.8, 4.4};
    const qmomkit::AbscissaBounds bounds{-10.0, 10.0};
    for (const bool adaptive : {false, true})
    {
        qmomkit::QmomOptions options;
        options.adaptive = adaptive;

        const qmomkit::QmomResult supported = qmomkit::compute_qmom(moments, 2, options, bounds);
        assert(supported.ok());
        assert(supported.weights.size() == 2);

        for (const std::size_t requested :
             {std::size_t{3}, std::size_t{64}, std::numeric_limits<std::size_t>::max() / 2 + 1, std::numeric_limits<std::size_t>::max()})
        {
            const qmomkit::QmomResult unbounded = qmomkit::compute_qmom(moments, requested, options);
            const qmomkit::QmomResult bounded = qmomkit::compute_qmom(moments, requested, options, bounds);
            for (const qmomkit::QmomResult* const result : {&unbounded, &bounded})
            {
                assert(!result->ok());
                assert(result->status == qmomkit::QmomStatus::invalid_input);
                assert(result->weights.empty());
            }
            const bool overflows = requested > std::numeric_limits<std::size_t>::max() / 2;
            assert(bounded.message == unbounded.message);
            assert(bounded.message ==
                   (overflows ? "max_nodes is too large to derive the required moment count." : "Moment sequence has fewer than 2*n entries."));
        }
    }
}

void test_qmom_recovers_bounded_rule_of_large_magnitude()
{
    const std::vector<double> abscissas{2.0e50, 5.0e50};
    const std::vector<double> weights{0.25, 0.75};
    const qmomkit::MomentSet moments(raw_moments_1d(abscissas, weights, 4));

    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, 2, {}, qmomkit::AbscissaBounds{1.0e50, 6.0e50});

    assert(result.ok());
    assert(result.abscissas.size() == abscissas.size());
    for (std::size_t node = 0; node < abscissas.size(); ++node)
    {
        assert_relative_near(result.abscissas[node], abscissas[node], 1.0e-12);
        assert_near(result.weights[node], weights[node], 1.0e-12);
    }
}

void test_qmom_is_unit_covariant()
{
    const std::vector<double> reference_abscissas{-0.5, 1.75};
    const std::vector<double> weights{0.25, 0.75};
    const qmomkit::QmomOptions options{false};
    const qmomkit::QmomResult reference = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d(reference_abscissas, weights, 4)), 2, options);
    assert(reference.ok());

    for (const double factor : std::vector<double>{1.0e-9, 1.0e-6, 1.0, 1.0e3, 1.0e6})
    {
        std::vector<double> scaled_abscissas = reference_abscissas;
        for (double& abscissa : scaled_abscissas)
        {
            abscissa *= factor;
        }
        const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d(scaled_abscissas, weights, 4)), 2, options);
        assert(result.ok());
        assert(result.weights.size() == reference.weights.size());
        for (std::size_t node = 0; node < result.weights.size(); ++node)
        {
            assert_near(result.weights[node], reference.weights[node], 1.0e-12);
            assert_near(result.abscissas[node] / factor, reference.abscissas[node], 1.0e-12);
        }
    }

    qmomkit::QmomOptions adaptive_options = options;
    adaptive_options.adaptive = true;
    const qmomkit::QmomResult adaptive_reference = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d({0.25}, {1.0}, 4)), 2, adaptive_options);
    assert(adaptive_reference.ok());
    assert(adaptive_reference.weights.size() == 1);
    for (const double factor : std::vector<double>{1.0e-9, 1.0e6})
    {
        const qmomkit::QmomResult adaptive_result =
            qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d({0.25 * factor}, {1.0}, 4)), 2, adaptive_options);
        assert(adaptive_result.ok());
        assert(adaptive_result.weights.size() == adaptive_reference.weights.size());
        assert_near(adaptive_result.abscissas.front() / factor, adaptive_reference.abscissas.front(), 1.0e-12);
    }
}

void test_fixed_cqmom_refuses_point_mass_branch_in_every_unit()
{
    // The second x-branch is a point mass in y, so its conditional moments
    // support one y-node, not the two a fixed count requires. The refusal must
    // not depend on the sign roundoff gives the degenerate recurrence entry,
    // which changes with the unit of y.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.25, 0.25};
    for (const double unit : {1.0e-6, 1.0e-3, 0.37, 1.0, 3.0, 1.0e3, 1.0e6})
    {
        const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5 * unit, 0.8 * unit, 0.3 * unit, 0.3 * unit}};
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(
            cqmom_case_from_tensor(hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, {}, qmomkit::QmomOptions{false}));
        assert(!result.ok());
        assert(result.status == qmomkit::CqmomStatus::qmom_failed);
        assert(result.failed_phase_space_dimension == 1);
    }
}

void test_one_dimensional_cqmom_agrees_with_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-1.0, 2.0}, {0.4, 0.6}, 4);
    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2);

    const qmomkit::PhaseSpaceHierarchy hierarchy({2});
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(qmomkit::CqmomCase{
        qmomkit::extract_cqmom_moments({{4}, moments}, qmomkit::create_cqmom_moment_layout(hierarchy)), {}, qmomkit::QmomOptions{false}});

    assert(qmom.ok());
    assert(cqmom.ok());
    assert(cqmom.quadrature.number_of_dimensions() == 1);
    assert(cqmom.quadrature.total_nodes() == 2);

    for (std::size_t node = 0; node < 2; ++node)
    {
        assert_near(cqmom.quadrature.weight({node}), qmom.weights[node], 1.0e-13);
        assert_near(cqmom.quadrature.abscissa(0, {node}), qmom.abscissas[node], 1.0e-13);
    }
}

void test_cqmom_public_api_two_dimensional_manufactured_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};

    const std::vector<double> x{-0.75, 1.25};
    const std::vector<double> y{0.2, 1.4};
    const std::vector<double> x_weights{0.35, 0.65};
    const std::vector<double> y_weights{0.8, 0.2};

    std::vector<double> weights(hierarchy.total_nodes(), 0.0);
    std::vector<std::vector<double>> abscissas(2, std::vector<double>(hierarchy.total_nodes(), 0.0));
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            const std::size_t flat = hierarchy.flatten_node_index({i, j});
            weights[flat] = x_weights[i] * y_weights[j];
            abscissas[0][flat] = x[i];
            abscissas[1][flat] = y[j];
        }
    }

    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);
    const auto cqmom_case =
        cqmom_case_from_tensor(hierarchy, {extents, moments}, {}, qmomkit::QmomOptions{false}, qmomkit::CqmomDiagnosticLevel::detailed);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    assert(result.ok());
    assert(result.status == qmomkit::CqmomStatus::success);
    assert(!result.diagnostics.empty());
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        assert(diagnostics.qmom_status == qmomkit::QmomStatus::success);
        assert(diagnostics.conditional_moment_status == qmomkit::ConditionalMomentStatus::success);
        if (diagnostics.conditioning_level + 1 < hierarchy.number_of_dimensions())
        {
            assert(diagnostics.conditional_moment_diagnostics.residual_norm.has_value());
            assert(diagnostics.conditional_moment_diagnostics.condition_estimate.has_value());
        }
    }

    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const std::size_t flat = hierarchy.flatten_node_index(node_index);
        assert_near(result.quadrature.weight(node_index), weights[flat], 1.0e-12);
        assert_near(result.quadrature.abscissa(0, node_index), abscissas[0][flat], 1.0e-12);
        assert_near(result.quadrature.abscissa(1, node_index), abscissas[1][flat], 1.0e-12);
    }

    const std::vector<double> first_node = result.quadrature.node({0, 0});
    const std::vector<double> last_node = result.quadrature.node({1, 1});
    assert(first_node.size() == 2);
    assert(last_node.size() == 2);
    assert_near(first_node[0], x[0], 1.0e-12);
    assert_near(first_node[1], y[0], 1.0e-12);
    assert_near(last_node[0], x[1], 1.0e-12);
    assert_near(last_node[1], y[1], 1.0e-12);
}

void test_cqmom_recovers_coordinates_of_disparate_magnitude()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> x{2.0e50, 5.0e50};
    const std::vector<double> y{-4.0e-50, 3.0e-50};
    const std::vector<double> x_weights{0.25, 0.75};
    const std::vector<double> y_weights{0.6, 0.4};
    std::vector<double> weights(hierarchy.total_nodes(), 0.0);
    std::vector<std::vector<double>> abscissas(2, std::vector<double>(hierarchy.total_nodes(), 0.0));
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            const std::size_t flat = hierarchy.flatten_node_index({i, j});
            weights[flat] = x_weights[i] * y_weights[j];
            abscissas[0][flat] = x[i];
            abscissas[1][flat] = y[j];
        }
    }

    const qmomkit::QmomOptions options{false};
    auto cqmom_case = cqmom_case_from_tensor(hierarchy, {extents, raw_moments_from_grid(hierarchy, extents, weights, abscissas)}, {}, options,
                                             qmomkit::CqmomDiagnosticLevel::detailed);
    cqmom_case.abscissa_bounds = {qmomkit::AbscissaBounds{1.0e50, 6.0e50}, qmomkit::AbscissaBounds{-5.0e-50, 4.0e-50}};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    const std::vector<double> reference_x{2.0, 5.0};
    const std::vector<double> reference_y{-4.0, 3.0};
    std::vector<std::vector<double>> reference_abscissas(2, std::vector<double>(hierarchy.total_nodes(), 0.0));
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            const std::size_t flat = hierarchy.flatten_node_index({i, j});
            reference_abscissas[0][flat] = reference_x[i];
            reference_abscissas[1][flat] = reference_y[j];
        }
    }
    const auto reference_cqmom_case =
        cqmom_case_from_tensor(hierarchy, {extents, raw_moments_from_grid(hierarchy, extents, weights, reference_abscissas)}, {}, options,
                               qmomkit::CqmomDiagnosticLevel::detailed);
    const qmomkit::CqmomResult reference = qmomkit::compute_cqmom(reference_cqmom_case);

    assert(result.ok());
    assert(reference.ok());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const std::size_t flat = hierarchy.flatten_node_index(node_index);
        assert_near(result.quadrature.weight(node_index), weights[flat], 1.0e-12);
        assert_relative_near(result.quadrature.abscissa(0, node_index), abscissas[0][flat], 1.0e-11);
        assert_relative_near(result.quadrature.abscissa(1, node_index), abscissas[1][flat], 1.0e-11);
        assert_near(result.quadrature.weight(node_index), reference.quadrature.weight(node_index), 1.0e-12);
        assert_relative_near(result.quadrature.abscissa(0, node_index) / 1.0e50, reference.quadrature.abscissa(0, node_index), 1.0e-11);
        assert_relative_near(result.quadrature.abscissa(1, node_index) / 1.0e-50, reference.quadrature.abscissa(1, node_index), 1.0e-11);
    }
    assert(!result.diagnostics.empty());
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (!diagnostics.bounded_support_interval_present)
        {
            continue;
        }
        const qmomkit::AbscissaBounds expected = cqmom_case.abscissa_bounds[diagnostics.phase_space_dimension].value();
        assert_relative_near(diagnostics.bounded_support_interval.lower, expected.lower, 1.0e-12);
        assert_relative_near(diagnostics.bounded_support_interval.upper, expected.upper, 1.0e-12);
    }
}

void test_cqmom_public_api_reports_conditional_moment_failure()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0}, {2.0, -1.0}};
    qmomkit::ConditionalMomentOptions conditional_options;
    conditional_options.singular_tolerance = 1.0;
    const auto cqmom_case = cqmom_case_from_tensor(hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, conditional_options,
                                                   qmomkit::QmomOptions{false});

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::conditional_moment_failed);
    assert(result.failed_phase_space_dimension == 0);
    // The public contract makes the quadrature diagnostic-only whenever the
    // overall CQMOM result fails; callers must gate it on result.ok().
    assert(!result.message.empty());
    assert(result.message.find("Vandermonde") != std::string::npos);
    assert(!result.diagnostics.empty());
    const auto failed_diagnostic =
        std::find_if(result.diagnostics.begin(), result.diagnostics.end(), [](const qmomkit::CqmomDiagnostics& diagnostics)
                     { return diagnostics.conditional_moment_status == qmomkit::ConditionalMomentStatus::singular_system; });
    assert(failed_diagnostic != result.diagnostics.end());
    assert(failed_diagnostic->conditioning_level == 0);
    assert(failed_diagnostic->phase_space_dimension == 0);
    assert(failed_diagnostic->parent_path.empty());
    assert(!failed_diagnostic->message.empty());
}

void test_cqmom_public_api_exposes_adaptive_branch_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    const std::vector<double> moments = raw_moments_from_terminals(extents, weights, abscissas);

    const auto cqmom_case =
        cqmom_case_from_tensor(hierarchy, {extents, moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, qmomkit::CqmomDiagnosticLevel::basic);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    assert(result.ok());
    assert(result.quadrature.adaptive_branch_counts_enabled());
    assert(!result.quadrature.has_rectangular_storage());
    assert(result.quadrature.total_nodes() == 3);
    assert(result.quadrature.terminal_paths().size() == 3);

    bool saw_reduced_branch = false;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        assert(diagnostics.adaptive_branch_counts_enabled);
        if (diagnostics.conditioning_level == 1 && diagnostics.parent_path == qmomkit::ConditioningPath({1}))
        {
            saw_reduced_branch = true;
            assert(diagnostics.requested_node_count == 2);
            assert(diagnostics.actual_node_count == 1);
            assert(diagnostics.adaptive_node_count_reduced);
        }
    }
    assert(saw_reduced_branch);
}

void test_cqmom_adaptive_root_branch_uses_relative_rmin_only()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const double density_scale = 0.05;
    const std::vector<double> weights(4, 0.25 * density_scale);
    const std::vector<std::vector<double>> abscissas{{10.0, 10.0, 20.0, 20.0}, {-1.0, 1.0, -1.0, 1.0}};
    const qmomkit::QmomOptions qmom_options{true, 0.1, 1.0e-10, 1.0e-30};
    const auto cqmom_case = cqmom_case_from_tensor(hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, {}, qmom_options,
                                                   qmomkit::CqmomDiagnosticLevel::basic);

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);

    assert(result.ok());
    assert(result.quadrature.total_nodes() == 4);
    assert_near(std::accumulate(result.quadrature.weights().begin(), result.quadrature.weights().end(), 0.0), density_scale, 1.0e-14);

    bool saw_unreduced_root = false;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level == 0)
        {
            saw_unreduced_root = true;
            assert(diagnostics.actual_node_count == 2);
            assert(!diagnostics.adaptive_node_count_reduced);
        }
    }
    assert(saw_unreduced_root);
}

void test_cqmom_diagnostic_levels_control_expensive_diagnostics()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.5, 1.5}, {0.25, 2.0, 0.25, 2.0}};
    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);

    auto default_cqmom_case = cqmom_case_from_tensor(hierarchy, {extents, moments}, {}, qmomkit::QmomOptions{false});
    assert(default_cqmom_case.diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    qmomkit::internal::performance::CounterSession default_session;
    const qmomkit::CqmomResult default_result = [&]()
    {
        qmomkit::internal::performance::CounterScope scope(default_session);
        return qmomkit::compute_cqmom(default_cqmom_case);
    }();
    assert(default_result.ok());
    assert(default_result.diagnostics.empty());
    assert(default_session.snapshot().linear_algebra.condition_estimates == 0);

    qmomkit::CqmomCase detailed_cqmom_case = default_cqmom_case;
    detailed_cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    const qmomkit::CqmomResult detailed_result = qmomkit::compute_cqmom(detailed_cqmom_case);
    assert(detailed_result.ok());
    assert_same_quadrature(default_result.quadrature, detailed_result.quadrature, 1.0e-13);
    assert(!detailed_result.diagnostics.empty());
    assert(detailed_result.diagnostics.front().conditional_moment_diagnostics.residual_norm.has_value());
    assert(detailed_result.diagnostics.front().conditional_moment_diagnostics.condition_estimate.has_value());

    qmomkit::CqmomCase basic_cqmom_case = default_cqmom_case;
    basic_cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    const qmomkit::CqmomResult basic_result = qmomkit::compute_cqmom(basic_cqmom_case);
    assert(basic_result.ok());
    assert_same_quadrature(default_result.quadrature, basic_result.quadrature, 1.0e-13);
    assert(basic_result.diagnostics.size() == detailed_result.diagnostics.size());
    for (const qmomkit::CqmomDiagnostics& diagnostics : basic_result.diagnostics)
    {
        assert(diagnostics.qmom_status == qmomkit::QmomStatus::success);
        assert(diagnostics.conditional_moment_status == qmomkit::ConditionalMomentStatus::success);
        assert(!diagnostics.conditional_moment_diagnostics.residual_norm.has_value());
        assert(!diagnostics.conditional_moment_diagnostics.condition_estimate.has_value());
        assert(!diagnostics.conditional_moment_diagnostics.conditioning_warning);
    }

    qmomkit::CqmomCase none_cqmom_case = default_cqmom_case;
    none_cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;
    const qmomkit::CqmomResult none_result = qmomkit::compute_cqmom(none_cqmom_case);
    assert(none_result.ok());
    assert_same_quadrature(default_result.quadrature, none_result.quadrature, 1.0e-13);
    assert(none_result.diagnostics.empty());
}

void test_cqmom_none_diagnostics_preserves_failure_context()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    auto silent_cqmom_case = cqmom_case_from_tensor(hierarchy, {{4, 4}, std::vector<double>(16, 0.0)}, {}, qmomkit::QmomOptions{false});
    silent_cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(silent_cqmom_case);

    assert(!result.ok());
    assert(result.status == qmomkit::CqmomStatus::qmom_failed);
    assert(!result.message.empty());
    assert(!result.diagnostics.empty());
    assert(result.diagnostics.front().qmom_status != qmomkit::QmomStatus::success);
    assert(!result.diagnostics.front().message.empty());
}

void test_cqmom_settings_format()
{
    const std::string formatted = qmomkit::format_cqmom_settings(true, qmomkit::VandermondeSolver::dense, qmomkit::CqmomDiagnosticLevel::detailed);
    assert(formatted.find("adaptive = enabled") != std::string::npos);
    assert(formatted.find("Vandermonde solver = dense") != std::string::npos);
    assert(formatted.find("diagnostic level = detailed") != std::string::npos);

    const std::string pivoted_qr =
        qmomkit::format_cqmom_settings(false, qmomkit::VandermondeSolver::pivoted_qr, qmomkit::CqmomDiagnosticLevel::basic);
    assert(pivoted_qr.find("adaptive = disabled") != std::string::npos);
    assert(pivoted_qr.find("Vandermonde solver = pivoted_qr") != std::string::npos);
    assert(pivoted_qr.find("diagnostic level = basic") != std::string::npos);

    const std::string svd = qmomkit::format_cqmom_settings(false, qmomkit::VandermondeSolver::svd, qmomkit::CqmomDiagnosticLevel::basic);
    assert(svd.find("Vandermonde solver = svd") != std::string::npos);

    const std::string moore_penrose =
        qmomkit::format_cqmom_settings(false, qmomkit::VandermondeSolver::moore_penrose, qmomkit::CqmomDiagnosticLevel::basic);
    assert(moore_penrose.find("Vandermonde solver = moore_penrose") != std::string::npos);
}

void test_cqmom_diagnostics_default_bounded_support_placeholder()
{
    const qmomkit::CqmomDiagnostics diagnostics;
    assert(!diagnostics.bounded_support_checked);
    assert(!diagnostics.bounded_support_interval_present);
    assert_near(diagnostics.bounded_support_interval.lower, 0.0, 0.0);
    assert_near(diagnostics.bounded_support_interval.upper, 0.0, 0.0);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"moment_set_public_container_behavior", test_moment_set_public_container_behavior},
            {"cqmom_prepared_case_and_state_assignment_preserve_independent_values",
             test_cqmom_prepared_case_and_state_assignment_preserve_independent_values},
            {"cqmom_moment_layout_two_dimensional_reduced_set", test_cqmom_moment_layout_two_dimensional_reduced_set},
            {"cqmom_moment_layout_three_dimensional_recursive_indices", test_cqmom_moment_layout_three_dimensional_recursive_indices},
            {"cqmom_moment_layout_non_identity_conditioning_order", test_cqmom_moment_layout_non_identity_conditioning_order},
            {"cqmom_moment_layout_unequal_node_counts", test_cqmom_moment_layout_unequal_node_counts},
            {"extract_cqmom_moments_takes_required_tensor_values", test_extract_cqmom_moments_takes_required_tensor_values},
            {"cqmom_moment_data_lookup_reflects_public_value_updates", test_cqmom_moment_data_lookup_reflects_public_value_updates},
            {"extract_cqmom_moments_respects_non_identity_conditioning_order", test_extract_cqmom_moments_respects_non_identity_conditioning_order},
            {"extract_cqmom_moments_handles_unequal_node_counts", test_extract_cqmom_moments_handles_unequal_node_counts},
            {"extract_cqmom_moments_reports_missing_tensor_moment", test_extract_cqmom_moments_reports_missing_tensor_moment},
            {"create_cqmom_moment_data_converts_valid_selected_result", test_create_cqmom_moment_data_converts_valid_selected_result},
            {"create_cqmom_moment_data_rejects_failed_selected_result", test_create_cqmom_moment_data_rejects_failed_selected_result},
            {"create_cqmom_moment_data_rejects_index_order_mismatch", test_create_cqmom_moment_data_rejects_index_order_mismatch},
            {"create_cqmom_moment_data_rejects_missing_values", test_create_cqmom_moment_data_rejects_missing_values},
            {"create_cqmom_moment_data_rejects_extra_values", test_create_cqmom_moment_data_rejects_extra_values},
            {"create_cqmom_moment_data_rejects_non_finite_values", test_create_cqmom_moment_data_rejects_non_finite_values},
            {"create_cqmom_moment_data_from_independent_product_selected_moments",
             test_create_cqmom_moment_data_from_independent_product_selected_moments},
            {"cqmom_case_promotes_to_state_and_exposes_quadrature", test_cqmom_case_promotes_to_state_and_exposes_quadrature},
            {"cqmom_state_candidate_matches_initial_state", test_cqmom_state_candidate_matches_initial_state},
            {"cqmom_state_candidate_applies_coordinate_policies", test_cqmom_state_candidate_applies_coordinate_policies},
            {"cqmom_state_candidate_reports_failures_without_throwing", test_cqmom_state_candidate_reports_failures_without_throwing},
            {"cqmom_case_public_api_rejects_missing_values", test_cqmom_case_public_api_rejects_missing_values},
            {"cqmom_case_public_api_rejects_non_finite_values", test_cqmom_case_public_api_rejects_non_finite_values},
            {"cqmom_case_public_api_rejects_invalid_options", test_cqmom_case_public_api_rejects_invalid_options},
            {"cqmom_case_public_api_rejects_incompatible_layout", test_cqmom_case_public_api_rejects_incompatible_layout},
            {"cqmom_public_api_rejects_unrepresentable_required_moment_extent_before_storage",
             test_cqmom_public_api_rejects_unrepresentable_required_moment_extent_before_storage},
            {"qmom_public_api_matches_wheeler", test_qmom_public_api_matches_wheeler},
            {"qmom_state_candidate_matches_initial_state", test_qmom_state_candidate_matches_initial_state},
            {"qmom_state_candidate_applies_the_realizability_policy", test_qmom_state_candidate_applies_the_realizability_policy},
            {"qmom_state_candidate_accepts_indeterminate_moments_when_allowed", test_qmom_state_candidate_accepts_indeterminate_moments_when_allowed},
            {"qmom_state_candidate_reports_failures_without_throwing", test_qmom_state_candidate_reports_failures_without_throwing},
            {"qmom_rejects_node_count_with_unrepresentable_moment_requirement", test_qmom_rejects_node_count_with_unrepresentable_moment_requirement},
            {"qmom_rejects_zero_node_count", test_qmom_rejects_zero_node_count},
            {"qmom_rejects_fewer_than_2n_moments_with_or_without_bounds", test_qmom_rejects_fewer_than_2n_moments_with_or_without_bounds},
            {"qmom_recovers_bounded_rule_of_large_magnitude", test_qmom_recovers_bounded_rule_of_large_magnitude},
            {"qmom_is_unit_covariant", test_qmom_is_unit_covariant},
            {"one_dimensional_cqmom_agrees_with_qmom", test_one_dimensional_cqmom_agrees_with_qmom},
            {"cqmom_public_api_two_dimensional_manufactured_case", test_cqmom_public_api_two_dimensional_manufactured_case},
            {"cqmom_recovers_coordinates_of_disparate_magnitude", test_cqmom_recovers_coordinates_of_disparate_magnitude},
            {"fixed_cqmom_refuses_point_mass_branch_in_every_unit", test_fixed_cqmom_refuses_point_mass_branch_in_every_unit},
            {"cqmom_public_api_reports_conditional_moment_failure", test_cqmom_public_api_reports_conditional_moment_failure},
            {"cqmom_public_api_exposes_adaptive_branch_counts", test_cqmom_public_api_exposes_adaptive_branch_counts},
            {"cqmom_adaptive_root_branch_uses_relative_rmin_only", test_cqmom_adaptive_root_branch_uses_relative_rmin_only},
            {"cqmom_diagnostic_levels_control_expensive_diagnostics", test_cqmom_diagnostic_levels_control_expensive_diagnostics},
            {"cqmom_none_diagnostics_preserves_failure_context", test_cqmom_none_diagnostics_preserves_failure_context},
            {"cqmom_settings_format", test_cqmom_settings_format},
            {"cqmom_diagnostics_default_bounded_support_placeholder", test_cqmom_diagnostics_default_bounded_support_placeholder},
        });
}
