// States created from moments the host changed itself, with and without a
// realizability check.
//
// What is verified
//   That `QmomCase::state_candidate(...)` and
//   `CqmomPreparedCase::state_candidate(...)` turn moments into the state of
//   the quadrature they describe, and that a QMOM state stores the moments of
//   that quadrature in place of given moments the reconstruction does not
//   read and keeps the maximum node count of its case, so that it can grow
//   again. That a realizability check runs before the reconstruction and
//   decides by its support, and that a CQMOM check names the rejected
//   coordinate.
//
// How the reference is constructed
//   By hand, from the moments of chosen weighted points, as a conservative
//   host update produces when it mixes populations. The cases are set up with
//   other moments, which a candidate does not use; the QMOM cases are
//   adaptive with at most two nodes, the CQMOM cases have two nodes per
//   coordinate.
//   - QMOM: 0.4 at 0.5 and 0.6 at 1.5, the mixture of two populations of one
//     point each. Its m4 and m5 are changed by +10 % and -10 %, as by a host
//     that transports them inexactly; the reconstruction does not read them.
//     The first four moments are realizable on the real line and on [0, 2],
//     but not on [0, 1], where m1 <= m0 is necessary and m1 = 1.1.
//   - One point at 0.5. A source with the constant rate M(target) - M(start)
//     then moves the weights linearly towards 0.5 at 0.5 and 0.5 at 1.5,
//     which SSP-RK3 integrates exactly, so after one step of 0.25 the state is
//     0.875 at 0.5 and 0.125 at 1.5.
//   - CQMOM: 0.15 at (0.5, 1), 0.2 at (0.5, 3), 0.3 at (1.5, 0.5) and 0.35 at
//     (1.5, 2), the equal mixture of two populations of two points each. The
//     marginal moments of xi_0 are realizable on [0, 2]; those of xi_1 are
//     not, since a measure on [0, 2] has 2 m1 - m2 >= 0 and here
//     2 m1 - m2 = -0.225.
//
// Error metric and pass criteria
//   Statuses, verdicts, node counts and the rejected coordinate must match
//   exactly; a rejection must leave the reconstruction status unset. Weights,
//   positions and stored moments must match the reference to 1.0e-12,
//   relative to max(|reference|, 1); observed errors are below 1e-14, the
//   rounding of the reconstruction.
//
// Reproduce
//   ctest -R verification_phase_space_evolution_state_candidate_
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

// The QMOM mixture: 0.4 at 0.5 and 0.6 at 1.5.
std::vector<double> mixture_weights()
{
    return {0.4, 0.6};
}

std::vector<double> mixture_nodes()
{
    return {0.5, 1.5};
}

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& nodes, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

qmomkit::QmomCase two_node_case()
{
    qmomkit::QmomCase qmom_case;
    qmom_case.moments = qmomkit::MomentSet({1.0, 0.0, 1.0, 0.0});
    qmom_case.node_count = 2;
    qmom_case.qmom_options.adaptive = true;
    return qmom_case;
}

qmomkit::RealizabilityPolicy real_line_policy()
{
    qmomkit::RealizabilityPolicy policy;
    policy.check_mode = qmomkit::RealizabilityCheckMode::FinalOnly;
    return policy;
}

qmomkit::RealizabilityPolicy interval_policy(const double upper)
{
    qmomkit::RealizabilityPolicy policy = real_line_policy();
    policy.support = qmomkit::MomentSupport::bounded_interval;
    policy.interval = qmomkit::BoundedInterval{0.0, upper};
    return policy;
}

void verify_quadrature(const qmomkit::QmomResult& quadrature, const std::vector<double>& weights, const std::vector<double>& nodes)
{
    QMOMKIT_VERIFY_CHECK(quadrature.abscissas.size() == nodes.size());
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas[node], nodes[node], 0.0, tolerance);
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights[node], weights[node], 0.0, tolerance);
    }
}

// The state of a candidate that must have succeeded.
template <typename Candidate> const auto& accepted_state(const Candidate& candidate)
{
    if (!candidate.ok() || !candidate.state.has_value())
    {
        qmomkit::verification::fail_check("candidate.ok() && candidate.state.has_value()", __FILE__, __LINE__, candidate.message);
    }
    return *candidate.state;
}

void verify_qmom_state_holds_the_quadrature_of_the_moments()
{
    const std::vector<double> exact = raw_moments(mixture_weights(), mixture_nodes(), 6);
    std::vector<double> transported = exact;
    transported[4] *= 1.1;
    transported[5] *= 0.9;

    const qmomkit::QmomStateCandidate candidate = two_node_case().state_candidate(qmomkit::MomentSet(transported));
    const qmomkit::QmomState& state = accepted_state(candidate);
    verify_quadrature(state.quadrature(), mixture_weights(), mixture_nodes());
    const std::vector<double>& stored = state.moment_values();
    QMOMKIT_VERIFY_CHECK(stored.size() == exact.size());
    for (std::size_t order = 0; order < exact.size(); ++order)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(stored[order], exact[order], 0.0, tolerance);
    }
}

void verify_qmom_policy_support_decides_the_verdict()
{
    const qmomkit::MomentSet moments(raw_moments(mixture_weights(), mixture_nodes(), 4));
    const qmomkit::QmomCase qmom_case = two_node_case();
    for (const qmomkit::RealizabilityPolicy& policy : {real_line_policy(), interval_policy(2.0)})
    {
        const qmomkit::QmomStateCandidate candidate = qmom_case.state_candidate(moments, policy);
        verify_quadrature(accepted_state(candidate).quadrature(), mixture_weights(), mixture_nodes());
        QMOMKIT_VERIFY_CHECK(candidate.realizability.has_value() &&
                             candidate.realizability->status == qmomkit::RealizabilityStatus::strictly_realizable);
    }

    const qmomkit::QmomStateCandidate rejected = qmom_case.state_candidate(moments, interval_policy(1.0));
    QMOMKIT_VERIFY_CHECK(rejected.status == qmomkit::StateCandidateStatus::realizability_rejected);
    QMOMKIT_VERIFY_CHECK(!rejected.state.has_value());
    QMOMKIT_VERIFY_CHECK(rejected.realizability.has_value() && rejected.realizability->status == qmomkit::RealizabilityStatus::not_realizable);
    // The check runs first, so a rejection leaves the reconstruction unrun.
    QMOMKIT_VERIFY_CHECK(!rejected.reconstruction_status.has_value());
}

void verify_qmom_state_keeps_the_maximum_node_count()
{
    const std::vector<double> start = raw_moments({1.0}, {0.5}, 4);
    const std::vector<double> target = raw_moments({0.5, 0.5}, mixture_nodes(), 4);
    std::vector<double> rate(start.size());
    for (std::size_t order = 0; order < rate.size(); ++order)
    {
        rate[order] = target[order] - start[order];
    }

    const qmomkit::QmomStateCandidate candidate = two_node_case().state_candidate(qmomkit::MomentSet(start));
    const qmomkit::QmomState& state = accepted_state(candidate);
    QMOMKIT_VERIFY_CHECK(state.active_node_count() == 1);

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::QmomMomentSource constant_rate = [&rate](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{rate}; };
    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 0.25, zero_velocity, constant_rate);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.state.active_node_count() == 2);
    verify_quadrature(result.state.quadrature(), {0.875, 0.125}, mixture_nodes());
}

struct WeightedPoint
{
    double weight = 0.0;
    std::array<double, 2> position{};
};

// Equal mixture of {0.3 at (0.5, 1), 0.7 at (1.5, 2)} and
// {0.4 at (0.5, 3), 0.6 at (1.5, 0.5)}, in the order of the reconstruction:
// by xi_0, then by xi_1 within each xi_0 node.
std::vector<WeightedPoint> cqmom_mixture()
{
    return {{0.15, {0.5, 1.0}}, {0.2, {0.5, 3.0}}, {0.3, {1.5, 0.5}}, {0.35, {1.5, 2.0}}};
}

qmomkit::CqmomMomentData cqmom_moments(const std::vector<WeightedPoint>& points)
{
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    std::vector<double> values;
    for (const qmomkit::MultiIndexView index : layout.required_indices)
    {
        double value = 0.0;
        for (const WeightedPoint& point : points)
        {
            value += point.weight * std::pow(point.position[0], static_cast<int>(index[0])) * std::pow(point.position[1], static_cast<int>(index[1]));
        }
        values.push_back(value);
    }
    return qmomkit::CqmomMomentData{std::move(layout), std::move(values)};
}

// The case is set up with other moments, which a candidate does not use.
qmomkit::CqmomPreparedCase cqmom_case()
{
    return qmomkit::CqmomPreparedCase(qmomkit::CqmomCase{cqmom_moments({{0.3, {0.5, 1.0}}, {0.7, {1.5, 2.0}}})});
}

void verify_cqmom_state_holds_the_quadrature_of_the_moments()
{
    const qmomkit::CqmomMomentData moments = cqmom_moments(cqmom_mixture());
    const qmomkit::CqmomStateCandidate candidate = cqmom_case().state_candidate(moments);
    const qmomkit::CqmomState& state = accepted_state(candidate);

    const std::vector<WeightedPoint> points = cqmom_mixture();
    const qmomkit::CqmomQuadrature& quadrature = state.quadrature();
    QMOMKIT_VERIFY_CHECK(quadrature.total_nodes() == points.size());
    for (std::size_t node = 0; node < points.size(); ++node)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.weights()[node], points[node].weight, 0.0, tolerance);
        for (std::size_t coordinate = 0; coordinate < 2; ++coordinate)
        {
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(quadrature.abscissas()[coordinate][node], points[node].position.at(coordinate), 0.0, tolerance);
        }
    }
    const std::vector<double>& stored = state.moment_values();
    QMOMKIT_VERIFY_CHECK(stored.size() == moments.values.size());
    for (std::size_t entry = 0; entry < stored.size(); ++entry)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(stored[entry], moments.values[entry], 0.0, tolerance);
    }
}

void verify_cqmom_policy_names_the_rejected_coordinate()
{
    const qmomkit::CqmomMomentData moments = cqmom_moments(cqmom_mixture());
    const qmomkit::CqmomPreparedCase prepared_case = cqmom_case();

    const qmomkit::CqmomStateCandidate rejected = prepared_case.state_candidate(moments, {interval_policy(2.0), interval_policy(2.0)});
    QMOMKIT_VERIFY_CHECK(rejected.status == qmomkit::StateCandidateStatus::realizability_rejected);
    QMOMKIT_VERIFY_CHECK(!rejected.state.has_value());
    QMOMKIT_VERIFY_CHECK(rejected.rejected_coordinate == std::optional<std::size_t>{1});
    QMOMKIT_VERIFY_CHECK(rejected.realizability.has_value() && rejected.realizability->status == qmomkit::RealizabilityStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(!rejected.reconstruction_status.has_value());

    const qmomkit::CqmomStateCandidate accepted = prepared_case.state_candidate(moments, {interval_policy(2.0), std::nullopt});
    QMOMKIT_VERIFY_CHECK(accepted.ok() && accepted.state.has_value());
    QMOMKIT_VERIFY_CHECK(!accepted.rejected_coordinate.has_value());
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"qmom_state_holds_the_quadrature_of_the_moments", verify_qmom_state_holds_the_quadrature_of_the_moments},
            {"qmom_policy_support_decides_the_verdict", verify_qmom_policy_support_decides_the_verdict},
            {"qmom_state_keeps_the_maximum_node_count", verify_qmom_state_keeps_the_maximum_node_count},
            {"cqmom_state_holds_the_quadrature_of_the_moments", verify_cqmom_state_holds_the_quadrature_of_the_moments},
            {"cqmom_policy_names_the_rejected_coordinate", verify_cqmom_policy_names_the_rejected_coordinate},
        });
}
