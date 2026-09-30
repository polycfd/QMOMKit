/*
 * WHAT IS TESTED
 *   That CQMOM reconstruction and phase-space evolution produce exactly the
 *   results recorded here, across a matrix spanning phase-space dimension,
 *   conditioning order, node counts, adaptive branch counts, bounded coordinate
 *   support, every Vandermonde solver policy, and all three splitting schemes.
 *
 * HOW IT IS TESTED
 *   Two independent checks per static case.
 *
 *   1. Analytical recovery. Every case is built from an explicit list of
 *      quadrature nodes, and its raw moments are computed from those nodes.
 *      That node list is the only quadrature reproducing those moments for the
 *      configured node counts, so reconstruction must return the nodes the case
 *      was built from. This is an analytical check: it does not depend on any
 *      recorded value, and would catch a wrong result even if the recorded
 *      values were themselves wrong.
 *
 *   2. Prepared-path agreement. The same case is reconstructed a second time
 *      through CqmomCase, which reuses a reconstruction plan built once from
 *      the hierarchy instead of deriving it per call. Both routes must return
 *      the identical quadrature, down to the last bit.
 *
 *   3. Recorded values. The accepted weights and abscissas are compared against
 *      the values recorded below, which pins the result far more tightly than
 *      recovery alone: it also fixes node ordering, conditioning-path ordering,
 *      which branches adaptive reconstruction kept, and the last bit of every
 *      value.
 *
 *   Evolution cases carry only the recorded check, because once the state has
 *   been advanced it is no longer the node list the case was built from.
 *
 * PASS CRITERIA
 *   Recovery: every node the case was built from is matched by an accepted
 *   node, with weight and coordinates agreeing to 1e-12 absolute.
 *   Recorded values: status, node count and every weight and abscissa compare
 *   bit-exactly against the recorded value.
 *
 *   The recorded comparison is deliberately exact rather than using a
 *   tolerance. Its purpose is to prove that a refactor changed nothing, and a
 *   tolerance would hide exactly the drift it exists to detect. A failure
 *   therefore reports how many representable floating-point values separate the
 *   produced and recorded numbers, so that a real numerical change (many steps)
 *   is immediately distinguishable from a compiler or platform difference (one
 *   or two steps).
 *
 * WHEN IT IS TESTED
 *   Debug and Release. The checks are explicit and do not rely on assert().
 */

#include "../verification/verification_helpers.hpp"
#include "cqmom_case_matrix.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/phase_space_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{

constexpr double node_recovery_tolerance = 1.0e-12;

struct ReferenceQuadrature
{
    const char* name;
    int status;
    std::vector<double> weights;
    // Abscissas in dimension-major order: all nodes of dimension 0, then all
    // nodes of dimension 1, and so on.
    std::vector<double> abscissas;
};

#include "cqmom_reference_values.inc"

[[nodiscard]] const ReferenceQuadrature& reference_for(const std::string& name)
{
    for (const ReferenceQuadrature& reference : reference_quadratures())
    {
        if (name == reference.name)
        {
            return reference;
        }
    }
    std::cerr << "No recorded reference for case '" << name << "'.\n";
    std::exit(EXIT_FAILURE);
}

// How many representable floating-point values lie between the two numbers, so
// that a failure says whether the result moved by a rounding step or by
// something that matters.
[[nodiscard]] long double representable_values_between(const double actual, const double expected)
{
    if (actual == expected)
    {
        return 0.0L;
    }
    if (!std::isfinite(actual) || !std::isfinite(expected))
    {
        return std::numeric_limits<long double>::infinity();
    }
    std::int64_t actual_bits = 0;
    std::int64_t expected_bits = 0;
    std::memcpy(&actual_bits, &actual, sizeof(actual_bits));
    std::memcpy(&expected_bits, &expected, sizeof(expected_bits));
    if ((actual_bits < 0) != (expected_bits < 0))
    {
        return std::numeric_limits<long double>::infinity();
    }
    return static_cast<long double>(std::llabs(actual_bits - expected_bits));
}

void check_recorded_value(const std::string& case_name, const char* quantity, const std::size_t position, const double actual, const double expected)
{
    if (actual == expected)
    {
        return;
    }
    std::ostringstream detail;
    detail << std::setprecision(17) << "case " << case_name << ", " << quantity << '[' << position << "]: produced " << actual << ", recorded "
           << expected << ", " << static_cast<double>(representable_values_between(actual, expected)) << " representable values apart";
    qmomkit::verification::fail_check("recorded value comparison", __FILE__, __LINE__, detail.str());
}

void compare_against_reference(const std::string& case_name, const qmomkit::CqmomQuadrature& quadrature, const int status)
{
    const ReferenceQuadrature& reference = reference_for(case_name);
    QMOMKIT_VERIFY_CHECK(status == reference.status);
    QMOMKIT_VERIFY_CHECK(quadrature.weights().size() == reference.weights.size());

    for (std::size_t node = 0; node < reference.weights.size(); ++node)
    {
        check_recorded_value(case_name, "weight", node, quadrature.weights()[node], reference.weights[node]);
    }

    const std::size_t dimensions = quadrature.number_of_dimensions();
    QMOMKIT_VERIFY_CHECK(dimensions * reference.weights.size() == reference.abscissas.size());
    std::size_t position = 0;
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            check_recorded_value(case_name, "abscissa", position, quadrature.abscissas()[dimension][node], reference.abscissas[position]);
            ++position;
        }
    }
}

/*
 * The nodes a case was built from are the only quadrature reproducing its
 * moments, so each of them must appear among the accepted nodes. They are
 * matched by nearest position rather than by index, because node ordering is
 * already fixed by the recorded comparison above.
 *
 * Each accepted node may be claimed only once. Without that, a reconstruction
 * that returned the same node twice and dropped another could still match every
 * node the case was built from, provided the dropped one happened to lie close
 * to the duplicate.
 */
void check_recovers_generating_nodes(const qmomkit::case_matrix::StaticCase& case_, const qmomkit::CqmomQuadrature& quadrature)
{
    const std::size_t dimensions = quadrature.number_of_dimensions();
    QMOMKIT_VERIFY_CHECK(quadrature.weights().size() == case_.nodes.size());
    std::vector<bool> claimed(quadrature.weights().size(), false);

    for (const qmomkit::case_matrix::QuadratureNode& expected : case_.nodes)
    {
        std::size_t closest = 0;
        double closest_distance = std::numeric_limits<double>::infinity();
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            if (claimed[node])
            {
                continue;
            }
            double distance = 0.0;
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            {
                const double offset = quadrature.abscissas()[dimension][node] - expected.coordinate[dimension];
                distance += offset * offset;
            }
            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest = node;
            }
        }

        QMOMKIT_VERIFY_CHECK(closest_distance < std::numeric_limits<double>::infinity());
        claimed[closest] = true;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            QMOMKIT_VERIFY_NEAR(quadrature.abscissas()[dimension][closest], expected.coordinate[dimension], node_recovery_tolerance);
        }
        QMOMKIT_VERIFY_NEAR(quadrature.weights()[closest], expected.weight, node_recovery_tolerance);
    }
}

/*
 * Reconstructing through CqmomCase reuses a reconstruction plan built once from
 * the hierarchy, where the explicit case derives one per call. The two routes
 * must agree exactly: the plan is bookkeeping about which moments are needed
 * and in what order, and does not enter the arithmetic.
 */
void check_prepared_path_agrees(const qmomkit::CqmomResult& direct, const qmomkit::CqmomResult& prepared)
{
    QMOMKIT_VERIFY_CHECK(direct.status == prepared.status);
    QMOMKIT_VERIFY_CHECK(direct.quadrature.weights().size() == prepared.quadrature.weights().size());
    QMOMKIT_VERIFY_CHECK(direct.quadrature.number_of_dimensions() == prepared.quadrature.number_of_dimensions());

    for (std::size_t node = 0; node < direct.quadrature.weights().size(); ++node)
    {
        QMOMKIT_VERIFY_CHECK(direct.quadrature.weights()[node] == prepared.quadrature.weights()[node]);
        for (std::size_t dimension = 0; dimension < direct.quadrature.number_of_dimensions(); ++dimension)
        {
            QMOMKIT_VERIFY_CHECK(direct.quadrature.abscissas()[dimension][node] == prepared.quadrature.abscissas()[dimension][node]);
        }
    }
}

void run_static_case(const qmomkit::case_matrix::StaticCase& case_)
{
    const qmomkit::CqmomLeanCase explicit_case = qmomkit::case_matrix::lean_case_of(case_);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(explicit_case);
    QMOMKIT_VERIFY_CHECK(result.ok());

    check_recovers_generating_nodes(case_, result.quadrature);
    check_prepared_path_agrees(result, qmomkit::compute_cqmom(qmomkit::CqmomCase(explicit_case)));
    compare_against_reference(case_.name, result.quadrature, static_cast<int>(result.status));
}

void run_evolution_case(const qmomkit::case_matrix::EvolutionCase& case_)
{
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::case_matrix::matrix_velocity();
    const qmomkit::CqmomMomentSource source = qmomkit::case_matrix::matrix_source();

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = case_.scheme;

    qmomkit::CqmomState state = qmomkit::CqmomCase(qmomkit::case_matrix::lean_case_of(case_.reconstruction)).initial_state();
    double time = 0.0;
    int status = 0;
    for (std::size_t step = 0; step < case_.step_count; ++step)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult evolved =
            qmomkit::advance_cqmom_phase_space(state, time, case_.time_step, velocity, source, options);
        status = static_cast<int>(evolved.status);
        QMOMKIT_VERIFY_CHECK(evolved.ok());
        state = evolved.state;
        time = evolved.final_time;
    }

    compare_against_reference(case_.reconstruction.name, state.quadrature(), status);
}

} // namespace

int main(const int argc, char** argv)
{
    const std::string selected_case = argc == 2 ? argv[1] : "";
    const auto selected = [&selected_case](const std::string& name) { return selected_case.empty() || selected_case == name; };
    bool matched = selected_case.empty();

    for (const qmomkit::case_matrix::StaticCase& case_ : qmomkit::case_matrix::static_cases())
    {
        if (selected(case_.name))
        {
            matched = true;
            run_static_case(case_);
        }
    }

    for (const qmomkit::case_matrix::EvolutionCase& case_ : qmomkit::case_matrix::evolution_cases())
    {
        if (selected(case_.reconstruction.name))
        {
            matched = true;
            run_evolution_case(case_);
        }
    }

    if (!matched)
    {
        std::cerr << "Unknown CQMOM reference-value case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
