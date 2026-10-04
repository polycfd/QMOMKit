#pragma once

/*
 * CQMOM case matrix for the reference-value regression test.
 *
 * The matrix is defined here, separately from the test that runs it, so that
 * the cases and the values recorded for them stay independent of each other.
 *
 * Every case is built from an explicit list of quadrature nodes, and the raw
 * moments are computed from those nodes. Each case is therefore realizable by
 * construction, and its moments are exact to round-off rather than sampled.
 */

#include "qmomkit/abscissa_bounds.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/phase_space_evolution.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit::case_matrix
{

// One quadrature node: a weight and the phase-space position it sits at.
struct QuadratureNode
{
    double weight = 0.0;
    std::vector<double> coordinate;
};

// A static reconstruction case: moments in, quadrature out.
struct StaticCase
{
    std::string name;
    std::vector<std::size_t> node_counts;
    std::vector<std::size_t> conditioning_order;
    std::vector<QuadratureNode> nodes;
    bool adaptive_branch_counts = false;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
    // Empty leaves every coordinate unconstrained.
    CqmomAbscissaBounds coordinate_support;
};

// A phase-space evolution case: the same reconstruction driven for a few steps
// under one splitting scheme.
struct EvolutionCase
{
    StaticCase reconstruction;
    PhaseSpaceEvolutionScheme scheme = PhaseSpaceEvolutionScheme::strang;
    double time_step = 0.0;
    std::size_t step_count = 0;
};

[[nodiscard]] std::vector<StaticCase> static_cases();
[[nodiscard]] std::vector<EvolutionCase> evolution_cases();

// The evolution callbacks are deliberately simple, deterministic and
// dimension-generic: the point of these cases is the reconstruction and the
// splitting, not the physics.
[[nodiscard]] PhaseSpaceVelocity matrix_velocity();
[[nodiscard]] CqmomMomentSource matrix_source();

[[nodiscard]] CqmomCase cqmom_case_of(const StaticCase& case_);

} // namespace qmomkit::case_matrix
