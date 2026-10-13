#pragma once

#include "qmomkit/phase_space_evolution.hpp"

#include <iosfwd>

namespace cqmom_2d_adaptive_node_counts
{

/**
 * Prints one node-count change that QMOMKit reported for a branch: the time,
 * the parent path of the branch, and its active node counts before and after.
 */
void print_transition(std::ostream& stream, double time, const qmomkit::CqmomNodeCountTransition& transition);

} // namespace cqmom_2d_adaptive_node_counts
