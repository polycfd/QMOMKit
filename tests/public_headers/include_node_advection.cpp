#include <cstddef>
#include <functional>
#include <qmomkit/node_advection.hpp>
#include <type_traits>
#include <vector>

// The phase-space velocity contract is the library's hottest public callback,
// so both its packed shape and the per-node convenience shape are pinned here.
static_assert(
    std::is_same_v<qmomkit::PhaseSpaceVelocity, std::function<void(double, const std::vector<double>&, std::vector<double>&, std::size_t)>>);
static_assert(std::is_same_v<qmomkit::PerNodePhaseSpaceVelocity, std::function<std::vector<double>(double, const std::vector<double>&)>>);

using MakePhaseSpaceVelocity = qmomkit::PhaseSpaceVelocity (*)(qmomkit::PerNodePhaseSpaceVelocity);
static_assert(std::is_same_v<decltype(&qmomkit::make_phase_space_velocity), MakePhaseSpaceVelocity>);

// Node advection accepts the packed velocity through every quadrature form.
using QmomQuadratureAdvance = qmomkit::QmomNodeAdvectionResult (*)(const qmomkit::Quadrature1D&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                                   const qmomkit::NodeAdvectionOptions&);
using QmomResultAdvance = qmomkit::QmomNodeAdvectionResult (*)(const qmomkit::QmomResult&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                               const qmomkit::NodeAdvectionOptions&);
using CqmomAdvance = qmomkit::CqmomNodeAdvectionResult (*)(const qmomkit::CqmomQuadrature&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                           const qmomkit::NodeAdvectionOptions&);

[[maybe_unused]] QmomQuadratureAdvance qmom_quadrature_advance = qmomkit::advance_nodes;
[[maybe_unused]] QmomResultAdvance qmom_result_advance = qmomkit::advance_nodes;
[[maybe_unused]] CqmomAdvance cqmom_advance = qmomkit::advance_nodes;
