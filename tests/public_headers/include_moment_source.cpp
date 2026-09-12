#include <qmomkit/moment_source.hpp>
#include <type_traits>

static_assert(std::is_default_constructible_v<qmomkit::MomentBoundaryRate>);
static_assert(std::is_default_constructible_v<qmomkit::MomentBoundaryContribution>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomMomentBoundaryContribution>);

using CqmomTerminalNodeVisitor = void (qmomkit::CqmomSourceContext::*)(const std::function<void(double, const std::vector<double>&)>&) const;
using CqmomGeneralizedMoment = double (qmomkit::CqmomSourceContext::*)(const std::vector<double>&) const;
using QmomGeneralizedMoment = double (qmomkit::QmomSourceContext::*)(double) const;

static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::for_each_terminal_node), CqmomTerminalNodeVisitor>);
static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::generalized_moment), CqmomGeneralizedMoment>);
static_assert(std::is_same_v<decltype(&qmomkit::QmomSourceContext::generalized_moment), QmomGeneralizedMoment>);
