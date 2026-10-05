#include <qmomkit/moment_source.hpp>
#include <type_traits>

static_assert(std::is_default_constructible_v<qmomkit::QmomMomentRate>);
static_assert(std::is_default_constructible_v<qmomkit::QmomMomentBoundaryContribution>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomMomentBoundaryContribution>);

using CqmomTerminalNodeVisitor = void (qmomkit::CqmomSourceContext::*)(const std::function<void(double, const std::vector<double>&)>&) const;
using CqmomGeneralizedMoment = double (qmomkit::CqmomSourceContext::*)(const std::vector<double>&) const;
using QmomGeneralizedMoment = double (qmomkit::QmomSourceContext::*)(double) const;

static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::for_each_terminal_node), CqmomTerminalNodeVisitor>);
static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::generalized_moment), CqmomGeneralizedMoment>);
static_assert(std::is_same_v<decltype(&qmomkit::QmomSourceContext::generalized_moment), QmomGeneralizedMoment>);

// The CQMOM source context hands out rows of a flat index table rather than
// owning indices, so the accessor shapes are pinned here.
using CqmomIndices = const qmomkit::MultiIndexTable& (qmomkit::CqmomSourceContext::*)() const noexcept;
using CqmomIndexAt = qmomkit::MultiIndexView (qmomkit::CqmomSourceContext::*)(std::size_t) const;
using CqmomMomentAt = double (qmomkit::CqmomSourceContext::*)(qmomkit::MultiIndexView) const;

static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::indices), CqmomIndices>);
static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::index), CqmomIndexAt>);
static_assert(std::is_same_v<decltype(&qmomkit::CqmomSourceContext::moment), CqmomMomentAt>);

// The contexts refer to the caller's inputs; binding a temporary to any of them
// must not compile.
static_assert(std::is_constructible_v<qmomkit::CqmomSourceContext, double, const qmomkit::PhaseSpaceHierarchy&, const qmomkit::MultiIndexTable&,
                                      const std::vector<double>&>);
static_assert(!std::is_constructible_v<qmomkit::CqmomSourceContext, double, qmomkit::PhaseSpaceHierarchy&&, const qmomkit::MultiIndexTable&,
                                       const std::vector<double>&>);
static_assert(!std::is_constructible_v<qmomkit::CqmomSourceContext, double, const qmomkit::PhaseSpaceHierarchy&, qmomkit::MultiIndexTable&&,
                                       const std::vector<double>&>);
static_assert(!std::is_constructible_v<qmomkit::CqmomSourceContext, double, const qmomkit::PhaseSpaceHierarchy&, const qmomkit::MultiIndexTable&,
                                       std::vector<double>&&>);
static_assert(
    !std::is_constructible_v<qmomkit::CqmomSourceContext, double, qmomkit::PhaseSpaceHierarchy&&, qmomkit::MultiIndexTable&&, std::vector<double>&&>);
static_assert(std::is_constructible_v<qmomkit::QmomSourceContext, double, const std::vector<double>&>);
static_assert(!std::is_constructible_v<qmomkit::QmomSourceContext, double, std::vector<double>&&>);
