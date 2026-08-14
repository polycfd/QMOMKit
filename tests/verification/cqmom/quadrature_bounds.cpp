#include "../verification_helpers.hpp"
#include "qmomkit/cqmom.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& nodes)
{
    std::vector<double> moments(4, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

qmomkit::CqmomDenseCase bounded_case(const std::vector<double>& moments, const qmomkit::BoundedInterval interval)
{
    qmomkit::CqmomDenseCase case_{qmomkit::PhaseSpaceHierarchy({2}), {{4}, moments}};
    case_.coordinate_support = {interval};
    return case_;
}

void verify_accepted_nodes_inside_interval()
{
    const std::vector<double> weights{0.4, 0.6};
    const std::vector<double> nodes{2.25, 4.5};
    const qmomkit::BoundedInterval interval{2.0, 5.0};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(bounded_case(raw_moments(weights, nodes), interval));

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights().size() == weights.size());
    QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas().size() == 1);
    QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[0].size() == nodes.size());
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(result.quadrature.weights()[node], weights[node], tolerance);
        QMOMKIT_VERIFY_NEAR(result.quadrature.abscissas()[0][node], nodes[node], tolerance);
        QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[0][node] >= interval.lower);
        QMOMKIT_VERIFY_CHECK(result.quadrature.abscissas()[0][node] <= interval.upper);
    }
}

void verify_rejected_nodes_outside_interval()
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<double> nodes{0.25, 1.25};
    const qmomkit::BoundedInterval interval{0.0, 1.0};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(bounded_case(raw_moments(weights, nodes), interval));

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.size() == 1);
    const qmomkit::CqmomDiagnostics& diagnostic = result.diagnostics.front();
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_support_checked);
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_nodes.size() == nodes.size());
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_weights.size() == weights.size());
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_nodes[node], nodes[node], tolerance);
        QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_weights[node], weights[node], tolerance);
    }
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_nodes[1] > interval.upper);
}

} // namespace

int main()
{
    verify_accepted_nodes_inside_interval();
    verify_rejected_nodes_outside_interval();
    return 0;
}
