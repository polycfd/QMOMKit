#include "../verification_helpers.hpp"
#include "qmomkit/qmom.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& nodes)
{
    std::vector<double> moments(2 * nodes.size(), 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
    }
    return moments;
}

void verify_fixed_bounded_qmom_accepts_nodes_inside_interval()
{
    const std::vector<double> weights{0.4, 0.6};
    const std::vector<double> nodes{2.25, 4.5};
    const qmomkit::MomentSet moments(raw_moments(weights, nodes));
    const qmomkit::AbscissaBounds bounds{2.0, 5.0};

    const qmomkit::QmomResult bounded = qmomkit::compute_qmom(moments, 2, {}, bounds);
    const qmomkit::QmomResult unbounded = qmomkit::compute_qmom(moments, 2);

    QMOMKIT_VERIFY_CHECK(bounded.ok());
    QMOMKIT_VERIFY_CHECK(unbounded.ok());
    QMOMKIT_VERIFY_CHECK(bounded.abscissas == unbounded.abscissas);
    QMOMKIT_VERIFY_CHECK(bounded.weights == unbounded.weights);
}

void verify_fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping()
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<double> nodes{0.25, 1.25};
    const qmomkit::QmomResult result =
        qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments(weights, nodes)), 2, {}, qmomkit::AbscissaBounds{0.0, 1.0});

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::QmomStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(result.abscissas.size() == nodes.size());
    QMOMKIT_VERIFY_NEAR(result.abscissas[1], nodes[1], tolerance);
    QMOMKIT_VERIFY_CHECK(result.abscissas[1] > 1.0);
}

void verify_adaptive_bounded_qmom_uses_highest_admissible_count()
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<double> nodes{0.25, 1.25};
    qmomkit::QmomOptions options;
    options.adaptive = true;

    const qmomkit::QmomResult result =
        qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments(weights, nodes)), 2, options, qmomkit::AbscissaBounds{0.0, 1.0});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.weights.size() == 1);
    QMOMKIT_VERIFY_CHECK(result.abscissas.size() == 1);
    QMOMKIT_VERIFY_NEAR(result.weights[0], 1.0, tolerance);
    QMOMKIT_VERIFY_NEAR(result.abscissas[0], 0.75, tolerance);
}

} // namespace

int main()
{
    verify_fixed_bounded_qmom_accepts_nodes_inside_interval();
    verify_fixed_bounded_qmom_rejects_nodes_outside_interval_without_clipping();
    verify_adaptive_bounded_qmom_uses_highest_admissible_count();
    return 0;
}
