// The adaptive QMOM path returning fewer nodes than requested when the data
// supports fewer.
//
// What is verified
//   That a moment sequence generated from two points is inverted to exactly
//   those two points even though three were requested, that the weights sum to
//   the zeroth moment, and that all six supplied moments are reproduced from
//   the two-node result.
//
// How the reference is constructed
//   Chosen, not computed: two positions and weights are picked and the six
//   moments are built from them. A two-point set is exactly determined by its
//   first four moments, so the remaining two are reproduced only if the
//   inversion found the right pair rather than any pair fitting a subset.
//
// Error metric and pass criteria
//   Maximum absolute deviation of positions, weights and rebuilt moments from
//   the chosen values, at 1.0e-12. The node count must be two, not three: a
//   third node would have to carry zero or negative weight.
//
// Reproduce
//   ctest -R verification_adaptive_public_qmom_reconstruction
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

std::vector<double> raw_moments(const std::vector<double>& weights, const std::vector<double>& abscissas, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], static_cast<int>(order));
        }
    }
    return moments;
}

void verify_adaptive_public_qmom_reconstruction()
{
    // The six supplied moments are independently constructed from this
    // two-node rule, although the public case requests up to three nodes.
    const std::vector<double> expected_weights{0.4, 0.6};
    const std::vector<double> expected_abscissas{-0.75, 1.25};
    const std::vector<double> moments = raw_moments(expected_weights, expected_abscissas, 6);

    qmomkit::QmomOptions options;
    options.adaptive = true;
    options.rmin = 1.0e-10;
    options.eabs = 1.0e-12;
    options.cutoff = 1.0e-30;
    const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::QmomCase{qmomkit::MomentSet(moments), 3, options});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.weights.size() == expected_weights.size());
    QMOMKIT_VERIFY_CHECK(result.abscissas.size() == expected_abscissas.size());
    for (std::size_t node = 0; node < expected_weights.size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(result.weights[node], expected_weights[node], tolerance);
        QMOMKIT_VERIFY_NEAR(result.abscissas[node], expected_abscissas[node], tolerance);
    }

    const qmomkit::Quadrature1D quadrature{result.abscissas, result.weights};
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(quadrature, moments.size());
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        QMOMKIT_VERIFY_NEAR(reconstructed[order], moments[order], tolerance);
    }
    QMOMKIT_VERIFY_NEAR(result.weights[0] + result.weights[1], moments[0], tolerance);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"adaptive_public_qmom_reconstruction", verify_adaptive_public_qmom_reconstruction},
                                             });
}
