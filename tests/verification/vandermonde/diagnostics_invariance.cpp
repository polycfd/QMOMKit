// Detailed conditional-moment diagnostics in other units and at other number
// densities.
//
// What is verified
//   That the condition estimate and the relative residual reported at
//   CqmomDiagnosticLevel::detailed describe the system the solvers decide on,
//   so that neither changes with the number density nor with a power-of-two
//   change of the unit of the parent coordinate, and that the warning does not
//   fire for a well-conditioned case in units where the unscaled weighted
//   matrix would exceed the default threshold.
//
// How the reference is constructed
//   Analytically. The 2D product case has parent nodes x = {1, 3}. Divided by
//   the power of two above the largest node, they are {0.25, 0.75} in every
//   power-of-two unit, so the scaled Vandermonde matrix is
//   [[1, 1], [0.25, 0.75]], with infinity norm 2 and inverse
//   [[1.5, -2], [-0.5, 2]], with infinity norm 3.5: the condition number is
//   7. The unscaled weighted matrix [[w1, w2], [w1 x1, w2 x2]] in units
//   2^40 times larger has the closed-form infinity-norm condition number
//   computed below, about 9e12.
//
// Error metric and pass criteria
//   In the original units the condition estimate matches 7 to a relative
//   1e-13, the rounding of the nodes Wheeler reconstructs. Under a
//   power-of-two change of number density or unit, every operation scales
//   exactly, so it equals that value exactly; under another number density,
//   whose Wheeler nodes round differently, to a relative 1e-12. The relative
//   residual is at most 1e-14, and the warning is off.
//
// Reproduce
//   ctest -R "^verification_vandermonde_diagnostics_"
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

constexpr double analytical_condition = 7.0;

// Moment data of the terminal quadrature x = {1, 3} with weights {0.4, 0.6},
// times y = {-1, 2} with weights {0.25, 0.75}, with x in units `x_scale` and
// every weight multiplied by `number_density`.
qmomkit::CqmomMomentData product_moment_data(const double x_scale, const double number_density)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<double> x{1.0 * x_scale, 3.0 * x_scale};
    const std::vector<double> x_weights{0.4, 0.6};
    const std::vector<double> y{-1.0, 2.0};
    const std::vector<double> y_weights{0.25, 0.75};
    std::vector<double> values;
    for (const qmomkit::MultiIndex& index : layout.required_indices.to_multi_indices())
    {
        double moment = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            for (std::size_t j = 0; j < y.size(); ++j)
            {
                moment += number_density * x_weights[i] * y_weights[j] * std::pow(x[i], static_cast<int>(index[0])) *
                          std::pow(y[j], static_cast<int>(index[1]));
            }
        }
        values.push_back(moment);
    }
    return qmomkit::CqmomMomentData{layout, values};
}

// The diagnostics of the root branch, the only branch with a conditional solve.
qmomkit::ConditionalMomentDiagnostics root_diagnostics(const double x_scale, const double number_density)
{
    qmomkit::CqmomCase cqmom_case{product_moment_data(x_scale, number_density)};
    cqmom_case.qmom_options.adaptive = false;
    cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_case);
    QMOMKIT_VERIFY_CHECK(result.ok());
    std::size_t solved_branches = 0;
    qmomkit::ConditionalMomentDiagnostics root;
    for (const qmomkit::CqmomDiagnostics& branch : result.diagnostics)
    {
        if (branch.conditional_moment_diagnostics.condition_estimate.has_value())
        {
            ++solved_branches;
            root = branch.conditional_moment_diagnostics;
        }
    }
    QMOMKIT_VERIFY_CHECK(solved_branches == 1);
    return root;
}

// The diagnostics in the original units, checked against the analytical value.
double reference_condition()
{
    const qmomkit::ConditionalMomentDiagnostics reference = root_diagnostics(1.0, 1.0);
    QMOMKIT_VERIFY_CHECK(reference.condition_estimate.has_value());
    QMOMKIT_VERIFY_NEAR(*reference.condition_estimate, analytical_condition, 1.0e-13 * analytical_condition);
    return *reference.condition_estimate;
}

void check_diagnostics(const qmomkit::ConditionalMomentDiagnostics& diagnostics, const double expected_condition, const double relative_tolerance)
{
    QMOMKIT_VERIFY_CHECK(diagnostics.condition_estimate.has_value() && diagnostics.residual_norm.has_value());
    QMOMKIT_VERIFY_NEAR(*diagnostics.condition_estimate, expected_condition, relative_tolerance * expected_condition);
    QMOMKIT_VERIFY_CHECK(*diagnostics.residual_norm <= 1.0e-14);
    QMOMKIT_VERIFY_CHECK(!diagnostics.conditioning_warning);
}

void verify_number_density_does_not_change_diagnostics()
{
    const double reference = reference_condition();
    for (const double number_density : {1.0, std::ldexp(1.0, -60), std::ldexp(1.0, 60)})
    {
        check_diagnostics(root_diagnostics(1.0, number_density), reference, 0.0);
    }
    check_diagnostics(root_diagnostics(1.0, 3.7e-9), reference, 1.0e-12);
}

void verify_power_of_two_unit_does_not_change_diagnostics()
{
    const double reference = reference_condition();
    const double x_scale = std::ldexp(1.0, 40);
    check_diagnostics(root_diagnostics(x_scale, 1.0), reference, 0.0);
    check_diagnostics(root_diagnostics(std::ldexp(1.0, -40), 1.0), reference, 0.0);

    // The unscaled weighted matrix in these units is far above the default
    // warning threshold, although the system is as well conditioned as in the
    // original units.
    const double w1 = 0.4;
    const double w2 = 0.6;
    const double x1 = 1.0 * x_scale;
    const double x2 = 3.0 * x_scale;
    const double determinant = w1 * w2 * (x2 - x1);
    const double matrix_norm = std::max(w1 + w2, w1 * x1 + w2 * x2);
    const double inverse_norm = std::max(w2 * x2 + w2, w1 * x1 + w1) / determinant;
    QMOMKIT_VERIFY_CHECK(matrix_norm * inverse_norm > qmomkit::ConditionalMomentOptions{}.conditioning_warning_threshold);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"number_density_does_not_change_diagnostics", verify_number_density_does_not_change_diagnostics},
            {"power_of_two_unit_does_not_change_diagnostics", verify_power_of_two_unit_does_not_change_diagnostics},
        });
}
