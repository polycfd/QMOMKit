// The public one-dimensional QMOM workflow against classical quadrature rules.
//
// What is verified
//   That loading a .qmom file and calling compute_qmom recovers the positions
//   and weights of a known rule, and that every configured raw moment can be
//   rebuilt from the recovered quadrature. Driving it through the configuration
//   loader rather than the inversion helper is deliberate: it checks the
//   workflow a user actually follows.
//
// How the reference is constructed
//   Not invented for this test. The expected positions and weights in
//   expected.hpp are the classical Gauss rules: Gauss-Legendre for a uniform
//   distribution, Gauss-Hermite for a normal one, Gauss-Laguerre for an
//   exponential one. The three-point normal rule, for instance, is positions
//   -sqrt(3), 0, +sqrt(3) with weights 1/6, 2/3, 1/6, which is published
//   independently of any software. The moments in each .qmom file are the raw
//   moments of that rule, so the inversion is being asked to recover a result
//   that was established long before QMOMKit existed.
//
// Error metric and pass criteria
//   Maximum absolute and relative deviation of positions, weights and rebuilt
//   moments from the reference, at 1.0e-12. An n-point rule is exactly
//   determined by its first 2n moments, so there is a single correct answer
//   rather than a family of acceptable ones.
//
// Reproduce
//   ctest -R verification_qmom_analytical
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "expected.hpp"
#include "qmomkit/configuration.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

int main(const int argc, char** argv)
{
    // Driver for the configured one-dimensional analytical suite. Each .qmom
    // file supplies raw moments generated from the named independent reference
    // in expected.hpp, so this executable verifies the public configuration
    // and compute_qmom workflow rather than only a helper-level inversion.
    if (argc != 3)
    {
        std::cerr << "Usage: qmom_analytical_verification <configuration.qmom> <case_identifier>\n";
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const std::string case_identifier = argv[2];
    const qmomkit::verification::QmomAnalyticalCase* const analytical_case = qmomkit::verification::find_qmom_analytical_case(case_identifier);
    if (analytical_case == nullptr)
    {
        std::cerr << "Unknown analytical case identifier '" << case_identifier << "'.\n";
        return 2;
    }

    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    const qmomkit::verification::QmomVerificationResult result = qmomkit::verification::verify_qmom_case(*analytical_case, configured.value());
    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!(result.passed()))
    {
        std::cerr << "QMOMKit one-dimensional QMOM analytical verification\n";
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Tolerance: " << std::scientific << qmomkit::verification::qmom_analytical_tolerance << " (absolute and relative)\n\n";
        std::cerr << std::left << std::setw(24) << "case" << std::right << std::setw(14) << "weight abs" << std::setw(14) << "node abs"
                  << std::setw(14) << "moment abs" << std::setw(14) << "max abs" << std::setw(14) << "max rel" << std::setw(9) << "status" << '\n';
        std::cerr << std::string(103, '-') << '\n';
        std::cerr << std::left << std::setw(24) << result.name << std::right << std::scientific << std::setprecision(4) << std::setw(14)
                  << result.weight_errors.maximum_absolute << std::setw(14) << result.abscissa_errors.maximum_absolute << std::setw(14)
                  << result.moment_errors.maximum_absolute << std::setw(14) << result.maximum_absolute_error() << std::setw(14)
                  << result.maximum_relative_error() << std::setw(9) << (result.passed() ? "PASS" : "FAIL") << '\n';
        if (!result.passed() && !result.message.empty())
        {
            std::cerr << result.message << '\n';
        }
    }

    return result.passed() ? 0 : 1;
}
