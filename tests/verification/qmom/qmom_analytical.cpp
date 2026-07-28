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
    std::cout << "QMOMKit one-dimensional QMOM analytical verification\n";
    std::cout << "Configuration: " << configuration_path << '\n';
    std::cout << "Tolerance: " << std::scientific << qmomkit::verification::qmom_analytical_tolerance << " (absolute and relative)\n\n";
    std::cout << std::left << std::setw(24) << "case" << std::right << std::setw(14) << "weight abs" << std::setw(14) << "node abs" << std::setw(14)
              << "moment abs" << std::setw(14) << "max abs" << std::setw(14) << "max rel" << std::setw(9) << "status" << '\n';
    std::cout << std::string(103, '-') << '\n';
    std::cout << std::left << std::setw(24) << result.name << std::right << std::scientific << std::setprecision(4) << std::setw(14)
              << result.weight_errors.maximum_absolute << std::setw(14) << result.abscissa_errors.maximum_absolute << std::setw(14)
              << result.moment_errors.maximum_absolute << std::setw(14) << result.maximum_absolute_error() << std::setw(14)
              << result.maximum_relative_error() << std::setw(9) << (result.passed() ? "PASS" : "FAIL") << '\n';
    if (!result.passed() && !result.message.empty())
    {
        std::cerr << result.message << '\n';
    }
    return result.passed() ? 0 : 1;
}
