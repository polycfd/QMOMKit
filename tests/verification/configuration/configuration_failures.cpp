#include "../verification_helpers.hpp"
#include "qmomkit/configuration.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace
{

enum class Loader
{
    qmom,
    cqmom
};

struct ExpectedFailure
{
    std::string filename;
    Loader loader = Loader::qmom;
    std::string message_fragment;
    bool line_required = true;
};

template <class T> bool has_message_fragment(const qmomkit::config::ConfigurationResult<T>& result, const std::string& fragment)
{
    return result.error().find(fragment) != std::string::npos;
}

template <class T> bool has_source_name(const qmomkit::config::ConfigurationResult<T>& result, const std::filesystem::path& path)
{
    return result.error().find(path.string()) != std::string::npos;
}

template <class T> bool has_reported_line(const qmomkit::config::ConfigurationResult<T>& result)
{
    return result.error().find(':') != std::string::npos;
}

template <class T>
void verify_failed_load(const qmomkit::config::ConfigurationResult<T>& result, const ExpectedFailure& expected, const std::filesystem::path& path)
{
    // Public loaders must not hand a user a partially constructed case when
    // configuration validation fails. The diagnostic string is the supported
    // way for applications to explain the bad file to the user.
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(!result);
    QMOMKIT_VERIFY_CHECK(!result.error().empty());
    if (!has_message_fragment(result, expected.message_fragment))
    {
        std::cerr << "Expected message fragment: " << expected.message_fragment << '\n';
        std::cerr << "  returned: " << result.error() << '\n';
    }
    QMOMKIT_VERIFY_CHECK(has_message_fragment(result, expected.message_fragment));
    QMOMKIT_VERIFY_CHECK(has_source_name(result, path));
    if (expected.line_required)
    {
        QMOMKIT_VERIFY_CHECK(has_reported_line(result));
    }
}

void verify_case(const ExpectedFailure& expected, const std::filesystem::path& directory)
{
    const std::filesystem::path path = directory / expected.filename;
    std::cerr << "Checking " << expected.filename << '\n';
    QMOMKIT_VERIFY_CHECK(std::filesystem::exists(path));

    if (expected.loader == Loader::qmom)
    {
        verify_failed_load(qmomkit::config::load_qmom(path), expected, path);
    }
    else
    {
        verify_failed_load(qmomkit::config::load_cqmom(path), expected, path);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc < 2 || argc > 2)
    {
        std::cerr << "Usage: configuration_failures_verification <failure-case.qmom>\n";
        return 2;
    }

    const std::filesystem::path case_directory = std::filesystem::path(argv[1]).parent_path();

    const std::vector<ExpectedFailure> failures{
        // A malformed section header is the representative syntax error. Users
        // should see a syntax failure tied to the offending file and line.
        {"malformed_syntax.qmom", Loader::qmom, "Malformed section header."},

        // Invalid booleans are common copy/edit mistakes. They must be rejected
        // instead of interpreted with truthy/falsy conventions.
        {"qmom_invalid_adaptive.qmom", Loader::qmom, "qmom.adaptive must be true or false."},

        // A QMOM case needs either explicit moments or an analytical
        // initialization. Without one, no meaningful case can be created.
        {"qmom_missing_moments.qmom", Loader::qmom, "Missing required key 'initial.type'", false},

        // Syntax is valid here, but the requested two-node inversion does not
        // have the required four raw moments. This should fail at case creation
        // time, before compute_qmom(...) is called.
        {"qmom_insufficient_moments.qmom", Loader::qmom, "requires at least 2 * node_count moments.", false},

        // Vandermonde solver selection is a user-facing CQMOM option. Unknown
        // names should not silently fall back to the default policy.
        {"cqmom_invalid_vandermonde_solver.qmom", Loader::cqmom, "vandermonde.solver must be"},

        // CQMOM branch-count adaptivity is a boolean workflow control separate
        // from one-dimensional Wheeler adaptivity.
        {"cqmom_invalid_adaptive_branch_counts.qmom", Loader::cqmom, "adaptive_branch_counts must be true or false."},

        // The public CQMOM loader cannot infer a phase-space hierarchy. Missing
        // node counts should be reported as missing required configuration.
        {"cqmom_missing_phase_space.qmom", Loader::cqmom, "Missing required key 'phase_space.node_counts'.", false},

        // A conditioning order must be a permutation of dimensions. Duplicate
        // dimensions describe an invalid mathematical CQMOM hierarchy.
        {"cqmom_invalid_conditioning_order.qmom", Loader::cqmom, "conditioning_order must be a permutation"},

        // Dense tensor moment input must be complete. Missing tensor entries
        // are rejected before the CQMOM case is exposed to the caller.
        {"cqmom_incomplete_tensor.qmom", Loader::cqmom, "moment tensor is incomplete", false},
    };

    bool matched = false;
    const std::string selected_filename = std::filesystem::path(argv[1]).filename().string();
    for (const ExpectedFailure& failure : failures)
    {
        if (failure.filename != selected_filename)
        {
            continue;
        }
        matched = true;
        verify_case(failure, case_directory);
    }
    if (!matched)
    {
        std::cerr << "Unknown configuration failure case '" << selected_filename << "'.\n";
        return 2;
    }

    std::cout << "QMOMKit configuration failure-path verification\n";
    std::cout << "Verified " << selected_filename << " through public load_qmom/load_cqmom workflows.\n";
    return 0;
}
