// Public configuration loaders rejecting malformed or incomplete .qmom files.
//
// What is verified
//   That each faulty configuration is refused with a specific error, and that
//   the loader returns no case at all rather than a partially constructed one
//   a caller might then use.
//
// How the reference is constructed
//   By construction rather than by computation: each case is a .qmom file
//   containing one deliberate fault - a malformed section header, an invalid
//   boolean, missing moments, a hierarchy the loader cannot infer, an unknown
//   Vandermonde solver, a moment set that cannot support the requested node
//   count. The expected error is written beside each case. There is no
//   numerical reference because no case should be produced.
//
// Pass criteria
//   The load must fail, the reported error must identify the intended cause,
//   and the result must hold no case. Checking failure alone would not
//   distinguish a correct diagnosis from an unrelated one.
//
// Reproduce
//   ctest -R verification_configuration
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/configuration.hpp"

#include <cctype>
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

template <class T> bool has_reported_line(const qmomkit::config::ConfigurationResult<T>& result, const std::filesystem::path& path)
{
    // A located error reads "<file>:<line>: <message>"; without a line it reads
    // "<file>: <message>". Require the file name, a colon, at least one digit,
    // and a second colon, so the line number itself must be present.
    const std::string message = result.error();
    const std::string location_prefix = path.string() + ':';
    const std::size_t prefix_position = message.find(location_prefix);
    if (prefix_position == std::string::npos)
    {
        return false;
    }
    const std::size_t digits_begin = prefix_position + location_prefix.size();
    std::size_t digits_end = digits_begin;
    while (digits_end < message.size() && std::isdigit(static_cast<unsigned char>(message[digits_end])) != 0)
    {
        ++digits_end;
    }
    return digits_end > digits_begin && digits_end < message.size() && message[digits_end] == ':';
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
        QMOMKIT_VERIFY_CHECK(has_reported_line(result, path));
    }
}

void verify_case(const ExpectedFailure& expected, const std::filesystem::path& directory)
{
    const std::filesystem::path path = directory / expected.filename;
    QMOMKIT_VERIFY_CHECK(std::filesystem::exists(path));

    if (expected.loader == Loader::qmom)
    {
        verify_failed_load(qmomkit::config::load_qmom(path), expected, path);
    }
    else
    {
        verify_failed_load(qmomkit::config::load_cqmom_prepared(path), expected, path);
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
        {"qmom_missing_moments.qmom", Loader::qmom, "QMOM configuration requires either '[moments] values'", false},

        // Syntax is valid here, but the requested two-node inversion does not
        // have the required four raw moments. This should fail at case creation
        // time, before compute_qmom(...) is called.
        {"qmom_insufficient_moments.qmom", Loader::qmom, "requires at least 2 * node_count moments.", false},

        // Vandermonde solver selection is a user-facing CQMOM option. Unknown
        // names should not silently fall back to the default policy.
        {"cqmom_invalid_vandermonde_solver.qmom", Loader::cqmom, "vandermonde.solver must be"},

        // The former branch-count switch is now [qmom] adaptive; old files name
        // the replacement instead of failing on an unknown key.
        {"cqmom_removed_adaptive_branch_counts.qmom", Loader::cqmom, "phase_space.adaptive_branch_counts was removed; use [qmom] adaptive."},

        // The public CQMOM loader cannot infer a phase-space hierarchy. Missing
        // node counts should be reported as missing required configuration.
        {"cqmom_missing_phase_space.qmom", Loader::cqmom, "Missing required key 'phase_space.node_counts'.", false},

        // A conditioning order must be a permutation of dimensions. Duplicate
        // dimensions describe an invalid mathematical CQMOM hierarchy.
        {"cqmom_invalid_conditioning_order.qmom", Loader::cqmom, "conditioning_order must be a permutation"},

        // The full tensor of raw moments must be complete. Missing tensor
        // entries are rejected before the CQMOM case is exposed to the caller.
        {"cqmom_incomplete_tensor.qmom", Loader::cqmom, "The full tensor of raw moments is incomplete; it is missing 1 of its 4 entries.", false},
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

    return 0;
}
