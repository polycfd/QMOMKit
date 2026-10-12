// The moment and quadrature writers, checked by reading their files back.
//
// What is verified
//   That MomentWriter and QuadratureWriter write a time row whose header names
//   the right columns in the right order and whose values are the ones handed
//   to the writer: QMOM moments and quadrature, CQMOM quadrature, a CqmomState
//   moment row, and requested CQMOM moments in their requested order.
//
// How the reference is constructed
//   Each row is parsed back and compared field by field with the values handed
//   to the writer, a round trip through the written text. The requested CQMOM
//   moments, which the writer reconstructs from the quadrature, are compared
//   instead with the moments of the analytical product rule behind
//   independent_product.qmom, evaluated here.
//
// Error metric and pass criteria
//   Header labels and indices must match exactly, values to 1.0e-13 absolute.
//   The round trip is observed to be exact, and the reconstructed moments
//   differ from the analytical ones by about 1e-15, so the bound is a safety
//   bound.
//
// Note on file output
//   The files written here are the subject under test, not diagnostics, so they
//   are produced on every run and read back. Like the other studies this one
//   prints nothing to the console on success; failures are reported through
//   QMOMKIT_VERIFY_* and an exit code.
//
// Reproduce
//   ctest -R verification_output
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/output.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

constexpr double output_time = 0.25;
constexpr double tolerance = 1.0e-13;

struct ParsedOutput
{
    std::vector<std::string> header_tokens;
    std::vector<std::vector<double>> rows;
};

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::QmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

std::vector<std::string> tokenize_header(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> tokens;
    std::string token;
    while (stream >> token)
    {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<double> tokenize_data(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<double> values;
    double value = 0.0;
    while (stream >> value)
    {
        values.push_back(value);
    }
    return values;
}

ParsedOutput read_output_table(const std::filesystem::path& path)
{
    QMOMKIT_VERIFY_CHECK(std::filesystem::exists(path));

    std::ifstream input(path);
    QMOMKIT_VERIFY_CHECK(input.good());

    std::string header_line;
    QMOMKIT_VERIFY_CHECK(static_cast<bool>(std::getline(input, header_line)));

    // Output files are documented as whitespace-delimited text with a commented
    // header, not CSV. CQMOM labels contain commas inside multi-indices, so the
    // row-level check focuses on data separation rather than label spelling.
    // Header fields may be centered within fixed-width columns.
    const std::size_t first_header_character = header_line.find_first_not_of(" \t");
    QMOMKIT_VERIFY_CHECK(first_header_character != std::string::npos);
    QMOMKIT_VERIFY_CHECK(header_line.compare(first_header_character, 6, "# time") == 0);

    ParsedOutput output{tokenize_header(header_line), {}};
    std::string data_line;
    while (std::getline(input, data_line))
    {
        QMOMKIT_VERIFY_CHECK(data_line.find(',') == std::string::npos);
        output.rows.push_back(tokenize_data(data_line));
        QMOMKIT_VERIFY_CHECK(output.header_tokens.size() == output.rows.back().size() + 1);
    }
    QMOMKIT_VERIFY_CHECK(!output.rows.empty());
    QMOMKIT_VERIFY_CHECK(output.header_tokens.size() >= 2);
    QMOMKIT_VERIFY_CHECK(output.header_tokens[0] == "#");
    QMOMKIT_VERIFY_CHECK(output.header_tokens[1] == "time");
    return output;
}

void verify_header_labels(const ParsedOutput& output, const std::vector<std::string>& labels)
{
    QMOMKIT_VERIFY_CHECK(output.header_tokens.size() == labels.size() + 2);
    for (std::size_t label = 0; label < labels.size(); ++label)
    {
        QMOMKIT_VERIFY_CHECK(output.header_tokens[label + 2] == labels[label]);
    }
}

std::filesystem::path output_path(const qmomkit::OutputOptions& options, const std::string& stem, const std::string& kind)
{
    return options.directory / (stem + "_" + kind + ".txt");
}

void remove_known_outputs(const std::filesystem::path& directory)
{
    for (const char* const filename : {"qmom_moments.txt", "qmom_quadrature.txt", "cqmom_moments.txt", "cqmom_quadrature.txt", "cqmom_state.txt"})
    {
        std::filesystem::remove(directory / filename);
    }
}

void write_qmom_outputs(const qmomkit::OutputOptions& options, const std::string& stem, const qmomkit::QmomState& state,
                        const qmomkit::QmomResult& result)
{
    // The writers read only the generalized moment orders of `OutputOptions`;
    // its directory places the files here.
    auto moment_writer = qmomkit::MomentWriter::for_qmom(output_path(options, stem, "moments"), state, options);
    moment_writer.write(output_time, state, result);
    auto quadrature_writer = qmomkit::QuadratureWriter::for_qmom(output_path(options, stem, "quadrature"), result.weights.size());
    quadrature_writer.write(output_time, result);
}

void write_cqmom_outputs(const qmomkit::OutputOptions& options, const std::string& stem, const qmomkit::CqmomResult& result)
{
    // Requested CQMOM moment output reconstructs exactly the requested raw
    // moments from terminal quadrature. The non-monotonic request makes the
    // output ordering contract visible in the file header and row.
    const qmomkit::RawMomentRequest requested_moments{{{1, 1}, {0, 0}, {2, 0}, {0, 1}}};
    auto moment_writer = qmomkit::MomentWriter::for_cqmom(output_path(options, stem, "moments"), requested_moments);
    moment_writer.write(output_time, result.quadrature);
    auto quadrature_writer = qmomkit::QuadratureWriter::for_cqmom(output_path(options, stem, "quadrature"), 2);
    quadrature_writer.write(output_time, result.quadrature);
}

void write_cqmom_state_output(const std::filesystem::path& path, const qmomkit::CqmomState& state, const qmomkit::CqmomQuadrature& quadrature)
{
    qmomkit::OutputOptions output_options;
    qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, state, output_options);
    writer.write(output_time, state, quadrature);
}

double analytical_product_moment(const qmomkit::MultiIndex& index)
{
    // Independent reference for independent_product.qmom. The x marginal has
    // nodes {1, 3} with weights {0.4, 0.6}; the y marginal has nodes {-1, 2}
    // with weights {0.25, 0.75}. Independence makes terminal weights products.
    struct Node
    {
        double weight;
        double x;
        double y;
    };
    const std::vector<Node> nodes{{0.4 * 0.25, 1.0, -1.0}, {0.4 * 0.75, 1.0, 2.0}, {0.6 * 0.25, 3.0, -1.0}, {0.6 * 0.75, 3.0, 2.0}};

    double moment = 0.0;
    for (const Node& node : nodes)
    {
        moment += node.weight * std::pow(node.x, static_cast<double>(index[0])) * std::pow(node.y, static_cast<double>(index[1]));
    }
    return moment;
}

void verify_qmom_moment_output(const std::filesystem::path& path, const qmomkit::MomentSet& moments)
{
    const ParsedOutput output = read_output_table(path);
    QMOMKIT_VERIFY_CHECK(output.rows.size() == 1);
    QMOMKIT_VERIFY_NEAR(output.rows[0][0], output_time, tolerance);

    std::vector<std::string> labels;
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        labels.push_back("M[" + std::to_string(order) + "]");
    }
    verify_header_labels(output, labels);

    // The moment row is expected to be the configured raw moment sequence in
    // public MomentSet order, preceded by the caller-supplied time.
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        QMOMKIT_VERIFY_NEAR(output.rows[0][order + 1], moments.values()[order], tolerance);
    }
}

void verify_qmom_quadrature_output(const std::filesystem::path& path, const qmomkit::QmomResult& result)
{
    const ParsedOutput output = read_output_table(path);
    const std::size_t node_count = result.weights.size();

    verify_header_labels(output, {"node_index", "weight", "x0"});
    QMOMKIT_VERIFY_CHECK(output.rows.size() == node_count);

    // QMOM quadrature output is active-node oriented: for each public node,
    // write one row with its node index, weight, and abscissa.
    for (std::size_t node = 0; node < node_count; ++node)
    {
        QMOMKIT_VERIFY_NEAR(output.rows[node][0], output_time, tolerance);
        QMOMKIT_VERIFY_NEAR(output.rows[node][1], static_cast<double>(node), 0.0);
        QMOMKIT_VERIFY_NEAR(output.rows[node][2], result.weights[node], tolerance);
        QMOMKIT_VERIFY_NEAR(output.rows[node][3], result.abscissas[node], tolerance);
    }
}

void verify_cqmom_moment_output(const std::filesystem::path& path)
{
    const ParsedOutput output = read_output_table(path);
    QMOMKIT_VERIFY_CHECK(output.rows.size() == 1);
    QMOMKIT_VERIFY_NEAR(output.rows[0][0], output_time, tolerance);
    const std::vector<qmomkit::MultiIndex> requested_moments{{1, 1}, {0, 0}, {2, 0}, {0, 1}};
    verify_header_labels(output, {"M[1,1]", "M[0,0]", "M[2,0]", "M[0,1]"});

    // The expected values are evaluated directly from the analytical product
    // quadrature, independent of QMOMKit's reconstruction routines.
    for (std::size_t moment = 0; moment < requested_moments.size(); ++moment)
    {
        QMOMKIT_VERIFY_NEAR(output.rows[0][moment + 1], analytical_product_moment(requested_moments[moment]), tolerance);
    }
}

void verify_cqmom_state_output(const std::filesystem::path& path, const qmomkit::CqmomState& state)
{
    const ParsedOutput output = read_output_table(path);
    QMOMKIT_VERIFY_CHECK(output.rows.size() == 1);
    QMOMKIT_VERIFY_NEAR(output.rows[0][0], output_time, tolerance);

    const std::vector<qmomkit::MultiIndex> indices = state.moment_indices();
    QMOMKIT_VERIFY_CHECK(indices.size() == state.moment_values().size());
    std::vector<std::string> labels;
    labels.reserve(indices.size());
    for (const qmomkit::MultiIndex& index : indices)
    {
        std::ostringstream label;
        label << "M[";
        for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
        {
            if (dimension != 0)
            {
                label << ',';
            }
            label << index[dimension];
        }
        label << ']';
        labels.push_back(label.str());
    }
    verify_header_labels(output, labels);
    for (std::size_t moment = 0; moment < state.moment_values().size(); ++moment)
    {
        QMOMKIT_VERIFY_NEAR(output.rows[0][moment + 1], state.moment_values()[moment], tolerance);
    }
}

void verify_cqmom_quadrature_output(const std::filesystem::path& path, const qmomkit::CqmomQuadrature& quadrature)
{
    const ParsedOutput output = read_output_table(path);
    const std::size_t node_count = quadrature.weights().size();

    verify_header_labels(output, {"node_index", "path_0", "path_1", "weight", "x0", "x1"});
    QMOMKIT_VERIFY_CHECK(output.rows.size() == node_count);

    // CQMOM quadrature output is terminal-node oriented: for each public
    // terminal node, write one row with its terminal path, weight, and
    // coordinates in phase-space dimension order.
    for (std::size_t node = 0; node < node_count; ++node)
    {
        QMOMKIT_VERIFY_NEAR(output.rows[node][0], output_time, tolerance);
        QMOMKIT_VERIFY_NEAR(output.rows[node][1], static_cast<double>(node), 0.0);
        QMOMKIT_VERIFY_NEAR(output.rows[node][2], static_cast<double>(quadrature.terminal_paths()[node].node_index_at_level(0)), 0.0);
        QMOMKIT_VERIFY_NEAR(output.rows[node][3], static_cast<double>(quadrature.terminal_paths()[node].node_index_at_level(1)), 0.0);
        QMOMKIT_VERIFY_NEAR(output.rows[node][4], quadrature.weights()[node], tolerance);
        QMOMKIT_VERIFY_NEAR(output.rows[node][5], quadrature.abscissas()[0][node], tolerance);
        QMOMKIT_VERIFY_NEAR(output.rows[node][6], quadrature.abscissas()[1][node], tolerance);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    // Output verification is intentionally workflow-level: compute public QMOM
    // and CQMOM quadratures, write moments and quadratures with the public
    // writers, and parse the files back as a downstream text consumer would.
    if (argc != 5)
    {
        std::cerr << "Usage: output_workflows_verification <qmom.qmom> <cqmom.qmom> <output-directory> <qmom|cqmom>\n";
        return 2;
    }

    const std::filesystem::path qmom_configuration = argv[1];
    const std::filesystem::path cqmom_configuration = argv[2];
    const std::filesystem::path output_directory = argv[3];
    const std::string selected_case = argv[4];
    if (selected_case != "qmom" && selected_case != "cqmom")
    {
        std::cerr << "Unknown output workflow case '" << selected_case << "'.\n";
        return 2;
    }

    std::filesystem::create_directories(output_directory);
    remove_known_outputs(output_directory);

    const auto configured_qmom = qmomkit::config::load_qmom(qmom_configuration);
    if (!configured_qmom.ok())
    {
        print_configuration_errors(configured_qmom);
        return 1;
    }
    const qmomkit::QmomResult qmom_result = qmomkit::compute_qmom(configured_qmom.value());
    QMOMKIT_VERIFY_CHECK(qmom_result.status == qmomkit::QmomStatus::success);

    const auto configured_cqmom = qmomkit::config::load_cqmom_prepared(cqmom_configuration);
    if (!configured_cqmom.ok())
    {
        print_configuration_errors(configured_cqmom);
        return 1;
    }
    const qmomkit::CqmomResult cqmom_result = qmomkit::compute_cqmom(configured_cqmom.value());
    QMOMKIT_VERIFY_CHECK(cqmom_result.status == qmomkit::CqmomStatus::success);

    const qmomkit::OutputOptions output_options{output_directory};
    if (selected_case == "qmom")
    {
        const qmomkit::QmomState qmom_state = configured_qmom.value().initial_state();
        write_qmom_outputs(output_options, "qmom", qmom_state, qmom_result);
        verify_qmom_moment_output(output_directory / "qmom_moments.txt", qmomkit::MomentSet(qmom_state.moment_values()));
        verify_qmom_quadrature_output(output_directory / "qmom_quadrature.txt", qmom_result);
    }
    else
    {
        write_cqmom_outputs(output_options, "cqmom", cqmom_result);
        verify_cqmom_moment_output(output_directory / "cqmom_moments.txt");
        verify_cqmom_quadrature_output(output_directory / "cqmom_quadrature.txt", cqmom_result.quadrature);

        const qmomkit::CqmomState state = configured_cqmom.value().initial_state();
        write_cqmom_state_output(output_directory / "cqmom_state.txt", state, cqmom_result.quadrature);
        verify_cqmom_state_output(output_directory / "cqmom_state.txt", state);
    }

    return 0;
}
