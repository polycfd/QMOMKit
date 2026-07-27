#include "../verification_helpers.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
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
    std::vector<double> values;
};

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::QmomCase>& configured)
{
    for (const qmomkit::config::ConfigurationError& error : configured.errors)
    {
        std::cerr << error.source_name;
        if (error.line != 0)
        {
            std::cerr << ':' << error.line;
        }
        std::cerr << ": " << error.message << '\n';
    }
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>& configured)
{
    for (const qmomkit::config::ConfigurationError& error : configured.errors)
    {
        std::cerr << error.source_name;
        if (error.line != 0)
        {
            std::cerr << ':' << error.line;
        }
        std::cerr << ": " << error.message << '\n';
    }
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
    std::string data_line;
    QMOMKIT_VERIFY_CHECK(static_cast<bool>(std::getline(input, header_line)));
    QMOMKIT_VERIFY_CHECK(static_cast<bool>(std::getline(input, data_line)));

    // Output files are documented as whitespace-delimited text with a commented
    // header, not CSV. CQMOM labels contain commas inside multi-indices, so the
    // row-level check focuses on data separation rather than label spelling.
    QMOMKIT_VERIFY_CHECK(header_line.rfind("# ", 0) == 0);
    QMOMKIT_VERIFY_CHECK(data_line.find(',') == std::string::npos);

    std::string trailing_line;
    QMOMKIT_VERIFY_CHECK(!static_cast<bool>(std::getline(input, trailing_line)));

    ParsedOutput output{tokenize_header(header_line), tokenize_data(data_line)};
    QMOMKIT_VERIFY_CHECK(output.header_tokens.size() == output.values.size() + 1);
    QMOMKIT_VERIFY_CHECK(output.header_tokens.size() >= 2);
    QMOMKIT_VERIFY_CHECK(output.header_tokens[0] == "#");
    QMOMKIT_VERIFY_CHECK(output.header_tokens[1] == "time");
    QMOMKIT_VERIFY_NEAR(output.values[0], output_time, tolerance);
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
    for (const std::string& filename :
         {"qmom_moments.txt", "qmom_quadrature.txt", "cqmom_moments.txt", "cqmom_quadrature.txt", "disabled_moments.txt", "disabled_quadrature.txt"})
    {
        std::filesystem::remove(directory / filename);
    }
}

void write_qmom_outputs(const qmomkit::OutputOptions& options, const std::string& stem, const qmomkit::QmomCase& qmom_case,
                        const qmomkit::QmomResult& result)
{
    // `OutputOptions` is a driver policy object. This representative workflow
    // shows how a user enables or suppresses file creation before constructing
    // the public writers.
    if (options.write_moments)
    {
        auto writer = qmomkit::MomentWriter::for_qmom(output_path(options, stem, "moments"), qmom_case.moments.size());
        writer.write(output_time, qmom_case.moments);
    }
    if (options.write_quadrature)
    {
        auto writer = qmomkit::QuadratureWriter::for_qmom(output_path(options, stem, "quadrature"), qmom_case.node_count);
        writer.write(output_time, result);
    }
}

void write_cqmom_outputs(const qmomkit::OutputOptions& options, const std::string& stem, const qmomkit::CqmomResult& result)
{
    // Requested CQMOM moment output reconstructs exactly the requested raw
    // moments from terminal quadrature. The non-monotonic request makes the
    // output ordering contract visible in the file header and row.
    const qmomkit::RawMomentRequest requested_moments{{{1, 1}, {0, 0}, {2, 0}, {0, 1}}};
    if (options.write_moments)
    {
        auto writer = qmomkit::MomentWriter::for_cqmom(output_path(options, stem, "moments"), requested_moments);
        writer.write(output_time, result.quadrature);
    }
    if (options.write_quadrature)
    {
        auto writer = qmomkit::QuadratureWriter::for_cqmom(output_path(options, stem, "quadrature"), 2, result.quadrature.weights().size());
        writer.write(output_time, result.quadrature);
    }
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
        QMOMKIT_VERIFY_NEAR(output.values[order + 1], moments.values()[order], tolerance);
    }
}

void verify_qmom_quadrature_output(const std::filesystem::path& path, const qmomkit::QmomResult& result)
{
    const ParsedOutput output = read_output_table(path);
    const std::size_t node_count = result.weights.size();

    std::vector<std::string> labels;
    for (std::size_t node = 0; node < node_count; ++node)
    {
        labels.push_back("w" + std::to_string(node));
        labels.push_back("x" + std::to_string(node));
    }
    verify_header_labels(output, labels);

    // QuadratureWriter preserves public QMOM node order and writes each
    // weight immediately followed by its abscissa.
    for (std::size_t node = 0; node < node_count; ++node)
    {
        QMOMKIT_VERIFY_NEAR(output.values[1 + 2 * node], result.weights[node], tolerance);
        QMOMKIT_VERIFY_NEAR(output.values[1 + 2 * node + 1], result.abscissas[node], tolerance);
    }
}

void verify_cqmom_moment_output(const std::filesystem::path& path)
{
    const ParsedOutput output = read_output_table(path);
    const std::vector<qmomkit::MultiIndex> requested_moments{{1, 1}, {0, 0}, {2, 0}, {0, 1}};
    verify_header_labels(output, {"M[1,1]", "M[0,0]", "M[2,0]", "M[0,1]"});

    // The expected values are evaluated directly from the analytical product
    // quadrature, independent of QMOMKit's reconstruction routines.
    for (std::size_t moment = 0; moment < requested_moments.size(); ++moment)
    {
        QMOMKIT_VERIFY_NEAR(output.values[moment + 1], analytical_product_moment(requested_moments[moment]), tolerance);
    }
}

void verify_cqmom_quadrature_output(const std::filesystem::path& path, const qmomkit::CqmomQuadrature& quadrature)
{
    const ParsedOutput output = read_output_table(path);
    const std::size_t node_count = quadrature.weights().size();

    std::vector<std::string> labels;
    for (std::size_t node = 0; node < node_count; ++node)
    {
        labels.push_back("w" + std::to_string(node));
        labels.push_back("x0_" + std::to_string(node));
        labels.push_back("x1_" + std::to_string(node));
    }
    verify_header_labels(output, labels);

    // CQMOM quadrature output is terminal-node oriented: for each public
    // terminal node, write its weight followed by coordinates in phase-space
    // dimension order.
    for (std::size_t node = 0; node < node_count; ++node)
    {
        const std::size_t offset = 1 + 3 * node;
        QMOMKIT_VERIFY_NEAR(output.values[offset], quadrature.weights()[node], tolerance);
        QMOMKIT_VERIFY_NEAR(output.values[offset + 1], quadrature.abscissas()[0][node], tolerance);
        QMOMKIT_VERIFY_NEAR(output.values[offset + 2], quadrature.abscissas()[1][node], tolerance);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    // Output verification is intentionally workflow-level: compute public QMOM
    // and CQMOM quadratures, write moments/quadrature with OutputOptions, and
    // parse the files back as a downstream text consumer would.
    if (argc != 4)
    {
        std::cerr << "Usage: output_workflows_verification <qmom.qmom> <cqmom.qmom> <output-directory>\n";
        return 2;
    }

    const std::filesystem::path qmom_configuration = argv[1];
    const std::filesystem::path cqmom_configuration = argv[2];
    const std::filesystem::path output_directory = argv[3];

    std::filesystem::create_directories(output_directory);
    remove_known_outputs(output_directory);

    const auto configured_qmom = qmomkit::config::load_qmom(qmom_configuration);
    if (!configured_qmom.ok())
    {
        print_configuration_errors(configured_qmom);
        return 1;
    }
    const qmomkit::QmomResult qmom_result = qmomkit::compute_qmom(*configured_qmom.value);
    QMOMKIT_VERIFY_CHECK(qmom_result.status == qmomkit::QmomStatus::success);

    const auto configured_cqmom = qmomkit::config::load_cqmom(cqmom_configuration);
    if (!configured_cqmom.ok())
    {
        print_configuration_errors(configured_cqmom);
        return 1;
    }
    const qmomkit::CqmomResult cqmom_result = qmomkit::compute_cqmom(*configured_cqmom.value);
    QMOMKIT_VERIFY_CHECK(cqmom_result.status == qmomkit::CqmomStatus::success);

    const qmomkit::OutputOptions disabled_options{output_directory, false, false, std::nullopt};
    write_qmom_outputs(disabled_options, "disabled", *configured_qmom.value, qmom_result);
    QMOMKIT_VERIFY_CHECK(!std::filesystem::exists(output_directory / "disabled_moments.txt"));
    QMOMKIT_VERIFY_CHECK(!std::filesystem::exists(output_directory / "disabled_quadrature.txt"));

    const qmomkit::OutputOptions enabled_options{output_directory, true, true, 0.5};
    write_qmom_outputs(enabled_options, "qmom", *configured_qmom.value, qmom_result);
    write_cqmom_outputs(enabled_options, "cqmom", cqmom_result);

    verify_qmom_moment_output(output_directory / "qmom_moments.txt", configured_qmom.value->moments);
    verify_qmom_quadrature_output(output_directory / "qmom_quadrature.txt", qmom_result);
    verify_cqmom_moment_output(output_directory / "cqmom_moments.txt");
    verify_cqmom_quadrature_output(output_directory / "cqmom_quadrature.txt", cqmom_result.quadrature);

    std::cout << "QMOMKit output writer workflow verification\n";
    std::cout << "QMOM configuration: " << qmom_configuration << '\n';
    std::cout << "CQMOM configuration: " << cqmom_configuration << '\n';
    std::cout << "Output directory: " << output_directory << '\n';
    std::cout << "Tolerance: " << std::scientific << tolerance << '\n';
    std::cout << "Verified moment and quadrature text output for public QMOM and CQMOM workflows.\n";
    return 0;
}
