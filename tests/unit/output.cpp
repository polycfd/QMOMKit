#include "qmomkit/output.hpp"

#include "internal/cqmom_case.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "test_case_dispatch.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

[[nodiscard]] std::filesystem::path temporary_output_path(const std::string& stem)
{
    std::filesystem::path path = std::filesystem::temp_directory_path() / (stem + ".txt");
    std::filesystem::remove(path);
    return path;
}

[[nodiscard]] std::string read_text_file(const std::filesystem::path& path)
{
    std::ifstream input(path);
    assert(input);
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

[[nodiscard]] std::string format_double(const double value)
{
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
}

[[nodiscard]] std::string format_scientific(const double value)
{
    std::ostringstream output;
    output << std::scientific << std::setprecision(17) << value;
    return output.str();
}

[[nodiscard]] std::string normalize_line(const std::string& line)
{
    std::istringstream input(line);
    std::ostringstream output;

    std::string token;
    bool first = true;

    while (input >> token)
    {
        if (!first)
            output << ' ';

        output << token;
        first = false;
    }
    return output.str();
}

void assert_text_ignore_whitespaces(const std::string& actual, const std::string& expected)
{
    std::istringstream actual_stream(actual);
    std::istringstream expected_stream(expected);

    std::string actual_line;
    std::string expected_line;

    while (std::getline(actual_stream, actual_line) && std::getline(expected_stream, expected_line))
        assert(normalize_line(actual_line) == normalize_line(expected_line));

    assert(!std::getline(actual_stream, actual_line));
    assert(!std::getline(expected_stream, expected_line));
}

[[nodiscard]] std::string cqmom_label(const qmomkit::MultiIndex& index)
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
    return label.str();
}

template <class Function> void assert_throws_invalid_argument(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
}

template <class Function> void assert_throws_invalid_argument_message(Function&& function, const std::string& expected_message)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument& error)
    {
        threw = true;
        assert(error.what() == expected_message);
    }
    assert(threw);
}

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

[[nodiscard]] std::vector<double> raw_moments_from_grid(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<std::size_t>& extents,
                                                        const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
        {
            const std::size_t flat_node = hierarchy.flatten_node_index(node_index);
            double contribution = weights[flat_node];
            for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][flat_node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

[[nodiscard]] std::vector<double> raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                                             const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

void test_qmom_all_moments()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_qmom_all_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_qmom(path, 4);
        writer.write(0.0, qmomkit::MomentSet({1.0, 1.0 / 3.0, 2.5, -4.0}));
        writer.write(0.1, qmomkit::MomentSet({2.0, 4.0, 6.0, 8.0}));
    }

    const std::string expected = "# time M[0] M[1] M[2] M[3]\n" + format_scientific(0.0) + " " + format_scientific(1.0) + " " +
                                 format_scientific(1.0 / 3.0) + " " + format_scientific(2.5) + " " + format_scientific(-4.0) + "\n" +
                                 format_scientific(0.1) + " " + format_scientific(2.0) + " " + format_scientific(4.0) + " " + format_scientific(6.0) +
                                 " " + format_scientific(8.0) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_qmom_selected_moments()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_qmom_selected_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_qmom(path, std::vector<std::size_t>{3, 1});
        writer.write(0.25, qmomkit::MomentSet({10.0, 20.0, 30.0, 40.0}));
    }

    const std::string expected =
        "# time M[3] M[1]\n" + format_scientific(0.25) + " " + format_scientific(40.0) + " " + format_scientific(20.0) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_cqmom_all_moments()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_all_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, {2, 3});
        writer.write(0.0, qmomkit::RawMomentTensor{{2, 3}, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}});
    }

    const std::string expected = "# time M[0,0] M[0,1] M[0,2] M[1,0] M[1,1] M[1,2]\n" + format_scientific(0.0) + " " + format_scientific(1.0) + " " +
                                 format_scientific(2.0) + " " + format_scientific(3.0) + " " + format_scientific(4.0) + " " + format_scientific(5.0) +
                                 " " + format_scientific(6.0) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_cqmom_selected_moments()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_selected_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, {2, 3}, {{1, 2}, {0, 1}});
        writer.write(0.5, qmomkit::RawMomentTensor{{2, 3}, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}});
    }

    const std::string expected =
        "# time M[1,2] M[0,1]\n" + format_scientific(0.5) + " " + format_scientific(6.0) + " " + format_scientific(2.0) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_cqmom_moment_data_default_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (std::size_t index = 0; index < layout.required_indices.size(); ++index)
    {
        values.push_back(10.0 + static_cast<double>(index));
    }

    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_moment_data_default_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, layout);
        writer.write(0.75, qmomkit::CqmomMomentData{layout, values});
    }

    std::string expected = "# time";
    for (const qmomkit::MultiIndex& moment_index : layout.required_indices)
    {
        expected += " " + cqmom_label(moment_index);
    }
    expected += "\n" + format_scientific(0.75);
    for (const double value : values)
    {
        expected += " " + format_scientific(value);
    }
    expected += "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_requested_cqmom_moments_are_reconstructed_from_quadrature()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, {}});
    assert(result.ok());

    const qmomkit::RawMomentRequest request{{{1, 0}, {0, 2}, {1, 1}}};
    const qmomkit::SelectedRawMomentResult reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, request);
    const std::filesystem::path path = temporary_output_path("qmomkit_requested_cqmom_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, request);
        writer.write(0.875, result.quadrature);
    }

    const std::string expected = "# time M[1,0] M[0,2] M[1,1]\n" + format_scientific(0.875) + " " + format_scientific(reconstructed.values[0]) + " " +
                                 format_scientific(reconstructed.values[1]) + " " + format_scientific(reconstructed.values[2]) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_requested_cqmom_moments_preserve_duplicate_columns()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, {}});
    assert(result.ok());

    const qmomkit::RawMomentRequest request{{{1, 0}, {1, 0}, {0, 0}}};
    const qmomkit::SelectedRawMomentResult reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, request);
    const std::filesystem::path path = temporary_output_path("qmomkit_requested_cqmom_duplicate_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, request);
        writer.write(1.0, result.quadrature);
    }

    const std::string expected = "# time M[1,0] M[1,0] M[0,0]\n" + format_scientific(1.0) + " " + format_scientific(reconstructed.values[0]) + " " +
                                 format_scientific(reconstructed.values[1]) + " " + format_scientific(reconstructed.values[2]) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_cqmom_moment_data_writer_rejects_layout_mismatch()
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    const qmomkit::CqmomMomentLayout wrong_layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 1}));

    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_moment_data_layout_mismatch");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, layout);
        assert_throws_invalid_argument_message(
            [&writer, &wrong_layout]
            { writer.write(0.0, qmomkit::CqmomMomentData{wrong_layout, std::vector<double>(wrong_layout.required_indices.size(), 1.0)}); },
            "CQMOM moment data output layout does not match the writer layout.");
    }
    std::filesystem::remove(path);
}

void test_requested_cqmom_writer_rejects_rank_mismatch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, {}});
    assert(result.ok());

    const std::filesystem::path path = temporary_output_path("qmomkit_requested_cqmom_rank_mismatch");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(path, qmomkit::RawMomentRequest{{{1}}});
        assert_throws_invalid_argument([&writer, &result] { writer.write(0.0, result.quadrature); });
    }
    std::filesystem::remove(path);
}

void test_writer_rejects_mismatched_inputs()
{
    const std::filesystem::path qmom_path = temporary_output_path("qmomkit_qmom_invalid_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_qmom(qmom_path, 2);
        assert_throws_invalid_argument([&writer] { writer.write(0.0, qmomkit::MomentSet({1.0})); });
    }
    std::filesystem::remove(qmom_path);

    const std::filesystem::path cqmom_path = temporary_output_path("qmomkit_cqmom_invalid_moments");
    {
        qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(cqmom_path, {2, 2});
        assert_throws_invalid_argument([&writer] { writer.write(0.0, qmomkit::RawMomentTensor{{2, 3}, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}}); });
    }
    std::filesystem::remove(cqmom_path);
}

void test_qmom_quadrature_output()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_qmom_quadrature");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_qmom(path, 2);
        writer.write(0.0, qmomkit::QmomResult{{-1.0, 2.0}, {0.25, 0.75}});
    }

    const std::string expected = "# time node_index weight x0\n" + format_scientific(0.0) + " " + format_scientific(0.0) + " " +
                                 format_scientific(0.25) + " " + format_scientific(-1.0) + "\n" + format_scientific(0.0) + " " +
                                 format_scientific(1.0) + " " + format_scientific(0.75) + " " + format_scientific(2.0) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_quadrature_1d_output()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_quadrature_1d");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_qmom(path, 2);
        writer.write(0.125, qmomkit::Quadrature1D{{-0.75, 1.25}, {0.4, 0.6}});
    }

    const std::string expected = "# time node_index weight x0\n" + format_scientific(0.125) + " " + format_scientific(0.0) + " " +
                                 format_scientific(0.4) + " " + format_scientific(-0.75) + "\n" + format_scientific(0.125) + " " +
                                 format_scientific(1.0) + " " + format_scientific(0.6) + " " + format_scientific(1.25) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_qmom_adaptive_quadrature_output_accepts_changing_node_counts()
{
    const std::filesystem::path path = temporary_output_path("qmomkit_qmom_adaptive_quadrature_history");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_qmom(path, 2);
        writer.write(0.0, qmomkit::QmomResult{{-1.0, 1.0}, {0.5, 0.5}});
        writer.write(0.05, qmomkit::QmomResult{{0.0}, {1.0}});
        writer.write(0.15, qmomkit::QmomResult{{-0.8944271909999159, 0.8944271909999159}, {0.5, 0.5}});
        writer.write(0.25, qmomkit::QmomResult{{0.0}, {1.0}});
    }

    const std::string expected =
        "# time node_index weight x0\n" + format_scientific(0.0) + " " + format_scientific(0.0) + " " + format_scientific(0.5) + " " +
        format_scientific(-1.0) + "\n" + format_scientific(0.0) + " " + format_scientific(1.0) + " " + format_scientific(0.5) + " " +
        format_scientific(1.0) + "\n" + format_scientific(0.05) + " " + format_scientific(0.0) + " " + format_scientific(1.0) + " " +
        format_scientific(0.0) + "\n" + format_scientific(0.15) + " " + format_scientific(0.0) + " " + format_scientific(0.5) + " " +
        format_scientific(-0.8944271909999159) + "\n" + format_scientific(0.15) + " " + format_scientific(1.0) + " " + format_scientific(0.5) + " " +
        format_scientific(0.8944271909999159) + "\n" + format_scientific(0.25) + " " + format_scientific(0.0) + " " + format_scientific(1.0) + " " +
        format_scientific(0.0) + "\n";
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_cqmom_quadrature_output()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.28, 0.07, 0.52, 0.13};
    const std::vector<std::vector<double>> abscissas{{-0.75, -0.75, 1.25, 1.25}, {0.2, 1.4, 0.2, 1.4}};
    const std::vector<double> moments = raw_moments_from_grid(hierarchy, extents, weights, abscissas);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, {}});

    assert(result.ok());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
    {
        const std::size_t flat_node = hierarchy.flatten_node_index(node_index);
        assert_near(result.quadrature.weights()[flat_node], weights[flat_node], 1.0e-12);
        assert_near(result.quadrature.abscissas()[0][flat_node], abscissas[0][flat_node], 1.0e-12);
        assert_near(result.quadrature.abscissas()[1][flat_node], abscissas[1][flat_node], 1.0e-12);
    }

    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_quadrature");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_cqmom(path, result.quadrature.number_of_dimensions());
        writer.write(0.5, result.quadrature);
    }

    std::string expected = "# time node_index path_0 path_1 weight x0 x1\n";
    for (std::size_t node = 0; node < result.quadrature.total_nodes(); ++node)
    {
        expected += format_scientific(0.5) + " " + format_scientific(static_cast<double>(node));
        for (std::size_t level = 0; level < result.quadrature.number_of_dimensions(); ++level)
        {
            expected += " " + format_scientific(static_cast<double>(result.quadrature.terminal_paths()[node].node_index_at_level(level)));
        }
        expected += " " + format_scientific(result.quadrature.weights()[node]);
        for (std::size_t dimension = 0; dimension < result.quadrature.number_of_dimensions(); ++dimension)
        {
            expected += " " + format_scientific(result.quadrature.abscissas()[dimension][node]);
        }
        expected += "\n";
    }
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_cqmom_adaptive_quadrature_output_accepts_changing_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};

    const std::vector<double> fixed_weights{0.28, 0.07, 0.52, 0.13};
    const std::vector<std::vector<double>> fixed_abscissas{{-0.75, -0.75, 1.25, 1.25}, {0.2, 1.4, 0.2, 1.4}};
    const std::vector<double> fixed_moments = raw_moments_from_grid(hierarchy, extents, fixed_weights, fixed_abscissas);
    const qmomkit::CqmomResult fixed = qmomkit::compute_cqmom({hierarchy, {extents, fixed_moments}, {}, {}});
    assert(fixed.ok());

    const std::vector<double> adaptive_weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> adaptive_abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    const std::vector<double> adaptive_moments = raw_moments_from_terminals(extents, adaptive_weights, adaptive_abscissas);
    const qmomkit::QmomOptions adaptive_qmom{true, 1.0e-10, 1.0e-10, 1.0e-12};
    const qmomkit::CqmomResult adaptive = qmomkit::compute_cqmom({hierarchy, {extents, adaptive_moments}, {}, adaptive_qmom, true});
    assert(adaptive.ok());
    assert(adaptive.quadrature.total_nodes() == 3);

    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_adaptive_quadrature_history");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_cqmom(path, 2);
        writer.write(0.0, fixed.quadrature);
        writer.write(0.5, adaptive.quadrature);
    }

    std::string expected = "# time node_index path_0 path_1 weight x0 x1\n";
    for (std::size_t node = 0; node < fixed.quadrature.total_nodes(); ++node)
    {
        expected += format_scientific(0.0) + " " + format_scientific(static_cast<double>(node)) + " " +
                    format_scientific(static_cast<double>(fixed.quadrature.terminal_paths()[node].node_index_at_level(0))) + " " +
                    format_scientific(static_cast<double>(fixed.quadrature.terminal_paths()[node].node_index_at_level(1))) + " " +
                    format_scientific(fixed.quadrature.weights()[node]) + " " + format_scientific(fixed.quadrature.abscissas()[0][node]) + " " +
                    format_scientific(fixed.quadrature.abscissas()[1][node]) + "\n";
    }
    for (std::size_t node = 0; node < adaptive.quadrature.total_nodes(); ++node)
    {
        expected += format_scientific(0.5) + " " + format_scientific(static_cast<double>(node)) + " " +
                    format_scientific(static_cast<double>(adaptive.quadrature.terminal_paths()[node].node_index_at_level(0))) + " " +
                    format_scientific(static_cast<double>(adaptive.quadrature.terminal_paths()[node].node_index_at_level(1))) + " " +
                    format_scientific(adaptive.quadrature.weights()[node]) + " " + format_scientific(adaptive.quadrature.abscissas()[0][node]) + " " +
                    format_scientific(adaptive.quadrature.abscissas()[1][node]) + "\n";
    }
    assert_text_ignore_whitespaces(read_text_file(path), expected);
    std::filesystem::remove(path);
}

void test_quadrature_writer_rejects_mismatched_inputs()
{
    const std::filesystem::path qmom_path = temporary_output_path("qmomkit_qmom_invalid_quadrature");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_qmom(qmom_path, 2);
        assert_throws_invalid_argument_message([&writer] { writer.write(0.0, qmomkit::QmomResult{{}, {}}); },
                                               "QMOM quadrature output requires at least one active node.");
        assert_throws_invalid_argument_message([&writer] { writer.write(0.0, qmomkit::QmomResult{{1.0}, {0.25, 0.75}}); },
                                               "QMOM quadrature output weight and abscissa counts must match.");
        assert_throws_invalid_argument_message([&writer] { writer.write(0.0, qmomkit::QmomResult{{-1.0, 0.0, 1.0}, {0.25, 0.5, 0.25}}); },
                                               "QMOM quadrature output active node count exceeds the writer maximum.");
    }
    std::filesystem::remove(qmom_path);

    const std::filesystem::path cqmom_path = temporary_output_path("qmomkit_cqmom_invalid_quadrature");
    {
        const qmomkit::PhaseSpaceHierarchy hierarchy({2});
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {{4}, {1.0, 0.0, 1.0, 0.0}}, {}, {}});
        assert(result.ok());

        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_cqmom(cqmom_path, 2);
        assert_throws_invalid_argument([&writer, &result] { writer.write(0.0, result.quadrature); });
    }
    std::filesystem::remove(cqmom_path);
}

void test_cqmom_quadrature_writer_rejects_malformed_terminal_storage()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1});
    const std::filesystem::path path = temporary_output_path("qmomkit_cqmom_malformed_terminal_quadrature");
    {
        qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_cqmom(path, 2);

        const qmomkit::CqmomQuadrature empty =
            qmomkit::internal::CqmomConstructionAccess::create_quadrature(hierarchy, {}, {{}, {}}, true, false, {});
        assert_throws_invalid_argument_message([&writer, &empty] { writer.write(0.0, empty); },
                                               "CQMOM quadrature output requires at least one active terminal node.");

        const qmomkit::CqmomQuadrature missing_path =
            qmomkit::internal::CqmomConstructionAccess::create_quadrature(hierarchy, {1.0}, {{-1.0}, {2.0}}, true, false, {});
        assert_throws_invalid_argument_message([&writer, &missing_path] { writer.write(0.0, missing_path); },
                                               "CQMOM quadrature output requires one terminal path per active terminal node.");

        const qmomkit::CqmomQuadrature wrong_depth = qmomkit::internal::CqmomConstructionAccess::create_quadrature(
            hierarchy, {1.0}, {{-1.0}, {2.0}}, true, false, {qmomkit::ConditioningPath({0})});
        assert_throws_invalid_argument_message([&writer, &wrong_depth] { writer.write(0.0, wrong_depth); },
                                               "CQMOM quadrature output terminal path depth does not match the writer dimension count.");

        const qmomkit::CqmomQuadrature wrong_abscissa_count = qmomkit::internal::CqmomConstructionAccess::create_quadrature(
            hierarchy, {1.0}, {{-1.0, 0.0}, {2.0}}, true, false, {qmomkit::ConditioningPath({0, 0})});
        assert_throws_invalid_argument_message([&writer, &wrong_abscissa_count] { writer.write(0.0, wrong_abscissa_count); },
                                               "CQMOM quadrature output abscissa node count is inconsistent.");

        const qmomkit::CqmomQuadrature nonfinite_weight = qmomkit::internal::CqmomConstructionAccess::create_quadrature(
            hierarchy, {std::numeric_limits<double>::infinity()}, {{-1.0}, {2.0}}, true, false, {qmomkit::ConditioningPath({0, 0})});
        assert_throws_invalid_argument_message([&writer, &nonfinite_weight] { writer.write(0.0, nonfinite_weight); },
                                               "CQMOM quadrature output weights must be finite.");

        const qmomkit::CqmomQuadrature nonfinite_abscissa = qmomkit::internal::CqmomConstructionAccess::create_quadrature(
            hierarchy, {1.0}, {{std::numeric_limits<double>::quiet_NaN()}, {2.0}}, true, false, {qmomkit::ConditioningPath({0, 0})});
        assert_throws_invalid_argument_message([&writer, &nonfinite_abscissa] { writer.write(0.0, nonfinite_abscissa); },
                                               "CQMOM quadrature output abscissas must be finite.");
    }
    std::filesystem::remove(path);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"qmom_all_moments", test_qmom_all_moments},
            {"qmom_selected_moments", test_qmom_selected_moments},
            {"cqmom_all_moments", test_cqmom_all_moments},
            {"cqmom_selected_moments", test_cqmom_selected_moments},
            {"cqmom_moment_data_default_moments", test_cqmom_moment_data_default_moments},
            {"requested_cqmom_moments_are_reconstructed_from_quadrature", test_requested_cqmom_moments_are_reconstructed_from_quadrature},
            {"requested_cqmom_moments_preserve_duplicate_columns", test_requested_cqmom_moments_preserve_duplicate_columns},
            {"cqmom_moment_data_writer_rejects_layout_mismatch", test_cqmom_moment_data_writer_rejects_layout_mismatch},
            {"requested_cqmom_writer_rejects_rank_mismatch", test_requested_cqmom_writer_rejects_rank_mismatch},
            {"writer_rejects_mismatched_inputs", test_writer_rejects_mismatched_inputs},
            {"qmom_quadrature_output", test_qmom_quadrature_output},
            {"quadrature_1d_output", test_quadrature_1d_output},
            {"qmom_adaptive_quadrature_output_accepts_changing_node_counts", test_qmom_adaptive_quadrature_output_accepts_changing_node_counts},
            {"cqmom_quadrature_output", test_cqmom_quadrature_output},
            {"cqmom_adaptive_quadrature_output_accepts_changing_node_counts", test_cqmom_adaptive_quadrature_output_accepts_changing_node_counts},
            {"quadrature_writer_rejects_mismatched_inputs", test_quadrature_writer_rejects_mismatched_inputs},
            {"cqmom_quadrature_writer_rejects_malformed_terminal_storage", test_cqmom_quadrature_writer_rejects_malformed_terminal_storage},
        });
}
