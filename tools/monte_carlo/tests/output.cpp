#include "monte_carlo/output.hpp"

#include "monte_carlo/moments.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"
#include "qmomkit/raw_moment_tensor.hpp"
#include "test_case_dispatch.hpp"
#include "test_checks.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Coverage for the Monte Carlo output module. The writers produce the artifacts
// a validation run is judged from, and `docs/developer/output_design.md` fixes
// their format: whitespace-delimited columns, one commented data header,
// optional metadata comments ahead of it, bracketed multi-index labels such as
// `M[0]` and `M[0,0]`, and 17-digit precision so a written value round-trips
// exactly.
//
// The moment-history factories and row writers are overload sets whose data
// argument must match the schema chosen at construction. These cases exercise
// each pairing, because a mismatched pairing is a compile-time contract that is
// only meaningful if the intended combinations are known to work.

namespace
{

std::filesystem::path scratch_path(const std::string& name)
{
    const std::filesystem::path path = std::filesystem::temp_directory_path() / ("qmomkit_mc_output_" + name + ".txt");
    std::filesystem::remove(path);
    return path;
}

std::vector<std::string> read_lines(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::string> split(const std::string& line)
{
    std::istringstream tokens(line);
    std::vector<std::string> fields;
    std::string token;
    while (tokens >> token)
    {
        fields.push_back(token);
    }
    return fields;
}

// Comment lines are right-aligned into their columns, so the marker is the
// first non-whitespace character rather than the first character.
bool is_comment(const std::string& line)
{
    for (const char character : line)
    {
        if (character != ' ' && character != '\t')
        {
            return character == '#';
        }
    }
    return false;
}

// The data header is the last comment line; anything before it is metadata.
std::size_t data_header_index(const std::vector<std::string>& lines)
{
    std::size_t header = 0;
    for (std::size_t index = 0; index < lines.size(); ++index)
    {
        if (is_comment(lines[index]))
        {
            header = index;
        }
    }
    return header;
}

bool contains_field(const std::vector<std::string>& fields, const std::string& wanted)
{
    for (const std::string& field : fields)
    {
        if (field == wanted)
        {
            return true;
        }
    }
    return false;
}

qmomkit::monte_carlo::Ensemble two_sample_ensemble()
{
    return qmomkit::monte_carlo::Ensemble{{1.0}, {3.0}};
}

void test_one_dimensional_moment_text_records_metadata_and_columns()
{
    const qmomkit::monte_carlo::Ensemble ensemble = two_sample_ensemble();
    const double number_density = 6.0;
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate = qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, 2);
    CHECK(estimate.ok());

    std::ostringstream output;
    qmomkit::monte_carlo::write_moment_text(output, estimate, ensemble.size(), 4242, number_density);
    const std::string text = output.str();
    CHECK(!text.empty());

    std::istringstream reader(text);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(reader, line))
    {
        lines.push_back(line);
    }
    CHECK(lines.size() >= 2);

    // The sample count, seed, and absolute number density are what make a
    // stochastic estimate reproducible, so each must survive into the file as
    // its own metadata comment rather than merely appearing somewhere.
    CHECK(contains_field(split(lines[0]), "samples=2"));
    CHECK(contains_field(split(lines[1]), "seed=4242"));
    CHECK(contains_field(split(lines[2]), "number_density=6"));

    const std::size_t header = data_header_index(lines);
    CHECK(is_comment(lines[header]));
    for (std::size_t index = 0; index < header; ++index)
    {
        CHECK(is_comment(lines[index]));
    }
    CHECK(lines.size() == header + 1 + estimate.moments.size());

    // Data rows are `alpha_0 value standard_error`, and m_0 is the number
    // density for an equally weighted ensemble.
    const std::vector<std::string> first_row = split(lines[header + 1]);
    CHECK(first_row.size() == 3);
    CHECK_CLOSE(std::stod(first_row[1]), number_density, 1.0e-12);
}

void test_tensor_moment_text_uses_phase_space_multi_index_columns()
{
    const qmomkit::monte_carlo::Ensemble ensemble{{1.0, 4.0}, {3.0, 8.0}};
    const double number_density = 2.0;
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, std::vector<std::size_t>{1, 2});
    CHECK(estimate.ok());

    std::ostringstream output;
    qmomkit::monte_carlo::write_moment_text(output, estimate, ensemble.size(), 7, number_density);

    std::istringstream reader(output.str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(reader, line))
    {
        lines.push_back(line);
    }
    const std::size_t header = data_header_index(lines);

    // Extents {2, 3} give six flat entries, each carrying its own multi-index.
    CHECK(estimate.moments.values.size() == 6);
    CHECK(lines.size() == header + 1 + estimate.moments.values.size());

    const std::vector<std::string> first_row = split(lines[header + 1]);
    CHECK(first_row.size() == 4);
    CHECK(first_row[0] == "0");
    CHECK(first_row[1] == "0");
    CHECK_CLOSE(std::stod(first_row[2]), number_density, 1.0e-12);
}

void test_one_dimensional_history_writer_round_trips_values()
{
    const std::filesystem::path path = scratch_path("history_1d");
    const qmomkit::MomentSet moments({2.0, 5.0, 17.0});
    const qmomkit::monte_carlo::Ensemble ensemble = two_sample_ensemble();
    const qmomkit::OutputOptions output_options;

    {
        qmomkit::MomentWriter writer = qmomkit::monte_carlo::make_moment_history_writer(path, moments.size(), output_options);
        qmomkit::monte_carlo::write_moment_history_row(writer, 0.25, moments, ensemble, 2.0, output_options);
        qmomkit::monte_carlo::write_moment_history_row(writer, 0.50, moments, ensemble, 2.0, output_options);
    }

    const std::vector<std::string> lines = read_lines(path);
    const std::size_t header = data_header_index(lines);
    CHECK(is_comment(lines[header]));

    const std::vector<std::string> header_fields = split(lines[header]);
    CHECK(contains_field(header_fields, "time"));
    CHECK(contains_field(header_fields, "M[0]"));
    CHECK(contains_field(header_fields, "M[2]"));

    CHECK(lines.size() == header + 3);
    const std::vector<std::string> row = split(lines[header + 1]);
    CHECK(row.size() == header_fields.size() - 1);

    // 17-digit precision is documented so a written value round-trips exactly.
    const std::vector<std::string> second = split(lines[header + 2]);
    CHECK(std::stod(second[0]) == 0.50);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        CHECK(std::stod(row[row.size() - moments.size() + order]) == moments[order]);
    }
    std::filesystem::remove(path);
}

void test_selected_tensor_history_writer_emits_requested_columns_only()
{
    const std::filesystem::path path = scratch_path("history_selected");
    const std::vector<std::size_t> extents{2, 2};
    const qmomkit::RawMomentTensor moments{extents, {1.0, 2.0, 3.0, 4.0}};
    const std::vector<qmomkit::MultiIndex> selected{{0, 0}, {1, 1}};
    const qmomkit::monte_carlo::Ensemble ensemble{{1.0, 1.0}, {2.0, 2.0}};
    const qmomkit::OutputOptions output_options;

    {
        qmomkit::MomentWriter writer = qmomkit::monte_carlo::make_moment_history_writer(path, extents, selected, output_options);
        qmomkit::monte_carlo::write_moment_history_row(writer, 0.0, moments, ensemble, 1.0, output_options);
    }

    const std::vector<std::string> lines = read_lines(path);
    const std::size_t header = data_header_index(lines);
    const std::vector<std::string> header_fields = split(lines[header]);

    // Asking for two of the four tensor entries must not widen the table to the
    // full rectangle.
    CHECK(contains_field(header_fields, "M[0,0]"));
    CHECK(contains_field(header_fields, "M[1,1]"));
    CHECK(!contains_field(header_fields, "M[0,1]"));
    CHECK(!contains_field(header_fields, "M[1,0]"));

    const std::vector<std::string> row = split(lines[header + 1]);
    CHECK(row.size() == selected.size() + 1);
    CHECK(std::stod(row[1]) == 1.0);
    CHECK(std::stod(row[2]) == 4.0);
    std::filesystem::remove(path);
}

void test_lean_cqmom_history_writer_follows_required_layout()
{
    const std::filesystem::path path = scratch_path("history_lean");
    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments{layout, std::vector<double>(layout.required_indices.size(), 1.5)};
    const qmomkit::monte_carlo::Ensemble ensemble{{1.0, 1.0}, {2.0, 2.0}};
    const qmomkit::OutputOptions output_options;

    {
        qmomkit::MomentWriter writer = qmomkit::monte_carlo::make_moment_history_writer(path, layout, output_options);
        qmomkit::monte_carlo::write_moment_history_row(writer, 1.0, moments, ensemble, 1.5, output_options);
    }

    const std::vector<std::string> lines = read_lines(path);
    const std::size_t header = data_header_index(lines);
    const std::vector<std::string> header_fields = split(lines[header]);
    const std::vector<std::string> row = split(lines[header + 1]);

    // One column per required moment, plus time. The layout decides the columns,
    // so the writer must not fall back to a dense rectangle.
    CHECK(header_fields.size() == layout.required_indices.size() + 2);
    CHECK(row.size() == layout.required_indices.size() + 1);
    CHECK(std::stod(row[0]) == 1.0);
    for (std::size_t index = 0; index < layout.required_indices.size(); ++index)
    {
        CHECK(std::stod(row[index + 1]) == 1.5);
    }
    std::filesystem::remove(path);
}

void test_generalized_columns_are_appended_when_requested()
{
    const std::filesystem::path bare = scratch_path("history_bare");
    const std::filesystem::path generalized = scratch_path("history_generalized");
    const qmomkit::MomentSet moments({2.0, 4.0});
    const qmomkit::monte_carlo::Ensemble ensemble = two_sample_ensemble();

    qmomkit::OutputOptions plain_options;
    qmomkit::OutputOptions generalized_options;
    generalized_options.generalized_moment_orders = {{0.5}, {1.5}};

    {
        qmomkit::MomentWriter writer = qmomkit::monte_carlo::make_moment_history_writer(bare, moments.size(), plain_options);
        qmomkit::monte_carlo::write_moment_history_row(writer, 0.0, moments, ensemble, 2.0, plain_options);
    }
    {
        qmomkit::MomentWriter writer = qmomkit::monte_carlo::make_moment_history_writer(generalized, moments.size(), generalized_options);
        qmomkit::monte_carlo::write_moment_history_row(writer, 0.0, moments, ensemble, 2.0, generalized_options);
    }

    const std::vector<std::string> plain_lines = read_lines(bare);
    const std::vector<std::string> generalized_lines = read_lines(generalized);
    const std::vector<std::string> plain_row = split(plain_lines[data_header_index(plain_lines) + 1]);
    const std::vector<std::string> generalized_row = split(generalized_lines[data_header_index(generalized_lines) + 1]);

    // Requesting generalized orders adds columns; it must not disturb the
    // standard ones that precede them.
    CHECK(generalized_row.size() == plain_row.size() + generalized_options.generalized_moment_orders.size());
    for (std::size_t field = 0; field < plain_row.size(); ++field)
    {
        CHECK(generalized_row[field] == plain_row[field]);
    }

    // Equal weights of 1.0 each: n0 / N * sum xi^0.5 = 2 / 2 * (1 + sqrt(3)).
    CHECK_CLOSE(std::stod(generalized_row[plain_row.size()]), 1.0 + std::sqrt(3.0), 1.0e-12);
    std::filesystem::remove(bare);
    std::filesystem::remove(generalized);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::monte_carlo::test::dispatch_test_case(
        argc, argv,
        {
            {"one_dimensional_moment_text_records_metadata_and_columns", test_one_dimensional_moment_text_records_metadata_and_columns},
            {"tensor_moment_text_uses_phase_space_multi_index_columns", test_tensor_moment_text_uses_phase_space_multi_index_columns},
            {"one_dimensional_history_writer_round_trips_values", test_one_dimensional_history_writer_round_trips_values},
            {"selected_tensor_history_writer_emits_requested_columns_only", test_selected_tensor_history_writer_emits_requested_columns_only},
            {"lean_cqmom_history_writer_follows_required_layout", test_lean_cqmom_history_writer_follows_required_layout},
            {"generalized_columns_are_appended_when_requested", test_generalized_columns_are_appended_when_requested},
        });
}
