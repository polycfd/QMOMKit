#include "qmomkit/configuration.hpp"

#include "internal/checked_size.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace qmomkit::config
{
namespace
{

struct Entry
{
    std::string value;
    std::size_t line = 0;
};

struct ParsedFile
{
    std::map<std::string, Entry> scalars;
    std::map<MultiIndex, Entry> moments;
    std::map<std::string, Entry> sections;
    std::map<std::string, std::vector<std::size_t>> section_lines;
    std::vector<ConfigurationError> errors;
    // Indexed sections such as "initial.mixture.component.0" that an
    // initialization parser read. Any other section or key under an indexed
    // prefix is unknown.
    std::set<std::string> consumed_indexed_sections;
};

constexpr std::string_view mixture_component_section_prefix = "initial.mixture.component.";
constexpr std::string_view independent_product_dimension_section_prefix = "initial.dimension.";
constexpr std::string_view realizability_section = "time_integration.realizability";
constexpr std::string_view realizability_key_prefix = "time_integration.realizability.";

// The realizability settings a parser accepts. QMOM reads the unnumbered
// "[time_integration.realizability]" policy and CQMOM the per-coordinate
// "[time_integration.realizability.N]" policies; each rejects the other's, which
// would otherwise be accepted and have no effect. Parsers that read neither,
// such as the splitting-scheme parser, accept both.
enum class RealizabilitySettings
{
    qmom,
    cqmom,
    any
};

bool accepts_unnumbered_realizability(const RealizabilitySettings settings)
{
    return settings != RealizabilitySettings::cqmom;
}

bool accepts_coordinate_realizability(const RealizabilitySettings settings)
{
    return settings != RealizabilitySettings::qmom;
}
constexpr std::string_view phase_space_evolution_scheme_key = "time_integration.scheme";
constexpr std::string_view node_advection_section = "time_integration.node_advection";
constexpr std::string_view node_advection_key_prefix = "time_integration.node_advection.";
constexpr std::string_view node_advection_scheme_key = "time_integration.node_advection.scheme";
constexpr std::string_view adaptive_section = "time_integration.node_advection.adaptive";
constexpr std::string_view adaptive_relative_tolerance_key = "time_integration.node_advection.adaptive.relative_tolerance";
constexpr std::string_view adaptive_absolute_tolerance_key = "time_integration.node_advection.adaptive.absolute_tolerance";
constexpr std::string_view adaptive_initial_step_key = "time_integration.node_advection.adaptive.initial_step";
constexpr std::string_view adaptive_max_step_key = "time_integration.node_advection.adaptive.max_step";
constexpr std::string_view adaptive_min_step_key = "time_integration.node_advection.adaptive.min_step";
constexpr std::string_view adaptive_max_steps_key = "time_integration.node_advection.adaptive.max_steps";
constexpr std::string_view output_section = "output";
constexpr std::string_view output_key_prefix = "output.";

std::string trim(const std::string& input)
{
    const std::size_t first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return {};
    }
    const std::size_t last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

void add_error(ParsedFile& parsed, const ConfigurationErrorCode code, const std::string& source_name, const std::size_t line,
               const std::string& message)
{
    parsed.errors.push_back({code, source_name, line, message});
}

bool parse_size(const std::string& text, std::size_t& value)
{
    if (text.empty() || text.front() == '-')
    {
        return false;
    }
    const char* const first = text.data();
    const char* const last = first + text.size();
    const auto result = std::from_chars(first, last, value);
    return result.ec == std::errc() && result.ptr == last;
}

// An index in a section or key name is accepted only in the canonical decimal
// form that std::to_string writes, "1" rather than "01": indexed sections are
// looked up by that form, so any other spelling would name a section that
// nothing reads.
bool parse_index(const std::string& text, std::size_t& index)
{
    return parse_size(text, index) && std::to_string(index) == text;
}

bool parse_coordinate_realizability_section(const std::string& section, std::size_t& coordinate)
{
    if (section.rfind(std::string(realizability_key_prefix), 0) != 0)
    {
        return false;
    }
    return parse_index(section.substr(realizability_key_prefix.size()), coordinate);
}

bool is_coordinate_realizability_key(const std::string& key)
{
    if (key.rfind(std::string(realizability_key_prefix), 0) != 0)
    {
        return false;
    }
    const std::size_t coordinate_end = key.find('.', realizability_key_prefix.size());
    if (coordinate_end == std::string::npos)
    {
        return false;
    }
    std::size_t coordinate = 0;
    return parse_index(key.substr(realizability_key_prefix.size(), coordinate_end - realizability_key_prefix.size()), coordinate);
}

bool is_configured_coordinate_realizability_key(const ParsedFile& parsed, const std::string& key)
{
    if (!is_coordinate_realizability_key(key))
    {
        return false;
    }
    const std::size_t coordinate_end = key.find('.', realizability_key_prefix.size());
    return parsed.sections.find(key.substr(0, coordinate_end)) != parsed.sections.end();
}

// Configuration numbers are plain decimal: an optional sign, digits with an
// optional fraction, and an optional exponent. Checking this form first keeps
// the accepted syntax the same on every standard library.
bool is_decimal_number(const std::string& text)
{
    const auto is_digit = [](const char character) { return std::isdigit(static_cast<unsigned char>(character)) != 0; };
    std::size_t position = 0;
    if (position < text.size() && (text[position] == '+' || text[position] == '-'))
    {
        ++position;
    }
    std::size_t mantissa_digits = 0;
    while (position < text.size() && is_digit(text[position]))
    {
        ++position;
        ++mantissa_digits;
    }
    if (position < text.size() && text[position] == '.')
    {
        ++position;
        while (position < text.size() && is_digit(text[position]))
        {
            ++position;
            ++mantissa_digits;
        }
    }
    if (mantissa_digits == 0)
    {
        return false;
    }
    if (position < text.size() && (text[position] == 'e' || text[position] == 'E'))
    {
        ++position;
        if (position < text.size() && (text[position] == '+' || text[position] == '-'))
        {
            ++position;
        }
        std::size_t exponent_digits = 0;
        while (position < text.size() && is_digit(text[position]))
        {
            ++position;
            ++exponent_digits;
        }
        if (exponent_digits == 0)
        {
            return false;
        }
    }
    return position == text.size();
}

// Parsed in the classic "C" locale, so the locale of a host application, which
// may use a decimal comma, cannot change how a configuration file is read.
// Values that overflow or fall below the normal range are rejected.
bool parse_double(const std::string& text, double& value)
{
    if (!is_decimal_number(text))
    {
        return false;
    }
    std::istringstream stream(text);
    stream.imbue(std::locale::classic());
    stream >> value;
    if (stream.fail() || stream.peek() != std::char_traits<char>::eof() || !std::isfinite(value))
    {
        return false;
    }
    return value == 0.0 || std::abs(value) >= std::numeric_limits<double>::min();
}

bool parse_size_list(const Entry& entry, std::vector<std::size_t>& values)
{
    std::istringstream stream(entry.value);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        std::size_t value = 0;
        if (!parse_size(trim(token), value))
        {
            return false;
        }
        values.push_back(value);
    }
    return !values.empty() && entry.value.back() != ',';
}

bool parse_double_list(const Entry& entry, std::vector<double>& values)
{
    std::istringstream stream(entry.value);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        double value = 0.0;
        if (!parse_double(trim(token), value))
        {
            return false;
        }
        values.push_back(value);
    }
    return !values.empty() && entry.value.back() != ',';
}

bool parse_optional_double_list(const Entry& entry, std::vector<std::optional<double>>& values)
{
    std::istringstream stream(entry.value);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        const std::string trimmed = trim(token);
        if (trimmed == "none")
        {
            values.push_back(std::nullopt);
            continue;
        }
        double value = 0.0;
        if (!parse_double(trimmed, value))
        {
            return false;
        }
        values.push_back(value);
    }
    return !values.empty() && entry.value.back() != ',';
}

bool validate_names(const ParsedFile& parsed, const std::size_t dimensions)
{
    const auto found = parsed.scalars.find("phase_space.names");
    if (found == parsed.scalars.end())
    {
        return true;
    }
    std::istringstream stream(found->second.value);
    std::string name;
    std::size_t count = 0;
    while (std::getline(stream, name, ','))
    {
        if (trim(name).empty())
        {
            return false;
        }
        ++count;
    }
    return count == dimensions && found->second.value.back() != ',';
}

bool parse_moment_key(const std::string& key, MultiIndex& index)
{
    if (key.size() < 4 || key.rfind("M[", 0) != 0 || key.back() != ']')
    {
        return false;
    }
    const std::string contents = key.substr(2, key.size() - 3);
    std::istringstream stream(contents);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        std::size_t component = 0;
        if (!parse_size(trim(token), component))
        {
            return false;
        }
        index.push_back(component);
    }
    return !index.empty() && contents.back() != ',';
}

std::string mixture_component_prefix(const std::string& initialization_prefix)
{
    return initialization_prefix + ".mixture.component.";
}

bool is_mixture_component_section(const std::string& section, const std::string& initialization_prefix)
{
    const std::string prefix = mixture_component_prefix(initialization_prefix);
    if (section.rfind(prefix, 0) != 0)
    {
        return false;
    }
    std::size_t component = 0;
    return parse_index(section.substr(prefix.size()), component);
}

bool is_independent_product_dimension_section(const std::string& section)
{
    if (section.rfind(std::string(independent_product_dimension_section_prefix), 0) != 0)
    {
        return false;
    }
    std::size_t dimension = 0;
    return parse_index(section.substr(independent_product_dimension_section_prefix.size()), dimension);
}

ParsedFile parse_file(const std::string_view text, const std::string& source_name)
{
    ParsedFile parsed;
    std::istringstream input{std::string(text)};
    std::string section;
    std::string line_text;
    std::size_t line_number = 0;
    while (std::getline(input, line_text))
    {
        ++line_number;
        const std::size_t comment = line_text.find_first_of("#;");
        if (comment != std::string::npos)
        {
            line_text.erase(comment);
        }
        line_text = trim(line_text);
        if (line_text.empty())
        {
            continue;
        }
        if (line_text.front() == '[')
        {
            if (line_text.back() != ']')
            {
                // A rejected header does not open a section: keys on subsequent
                // lines must not be silently attributed to whatever section was
                // open before this malformed line.
                section.clear();
                add_error(parsed, ConfigurationErrorCode::syntax_error, source_name, line_number, "Malformed section header.");
                continue;
            }
            const std::string candidate_section = trim(line_text.substr(1, line_text.size() - 2));
            if (candidate_section.empty())
            {
                section.clear();
                add_error(parsed, ConfigurationErrorCode::syntax_error, source_name, line_number, "Malformed section header.");
                continue;
            }
            section = candidate_section;
            parsed.sections.emplace(section, Entry{"", line_number});
            parsed.section_lines[section].push_back(line_number);
            continue;
        }

        const std::size_t equals = line_text.find('=');
        if (equals == std::string::npos)
        {
            add_error(parsed, ConfigurationErrorCode::syntax_error, source_name, line_number, "Expected key = value.");
            continue;
        }
        const std::string key = trim(line_text.substr(0, equals));
        const std::string value = trim(line_text.substr(equals + 1));
        if (key.empty() || value.empty())
        {
            add_error(parsed, ConfigurationErrorCode::syntax_error, source_name, line_number, "Configuration keys and values must not be empty.");
            continue;
        }

        MultiIndex moment_index;
        if (section == "moments" && parse_moment_key(key, moment_index))
        {
            if (!parsed.moments.emplace(moment_index, Entry{value, line_number}).second)
            {
                add_error(parsed, ConfigurationErrorCode::duplicate_entry, source_name, line_number, "Duplicate moment index '" + key + "'.");
            }
            continue;
        }
        if (section == "moments" && key.rfind("M[", 0) == 0)
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, line_number,
                      "Invalid moment index '" + key + "': indices must be comma-separated non-negative integers.");
            continue;
        }

        const std::string qualified_key = section.empty() ? key : section + "." + key;
        if (!parsed.scalars.emplace(qualified_key, Entry{value, line_number}).second)
        {
            add_error(parsed, ConfigurationErrorCode::duplicate_entry, source_name, line_number, "Duplicate key '" + qualified_key + "'.");
        }
    }
    return parsed;
}

const Entry* require_entry(ParsedFile& parsed, const std::string& key, const std::string& source_name)
{
    const auto found = parsed.scalars.find(key);
    if (found != parsed.scalars.end())
    {
        return &found->second;
    }
    add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0, "Missing required key '" + key + "'.");
    return nullptr;
}

const Entry* find_entry(const ParsedFile& parsed, const std::string& key)
{
    const auto found = parsed.scalars.find(key);
    if (found == parsed.scalars.end())
    {
        return nullptr;
    }
    return &found->second;
}

void reject_unknown_keys(ParsedFile& parsed, const std::vector<std::string>& allowed, const std::vector<std::string>& allowed_prefixes,
                         const RealizabilitySettings realizability, const std::string& source_name)
{
    for (const auto& item : parsed.scalars)
    {
        bool allowed_key = std::find(allowed.begin(), allowed.end(), item.first) != allowed.end();
        allowed_key =
            allowed_key || (accepts_coordinate_realizability(realizability) && is_configured_coordinate_realizability_key(parsed, item.first));
        for (const std::string& prefix : allowed_prefixes)
        {
            allowed_key = allowed_key || item.first.rfind(prefix, 0) == 0;
        }
        if (!allowed_key)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }
}

void reject_unknown_sections(ParsedFile& parsed, const std::vector<std::string>& allowed, const std::vector<std::string>& allowed_prefixes,
                             const RealizabilitySettings realizability, const std::string& source_name)
{
    for (const auto& item : parsed.sections)
    {
        bool allowed_section = std::find(allowed.begin(), allowed.end(), item.first) != allowed.end();
        std::size_t coordinate = 0;
        allowed_section =
            allowed_section || (accepts_coordinate_realizability(realizability) && parse_coordinate_realizability_section(item.first, coordinate));
        for (const std::string& prefix : allowed_prefixes)
        {
            allowed_section = allowed_section || item.first.rfind(prefix, 0) == 0;
        }
        if (!allowed_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
}

void validate_format_version(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* version = find_entry(parsed, "format_version");
    if (version != nullptr && version->value != "1")
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, version->line, "format_version must be 1.");
    }
}

std::optional<std::vector<std::size_t>> parse_phase_space_node_counts(ParsedFile& parsed, const Entry& entry, const std::string& source_name,
                                                                      const std::string& message)
{
    std::vector<std::size_t> node_counts;
    if (!parse_size_list(entry, node_counts))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line, message);
        return std::nullopt;
    }
    return node_counts;
}

void validate_phase_space_names(ParsedFile& parsed, const std::string& source_name, const std::size_t dimensions)
{
    if (!validate_names(parsed, dimensions))
    {
        const Entry& names = parsed.scalars.at("phase_space.names");
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, names.line,
                  "Dimension names must contain one nonempty name per phase-space dimension.");
    }
}

void validate_qmom_phase_space(ParsedFile& parsed, const std::string& source_name, const Entry& node_counts_entry,
                               const std::vector<std::size_t>& node_counts)
{
    // Formulation selection is phase-space driven: one dimension is QMOM.
    if (node_counts.size() != 1 || node_counts[0] == 0)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, node_counts_entry.line,
                  "QMOM node_counts must contain exactly one positive integer.");
        return;
    }
    validate_phase_space_names(parsed, source_name, 1);
}

void validate_cqmom_phase_space(ParsedFile& parsed, const std::string& source_name, const Entry& node_counts_entry,
                                const std::vector<std::size_t>& node_counts)
{
    // Multivariate phase space selects CQMOM; one-dimensional cases remain the
    // QMOM formulation rather than a degenerate CQMOM hierarchy.
    if (std::find(node_counts.begin(), node_counts.end(), 0) != node_counts.end())
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, node_counts_entry.line, "All node counts must be positive.");
    }
    if (node_counts.size() == 1)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, node_counts_entry.line,
                  "CQMOM configuration requires a multivariate phase space; one-dimensional phase space is interpreted as QMOM.");
    }
    if (!node_counts.empty())
    {
        validate_phase_space_names(parsed, source_name, node_counts.size());
    }
}

void validate_optional_qmom_phase_space(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* node_counts_entry = find_entry(parsed, "phase_space.node_counts");
    if (node_counts_entry == nullptr)
    {
        return;
    }

    const std::optional<std::vector<std::size_t>> node_counts =
        parse_phase_space_node_counts(parsed, *node_counts_entry, source_name, "QMOM node_counts must contain exactly one positive integer.");
    if (node_counts.has_value())
    {
        validate_qmom_phase_space(parsed, source_name, *node_counts_entry, *node_counts);
    }
}

void validate_optional_cqmom_phase_space(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* node_counts_entry = find_entry(parsed, "phase_space.node_counts");
    if (node_counts_entry == nullptr)
    {
        return;
    }

    const std::optional<std::vector<std::size_t>> node_counts =
        parse_phase_space_node_counts(parsed, *node_counts_entry, source_name, "node_counts must contain positive integers.");
    if (node_counts.has_value())
    {
        validate_cqmom_phase_space(parsed, source_name, *node_counts_entry, *node_counts);
    }
}

bool parse_required_double(ParsedFile& parsed, const std::string& key, const std::string& source_name, double& value)
{
    const Entry* entry = require_entry(parsed, key, source_name);
    if (entry == nullptr)
    {
        return false;
    }
    if (!parse_double(entry->value, value))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, "Key '" + key + "' must be a finite number.");
        return false;
    }
    return true;
}

bool parse_optional_double(ParsedFile& parsed, const std::string& key, const std::string& source_name, double& value)
{
    const Entry* entry = find_entry(parsed, key);
    if (entry == nullptr)
    {
        return true;
    }
    if (!parse_double(entry->value, value))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, "Key '" + key + "' must be a finite number.");
        return false;
    }
    return true;
}

bool parse_optional_bool(ParsedFile& parsed, const std::string& key, const std::string& source_name, bool& value)
{
    const Entry* entry = find_entry(parsed, key);
    if (entry == nullptr)
    {
        return true;
    }
    if (entry->value == "true")
    {
        value = true;
        return true;
    }
    if (entry->value == "false")
    {
        value = false;
        return true;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, "Key '" + key + "' must be true or false.");
    return false;
}

bool parse_optional_size(ParsedFile& parsed, const std::string& key, const std::string& source_name, std::size_t& value)
{
    const Entry* entry = find_entry(parsed, key);
    if (entry == nullptr)
    {
        return true;
    }
    if (!parse_size(entry->value, value))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, "Key '" + key + "' must be a non-negative integer.");
        return false;
    }
    return true;
}

bool parse_optional_positive_double(ParsedFile& parsed, const std::string& key, const std::string& source_name, double& value)
{
    const Entry* entry = find_entry(parsed, key);
    if (entry == nullptr)
    {
        return true;
    }
    double parsed_value = 0.0;
    if (!parse_double(entry->value, parsed_value) || parsed_value <= 0.0)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, key + " must be a positive finite number.");
        return false;
    }
    value = parsed_value;
    return true;
}

std::optional<std::size_t> infer_qmom_initial_maximum_order(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* node_counts_entry = require_entry(parsed, "phase_space.node_counts", source_name);
    if (node_counts_entry == nullptr)
    {
        return std::nullopt;
    }
    std::optional<std::vector<std::size_t>> node_counts =
        parse_phase_space_node_counts(parsed, *node_counts_entry, source_name, "QMOM node_counts must contain exactly one positive integer.");
    if (!node_counts.has_value())
    {
        return std::nullopt;
    }
    validate_qmom_phase_space(parsed, source_name, *node_counts_entry, *node_counts);
    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    if ((*node_counts)[0] > std::numeric_limits<std::size_t>::max() / 2)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, node_counts_entry->line,
                  "QMOM node_count is too large to infer initial.maximum_order.");
        return std::nullopt;
    }
    return 2 * (*node_counts)[0] - 1;
}

bool parse_optional_positive_size(ParsedFile& parsed, const std::string& key, const std::string& source_name, std::size_t& value)
{
    const Entry* entry = find_entry(parsed, key);
    if (entry == nullptr)
    {
        return true;
    }
    std::size_t parsed_value = 0;
    if (!parse_size(entry->value, parsed_value) || parsed_value == 0)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, key + " must be a positive integer.");
        return false;
    }
    value = parsed_value;
    return true;
}

std::optional<DistributionSpecification> parse_distribution(ParsedFile& parsed, const std::string& key_prefix, const std::string& source_name)
{
    const Entry* family_entry = require_entry(parsed, key_prefix + ".family", source_name);
    if (family_entry == nullptr)
    {
        return std::nullopt;
    }

    const std::string family = family_entry->value;
    if (family == "normal")
    {
        double mean = 0.0;
        double standard_deviation = 0.0;
        if (parse_required_double(parsed, key_prefix + ".mean", source_name, mean) &&
            parse_required_double(parsed, key_prefix + ".standard_deviation", source_name, standard_deviation))
        {
            return NormalDistribution{mean, standard_deviation};
        }
        return std::nullopt;
    }
    if (family == "lognormal" || family == "log-normal" || family == "log_normal")
    {
        double normal_mean = 0.0;
        double normal_standard_deviation = 0.0;
        if (parse_required_double(parsed, key_prefix + ".mu_log", source_name, normal_mean) &&
            parse_required_double(parsed, key_prefix + ".sigma_log", source_name, normal_standard_deviation))
        {
            return LogNormalDistribution{normal_mean, normal_standard_deviation};
        }
        return std::nullopt;
    }
    if (family == "beta")
    {
        double alpha = 0.0;
        double beta = 0.0;
        if (parse_required_double(parsed, key_prefix + ".alpha", source_name, alpha) &&
            parse_required_double(parsed, key_prefix + ".beta", source_name, beta))
        {
            return BetaDistribution{alpha, beta};
        }
        return std::nullopt;
    }
    if (family == "gamma")
    {
        double shape = 0.0;
        double scale = 0.0;
        if (parse_required_double(parsed, key_prefix + ".shape", source_name, shape) &&
            parse_required_double(parsed, key_prefix + ".scale", source_name, scale))
        {
            return GammaDistribution{shape, scale};
        }
        return std::nullopt;
    }
    if (family == "chi")
    {
        double degrees_of_freedom = 0.0;
        if (parse_required_double(parsed, key_prefix + ".degrees_of_freedom", source_name, degrees_of_freedom))
        {
            return ChiDistribution{degrees_of_freedom};
        }
        return std::nullopt;
    }
    if (family == "chi-squared" || family == "chi_squared" || family == "chisquared")
    {
        double degrees_of_freedom = 0.0;
        if (parse_required_double(parsed, key_prefix + ".degrees_of_freedom", source_name, degrees_of_freedom))
        {
            return ChiSquaredDistribution{degrees_of_freedom};
        }
        return std::nullopt;
    }
    if (family == "weibull")
    {
        double shape = 0.0;
        double scale = 0.0;
        if (parse_required_double(parsed, key_prefix + ".shape", source_name, shape) &&
            parse_required_double(parsed, key_prefix + ".scale", source_name, scale))
        {
            return WeibullDistribution{shape, scale};
        }
        return std::nullopt;
    }

    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, family_entry->line, "Unsupported distribution family '" + family + "'.");
    return std::nullopt;
}

std::optional<DistributionSpecification> parse_distribution(ParsedFile& parsed, const std::string& source_name)
{
    return parse_distribution(parsed, "initial", source_name);
}

bool component_key_allowed(const std::string& key, const std::string& prefix)
{
    static const std::vector<std::string> fields{"family", "number_density", "mean",  "standard_deviation", "mu_log", "sigma_log", "alpha",
                                                 "beta",   "shape",          "scale", "degrees_of_freedom"};
    if (key.rfind(prefix + ".", 0) != 0)
    {
        return false;
    }
    const std::string field = key.substr(prefix.size() + 1);
    return std::find(fields.begin(), fields.end(), field) != fields.end();
}

bool marginal_key_allowed(const std::string& key, const std::string& prefix)
{
    static const std::vector<std::string> fields{"type", "family", "mean",  "standard_deviation", "mu_log", "sigma_log", "alpha",
                                                 "beta", "shape",  "scale", "degrees_of_freedom"};
    if (key.rfind(prefix + ".", 0) != 0)
    {
        return false;
    }
    if (key.rfind(mixture_component_prefix(prefix), 0) == 0)
    {
        return true;
    }
    const std::string field = key.substr(prefix.size() + 1);
    return std::find(fields.begin(), fields.end(), field) != fields.end();
}

bool has_explicit_qmom_moments(const ParsedFile& parsed)
{
    return find_entry(parsed, "moments.values") != nullptr;
}

bool has_explicit_cqmom_moments(const ParsedFile& parsed)
{
    return find_entry(parsed, "moments.extents") != nullptr || !parsed.moments.empty();
}

std::vector<std::string> qmom_allowed_keys()
{
    return {"format_version",
            "phase_space.names",
            "phase_space.node_counts",
            "phase_space.support_lower",
            "phase_space.support_upper",
            "qmom.adaptive",
            "moments.values",
            "initial.type",
            "initial.family",
            "initial.number_density",
            "initial.maximum_order",
            "initial.mean",
            "initial.standard_deviation",
            "initial.mu_log",
            "initial.sigma_log",
            "initial.alpha",
            "initial.beta",
            "initial.shape",
            "initial.scale",
            "initial.degrees_of_freedom"};
}

void append_time_integration_allowed_keys(std::vector<std::string>& allowed, const RealizabilitySettings realizability)
{
    allowed.push_back(std::string(phase_space_evolution_scheme_key));
    allowed.push_back(std::string(node_advection_scheme_key));
    allowed.push_back(std::string(adaptive_relative_tolerance_key));
    allowed.push_back(std::string(adaptive_absolute_tolerance_key));
    allowed.push_back(std::string(adaptive_initial_step_key));
    allowed.push_back(std::string(adaptive_max_step_key));
    allowed.push_back(std::string(adaptive_min_step_key));
    allowed.push_back(std::string(adaptive_max_steps_key));
    if (!accepts_unnumbered_realizability(realizability))
    {
        return;
    }
    allowed.push_back("time_integration.realizability.check_mode");
    allowed.push_back("time_integration.realizability.accept_boundary");
    allowed.push_back("time_integration.realizability.support");
    allowed.push_back("time_integration.realizability.tolerance");
    allowed.push_back("time_integration.realizability.relative_tolerance");
    allowed.push_back("time_integration.realizability.lower");
    allowed.push_back("time_integration.realizability.upper");
}

void append_output_allowed_keys(std::vector<std::string>& allowed)
{
    allowed.push_back("output.directory");
    allowed.push_back("output.write_moments");
    allowed.push_back("output.write_quadrature");
    allowed.push_back("output.interval");
    allowed.push_back("output.generalized_moments");
}

std::vector<std::string> with_runtime_allowed_keys(std::vector<std::string> allowed, const RealizabilitySettings realizability)
{
    append_time_integration_allowed_keys(allowed, realizability);
    append_output_allowed_keys(allowed);
    return allowed;
}

std::vector<std::string> time_integration_sections(const RealizabilitySettings realizability)
{
    std::vector<std::string> sections{"time_integration", "time_integration.node_advection", "time_integration.node_advection.adaptive"};
    if (accepts_unnumbered_realizability(realizability))
    {
        sections.emplace_back(realizability_section);
    }
    return sections;
}

void reject_unknown_qmom_sections(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    std::vector<std::string> allowed = {"phase_space", "moments", "qmom", "initial", "output"};
    const std::vector<std::string> time_sections = time_integration_sections(RealizabilitySettings::qmom);
    allowed.insert(allowed.end(), time_sections.begin(), time_sections.end());
    if (explicit_moments)
    {
        reject_unknown_sections(parsed, allowed, {"initial."}, RealizabilitySettings::qmom, source_name);
        return;
    }
    reject_unknown_sections(parsed, allowed, {std::string(mixture_component_section_prefix)}, RealizabilitySettings::qmom, source_name);
}

void reject_unknown_qmom_keys(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    if (explicit_moments)
    {
        reject_unknown_keys(parsed,
                            with_runtime_allowed_keys({"format_version", "phase_space.names", "phase_space.node_counts", "phase_space.support_lower",
                                                       "phase_space.support_upper", "qmom.adaptive", "moments.values"},
                                                      RealizabilitySettings::qmom),
                            {"initial."}, RealizabilitySettings::qmom, source_name);
        return;
    }
    reject_unknown_keys(parsed, with_runtime_allowed_keys(qmom_allowed_keys(), RealizabilitySettings::qmom),
                        {std::string(mixture_component_section_prefix)}, RealizabilitySettings::qmom, source_name);
}

// Rejects a re-declared indexed section (e.g. "initial.dimension.N" or
// "initial.mixture.component.N"). Deliberately scoped to indexed sections only:
// ordinary sections (e.g. "[phase_space]") may legitimately be reopened to add
// keys incrementally, so re-declaring one of those is not itself an error.
void reject_duplicate_indexed_section(ParsedFile& parsed, const std::string& section, const std::string& source_name)
{
    const auto lines = parsed.section_lines.find(section);
    if (lines != parsed.section_lines.end() && lines->second.size() > 1)
    {
        add_error(parsed, ConfigurationErrorCode::duplicate_entry, source_name, lines->second[1], "Duplicate section '" + section + "'.");
    }
}

// Reports every section and key under an indexed prefix, such as
// "initial.mixture.component.", that does not belong to an indexed section the
// initialization actually read: a section for another initialization type, an
// index without a section of its own, or an index that is not canonical. A key
// belongs to the section that is its direct parent, so
// "initial.dimension.0.mixture.component.7.mean" needs component 7 itself.
void reject_unconsumed_indexed_entries(ParsedFile& parsed, const std::string& indexed_prefix, const std::string& source_name)
{
    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind(indexed_prefix, 0) == 0 && parsed.consumed_indexed_sections.count(item.first) == 0)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind(indexed_prefix, 0) != 0)
        {
            continue;
        }
        const std::string parent = item.first.substr(0, item.first.rfind('.'));
        if (parsed.consumed_indexed_sections.count(parent) == 0)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }
}

std::optional<DistributionMixture> parse_mixture(ParsedFile& parsed, const std::string& source_name,
                                                 const std::string& initialization_prefix = "initial")
{
    const std::string component_section_prefix = mixture_component_prefix(initialization_prefix);
    std::map<std::size_t, Entry> component_sections;
    for (const auto& item : parsed.sections)
    {
        if (!is_mixture_component_section(item.first, initialization_prefix))
        {
            continue;
        }
        reject_duplicate_indexed_section(parsed, item.first, source_name);
        std::size_t component = 0;
        parse_size(item.first.substr(component_section_prefix.size()), component);
        component_sections.emplace(component, item.second);
    }

    if (component_sections.empty())
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "Mixture initialization requires at least one '[" + component_section_prefix + "N]' section.");
        return std::nullopt;
    }

    std::vector<DistributionComponent> components;
    components.reserve(component_sections.size());
    std::size_t expected_component = 0;
    for (const auto& item : component_sections)
    {
        if (item.first != expected_component)
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, item.second.line,
                      "Mixture component indices must be contiguous starting at zero.");
            return std::nullopt;
        }

        const std::string prefix = component_section_prefix + std::to_string(item.first);
        parsed.consumed_indexed_sections.insert(prefix);
        for (const auto& scalar : parsed.scalars)
        {
            if (scalar.first.rfind(prefix + ".", 0) == 0 && !component_key_allowed(scalar.first, prefix))
            {
                add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, scalar.second.line, "Unknown key '" + scalar.first + "'.");
            }
        }

        std::optional<DistributionSpecification> distribution = parse_distribution(parsed, prefix, source_name);
        double number_density = 0.0;
        const bool parsed_density = parse_required_double(parsed, prefix + ".number_density", source_name, number_density);
        if (parsed_density && number_density < 0.0)
        {
            const Entry* entry = find_entry(parsed, prefix + ".number_density");
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? item.second.line : entry->line,
                      "Mixture component number_density must be finite and non-negative.");
        }

        if (distribution.has_value())
        {
            components.push_back({*distribution, number_density});
        }
        ++expected_component;
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return DistributionMixture{std::move(components)};
}

bool has_mixture_component_sections(const ParsedFile& parsed, const std::string& initialization_prefix = "initial")
{
    for (const auto& item : parsed.sections)
    {
        if (is_mixture_component_section(item.first, initialization_prefix))
        {
            return true;
        }
    }
    return false;
}

std::vector<std::string> cqmom_explicit_allowed_keys()
{
    return with_runtime_allowed_keys({"format_version", "phase_space.names", "phase_space.node_counts", "phase_space.conditioning_order",
                                      "phase_space.adaptive_branch_counts", "phase_space.support_lower", "phase_space.support_upper", "qmom.adaptive",
                                      "cqmom.diagnostic_level", "vandermonde.solver", "vandermonde.rcond", "moments.extents"},
                                     RealizabilitySettings::cqmom);
}

std::vector<std::string> cqmom_initial_allowed_keys()
{
    std::vector<std::string> allowed = cqmom_explicit_allowed_keys();
    allowed.push_back("initial.type");
    allowed.push_back("initial.number_density");
    allowed.push_back("initial.maximum_orders");
    return allowed;
}

void reject_unknown_cqmom_sections(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    std::vector<std::string> allowed = {"phase_space", "moments", "qmom", "cqmom", "vandermonde", "initial", "output"};
    const std::vector<std::string> time_sections = time_integration_sections(RealizabilitySettings::cqmom);
    allowed.insert(allowed.end(), time_sections.begin(), time_sections.end());
    if (explicit_moments)
    {
        reject_unknown_sections(parsed, allowed, {"initial."}, RealizabilitySettings::cqmom, source_name);
        return;
    }
    reject_unknown_sections(parsed, allowed, {std::string(independent_product_dimension_section_prefix)}, RealizabilitySettings::cqmom, source_name);
}

void reject_unknown_cqmom_keys(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    // `[qmom] adaptive` decides both the inversions and the branch counts. The
    // former key stays recognised so that a file using it gets this message
    // rather than a bare unknown-key error.
    if (const Entry* removed = find_entry(parsed, "phase_space.adaptive_branch_counts"); removed != nullptr)
    {
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, removed->line,
                  "phase_space.adaptive_branch_counts was removed; use [qmom] adaptive.");
    }
    if (explicit_moments)
    {
        reject_unknown_keys(parsed, cqmom_explicit_allowed_keys(), {"initial."}, RealizabilitySettings::cqmom, source_name);
        return;
    }
    reject_unknown_keys(parsed, cqmom_initial_allowed_keys(), {std::string(independent_product_dimension_section_prefix)},
                        RealizabilitySettings::cqmom, source_name);
}

std::optional<IndependentProductDistribution> parse_independent_product(ParsedFile& parsed, const std::string& source_name)
{
    std::map<std::size_t, Entry> dimension_sections;
    for (const auto& item : parsed.sections)
    {
        if (!is_independent_product_dimension_section(item.first))
        {
            continue;
        }
        reject_duplicate_indexed_section(parsed, item.first, source_name);
        std::size_t dimension = 0;
        parse_size(item.first.substr(independent_product_dimension_section_prefix.size()), dimension);
        dimension_sections.emplace(dimension, item.second);
    }

    if (dimension_sections.empty())
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "Independent product initialization requires at least one '[initial.dimension.N]' section.");
        return std::nullopt;
    }

    double number_density = 0.0;
    const bool parsed_density = parse_required_double(parsed, "initial.number_density", source_name, number_density);
    if (parsed_density && number_density < 0.0)
    {
        const Entry* entry = find_entry(parsed, "initial.number_density");
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? 0 : entry->line,
                  "initial.number_density must be finite and non-negative.");
    }

    std::vector<IndependentProductMarginal> marginals;
    marginals.reserve(dimension_sections.size());
    std::size_t expected_dimension = 0;
    for (const auto& item : dimension_sections)
    {
        if (item.first != expected_dimension)
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, item.second.line,
                      "Independent product dimension indices must be contiguous starting at zero.");
            return std::nullopt;
        }

        const std::string prefix = std::string(independent_product_dimension_section_prefix) + std::to_string(item.first);
        parsed.consumed_indexed_sections.insert(prefix);
        for (const auto& scalar : parsed.scalars)
        {
            if (scalar.first.rfind(prefix + ".", 0) == 0 && !marginal_key_allowed(scalar.first, prefix))
            {
                add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, scalar.second.line, "Unknown key '" + scalar.first + "'.");
            }
        }

        const Entry* marginal_type = find_entry(parsed, prefix + ".type");
        const bool has_components = has_mixture_component_sections(parsed, prefix);
        if (marginal_type != nullptr)
        {
            if (marginal_type->value != "mixture")
            {
                add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, marginal_type->line,
                          "Independent product marginal type must be 'mixture' when specified.");
            }
            else
            {
                std::optional<DistributionMixture> mixture = parse_mixture(parsed, source_name, prefix);
                if (mixture.has_value())
                {
                    marginals.push_back(*mixture);
                }
            }
        }
        else if (has_components)
        {
            add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, item.second.line,
                      "Independent product mixture marginal requires 'type = mixture'.");
        }
        else
        {
            std::optional<DistributionSpecification> distribution = parse_distribution(parsed, prefix, source_name);
            if (distribution.has_value())
            {
                marginals.push_back(*distribution);
            }
        }
        ++expected_dimension;
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return IndependentProductDistribution{std::move(marginals), number_density};
}

bool checked_product(const std::vector<std::size_t>& extents, std::size_t& product)
{
    const std::optional<std::size_t> value_count = internal::raw_moment_tensor_value_count(extents);
    if (!value_count.has_value())
    {
        return false;
    }
    product = *value_count;
    return true;
}

std::optional<RawMomentTensor> parse_explicit_tensor(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* extents_entry = require_entry(parsed, "moments.extents", source_name);
    std::vector<std::size_t> extents;
    if (extents_entry != nullptr && !parse_size_list(*extents_entry, extents))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, extents_entry->line, "extents must contain positive integers.");
    }

    std::size_t moment_count = 0;
    if (!extents.empty() && !checked_product(extents, moment_count))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, extents_entry == nullptr ? 0 : extents_entry->line,
                  "Moment extents are zero or their product overflows.");
    }

    // Every entry is checked, and the entries inside the extents are counted,
    // before the tensor is allocated: declared extents far beyond the supplied
    // entries are reported as an incomplete tensor rather than allocated. Each
    // index appears once, so the count tells whether the tensor is complete.
    std::size_t supplied_count = 0;
    for (const auto& item : parsed.moments)
    {
        if (item.first.size() != extents.size())
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, item.second.line,
                      "Moment index rank does not match the number of entries in moments.extents.");
            continue;
        }
        bool in_range = true;
        for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
        {
            in_range = in_range && item.first[dimension] < extents[dimension];
        }
        double value = 0.0;
        if (!in_range)
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, item.second.line,
                      "Moment index is outside the declared extents.");
        }
        else if (!parse_double(item.second.value, value))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, item.second.line, "Moment value must be a finite number.");
        }
        else
        {
            ++supplied_count;
        }
    }

    if (moment_count != 0 && supplied_count < moment_count)
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "The full tensor of raw moments is incomplete; it is missing " + std::to_string(moment_count - supplied_count) + " of its " +
                      std::to_string(moment_count) + " entries.");
    }
    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }

    std::vector<double> values(moment_count, 0.0);
    for (const auto& item : parsed.moments)
    {
        static_cast<void>(parse_double(item.second.value, values[flatten_multi_index(item.first, extents)]));
    }
    return RawMomentTensor{std::move(extents), std::move(values)};
}

std::string format_configuration_error(const ConfigurationError& error)
{
    std::ostringstream message;
    if (!error.source_name.empty())
    {
        message << error.source_name;
        if (error.line != 0)
        {
            message << ':' << error.line;
        }
        message << ": ";
    }
    message << error.message;
    return message.str();
}

std::string format_configuration_errors(const std::vector<ConfigurationError>& errors)
{
    if (errors.empty())
    {
        return "Configuration parsing failed.";
    }
    return format_configuration_error(errors.front());
}

template <class T> ConfigurationResult<T> failure_result(ParsedFile&& parsed)
{
    return ConfigurationResult<T>::failure(format_configuration_errors(parsed.errors));
}

template <class T> ConfigurationResult<T> failure_result(ConfigurationError error)
{
    return ConfigurationResult<T>::failure(format_configuration_error(error));
}

template <class T> ConfigurationResult<T> file_failure(const std::filesystem::path& filename)
{
    return failure_result<T>({ConfigurationErrorCode::file_error, filename.string(), 0, "Unable to open configuration file."});
}

[[nodiscard]] std::optional<std::string> validate_replacement_tensor_structure(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy)
{
    if (moments.extents.size() != hierarchy.number_of_dimensions())
    {
        return "The rank of the full tensor of raw moments must match the base CQMOM hierarchy.";
    }

    const std::optional<std::size_t> value_count = internal::raw_moment_tensor_value_count(moments.extents);
    if (!value_count.has_value())
    {
        return "The extents of the full tensor of raw moments must be positive and their product must not overflow.";
    }
    if (moments.values.size() != *value_count)
    {
        return "The value count of the full tensor of raw moments does not match its extents.";
    }
    if (!internal::all_finite(moments.values))
    {
        return "The values of the full tensor of raw moments must be finite.";
    }
    return std::nullopt;
}

[[nodiscard]] std::string format_moment_index(const MultiIndexView index)
{
    std::ostringstream message;
    message << '[';
    for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
    {
        if (dimension != 0)
        {
            message << ", ";
        }
        message << index[dimension];
    }
    message << ']';
    return message.str();
}

std::optional<std::string> read_configuration_file(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return std::nullopt;
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

struct ParsedCqmomCaseSettings
{
    PhaseSpaceHierarchy hierarchy;
    QmomOptions qmom_options;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
    double moore_penrose_rcond = 1.0e-14;
    CqmomAbscissaBounds coordinate_support;
};

bool has_realizability_section(const ParsedFile& parsed)
{
    if (parsed.sections.find(std::string(realizability_section)) != parsed.sections.end())
    {
        return true;
    }
    return std::any_of(parsed.scalars.begin(), parsed.scalars.end(),
                       [](const auto& item) { return item.first.rfind(std::string(realizability_key_prefix), 0) == 0; });
}

bool has_unnumbered_realizability_section(const ParsedFile& parsed)
{
    if (parsed.sections.find(std::string(realizability_section)) != parsed.sections.end())
    {
        return true;
    }
    return std::any_of(parsed.scalars.begin(), parsed.scalars.end(), [](const auto& item)
                       { return item.first.rfind(std::string(realizability_key_prefix), 0) == 0 && !is_coordinate_realizability_key(item.first); });
}

std::vector<std::pair<std::size_t, std::string>> coordinate_realizability_sections(const ParsedFile& parsed)
{
    std::vector<std::pair<std::size_t, std::string>> sections;
    for (const auto& item : parsed.sections)
    {
        std::size_t coordinate = 0;
        if (parse_coordinate_realizability_section(item.first, coordinate))
        {
            sections.emplace_back(coordinate, item.first);
        }
    }
    return sections;
}

bool is_realizability_key(const std::string& key)
{
    return key.rfind(std::string(realizability_key_prefix), 0) == 0;
}

const std::vector<std::string>& realizability_policy_fields()
{
    static const std::vector<std::string> fields{"check_mode", "accept_boundary", "support", "tolerance", "relative_tolerance", "lower", "upper"};
    return fields;
}

bool is_node_advection_key(const std::string& key)
{
    return key.rfind(std::string(node_advection_key_prefix), 0) == 0;
}

const std::vector<std::string>& node_advection_allowed_keys()
{
    static const std::vector<std::string> allowed{std::string(node_advection_scheme_key),       std::string(adaptive_relative_tolerance_key),
                                                  std::string(adaptive_absolute_tolerance_key), std::string(adaptive_initial_step_key),
                                                  std::string(adaptive_max_step_key),           std::string(adaptive_min_step_key),
                                                  std::string(adaptive_max_steps_key)};
    return allowed;
}

std::optional<QmomOptions> parse_qmom_options_from_parsed(ParsedFile& parsed, const std::string& source_name)
{
    QmomOptions options;
    const auto adaptive_entry = parsed.scalars.find("qmom.adaptive");
    if (adaptive_entry != parsed.scalars.end())
    {
        if (adaptive_entry->second.value == "true")
        {
            options.adaptive = true;
        }
        else if (adaptive_entry->second.value == "false")
        {
            options.adaptive = false;
        }
        else
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, adaptive_entry->second.line,
                      "qmom.adaptive must be true or false.");
        }
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return options;
}

void reject_unknown_qmom_option_keys(ParsedFile& parsed, const std::string& source_name)
{
    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind("qmom.", 0) == 0 && item.first != "qmom.adaptive")
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }
    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("qmom.", 0) == 0)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
}

std::optional<RealizabilityCheckMode> parse_realizability_check_mode_value(ParsedFile& parsed, const Entry& entry, const std::string& key,
                                                                           const std::string& source_name)
{
    if (entry.value == "none")
    {
        return RealizabilityCheckMode::None;
    }
    if (entry.value == "every_stage")
    {
        return RealizabilityCheckMode::EveryStage;
    }
    if (entry.value == "final_only")
    {
        return RealizabilityCheckMode::FinalOnly;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line, key + " must be none, every_stage, or final_only.");
    return std::nullopt;
}

std::optional<RealizabilityCheckMode> parse_realizability_check_mode(ParsedFile& parsed, const std::string& key, const std::string& source_name)
{
    const Entry* check_mode = find_entry(parsed, key);
    if (check_mode != nullptr)
    {
        return parse_realizability_check_mode_value(parsed, *check_mode, key, source_name);
    }

    return RealizabilityCheckMode::None;
}

std::optional<MomentSupport> parse_moment_support(ParsedFile& parsed, const std::string& key, const std::string& source_name)
{
    const Entry* entry = require_entry(parsed, key, source_name);
    if (entry == nullptr)
    {
        return std::nullopt;
    }
    if (entry->value == "real_line")
    {
        return MomentSupport::real_line;
    }
    if (entry->value == "nonnegative")
    {
        return MomentSupport::nonnegative;
    }
    if (entry->value == "bounded_interval")
    {
        return MomentSupport::bounded_interval;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line, key + " must be real_line, nonnegative, or bounded_interval.");
    return std::nullopt;
}

std::optional<PhaseSpaceEvolutionScheme> parse_phase_space_evolution_scheme_value(ParsedFile& parsed, const Entry& entry,
                                                                                  const std::string& source_name)
{
    if (entry.value == "strang")
    {
        return PhaseSpaceEvolutionScheme::strang;
    }
    if (entry.value == "advsrc")
    {
        return PhaseSpaceEvolutionScheme::advsrc;
    }
    if (entry.value == "srcadv")
    {
        return PhaseSpaceEvolutionScheme::srcadv;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line, "time_integration.scheme must be strang, advsrc, or srcadv.");
    return std::nullopt;
}

std::optional<NodeAdvectionScheme> parse_node_advection_scheme_value(ParsedFile& parsed, const Entry& entry, const std::string& source_name)
{
    if (entry.value == "rk23")
    {
        return NodeAdvectionScheme::rk23;
    }
    if (entry.value == "rk45")
    {
        return NodeAdvectionScheme::rk45;
    }
    if (entry.value == "ssp_rk3")
    {
        return NodeAdvectionScheme::ssp_rk3;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line,
              "time_integration.node_advection.scheme must be rk23, rk45, or ssp_rk3.");
    return std::nullopt;
}

std::optional<NodeAdvectionOptions> parse_node_advection_options_from_parsed(ParsedFile& parsed, const std::string& source_name)
{
    for (const auto& item : parsed.scalars)
    {
        if (!is_node_advection_key(item.first))
        {
            continue;
        }
        const std::vector<std::string>& allowed = node_advection_allowed_keys();
        if (std::find(allowed.begin(), allowed.end(), item.first) == allowed.end())
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }

    NodeAdvectionOptions options;
    const Entry* scheme = find_entry(parsed, std::string(node_advection_scheme_key));
    if (scheme != nullptr)
    {
        const std::optional<NodeAdvectionScheme> parsed_scheme = parse_node_advection_scheme_value(parsed, *scheme, source_name);
        if (parsed_scheme.has_value())
        {
            options.scheme = *parsed_scheme;
        }
    }

    parse_optional_positive_double(parsed, std::string(adaptive_relative_tolerance_key), source_name, options.adaptive_options.relative_tolerance);
    parse_optional_positive_double(parsed, std::string(adaptive_absolute_tolerance_key), source_name, options.adaptive_options.absolute_tolerance);
    parse_optional_positive_double(parsed, std::string(adaptive_initial_step_key), source_name, options.adaptive_options.initial_step);
    parse_optional_positive_double(parsed, std::string(adaptive_max_step_key), source_name, options.adaptive_options.max_step);
    parse_optional_positive_double(parsed, std::string(adaptive_min_step_key), source_name, options.adaptive_options.min_step);
    parse_optional_positive_size(parsed, std::string(adaptive_max_steps_key), source_name, options.adaptive_options.max_steps);

    const Entry* min_step = find_entry(parsed, std::string(adaptive_min_step_key));
    const Entry* max_step = find_entry(parsed, std::string(adaptive_max_step_key));
    if (min_step != nullptr && max_step != nullptr && parsed.errors.empty() && options.adaptive_options.min_step > options.adaptive_options.max_step)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, min_step->line,
                  "time_integration.node_advection.adaptive.min_step must not exceed max_step.");
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return options;
}

std::optional<RealizabilityPolicy> parse_realizability_policy_from_parsed(ParsedFile& parsed, const std::string& section,
                                                                          const std::string& source_name)
{
    const std::string key_prefix = section + ".";

    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind(key_prefix, 0) != 0 || (section == realizability_section && is_coordinate_realizability_key(item.first)))
        {
            continue;
        }
        const std::string field = item.first.substr(key_prefix.size());
        if (std::find(realizability_policy_fields().begin(), realizability_policy_fields().end(), field) == realizability_policy_fields().end())
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }

    const std::string check_mode_key = key_prefix + "check_mode";
    const std::string support_key = key_prefix + "support";
    const std::string accept_boundary_key = key_prefix + "accept_boundary";
    const std::string tolerance_key = key_prefix + "tolerance";
    const std::string relative_tolerance_key = key_prefix + "relative_tolerance";
    const std::string lower_key = key_prefix + "lower";
    const std::string upper_key = key_prefix + "upper";

    std::optional<RealizabilityCheckMode> check_mode = parse_realizability_check_mode(parsed, check_mode_key, source_name);
    if (!check_mode.has_value())
    {
        return std::nullopt;
    }

    RealizabilityPolicy policy;
    policy.check_mode = *check_mode;
    parse_optional_bool(parsed, accept_boundary_key, source_name, policy.accept_boundary);
    parse_optional_double(parsed, tolerance_key, source_name, policy.options.tolerance);
    parse_optional_bool(parsed, relative_tolerance_key, source_name, policy.options.relative_tolerance);
    if (policy.options.tolerance <= 0.0)
    {
        const Entry* entry = find_entry(parsed, tolerance_key);
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? 0 : entry->line,
                  tolerance_key + " must be positive.");
    }

    if (policy.check_mode != RealizabilityCheckMode::None || find_entry(parsed, support_key) != nullptr)
    {
        std::optional<MomentSupport> support = parse_moment_support(parsed, support_key, source_name);
        if (support.has_value())
        {
            policy.support = *support;
        }
    }

    if (policy.support == MomentSupport::bounded_interval)
    {
        double lower = 0.0;
        double upper = 0.0;
        const bool parsed_lower = parse_required_double(parsed, lower_key, source_name, lower);
        const bool parsed_upper = parse_required_double(parsed, upper_key, source_name, upper);
        if (parsed_lower && parsed_upper)
        {
            if (lower >= upper)
            {
                const Entry* entry = find_entry(parsed, lower_key);
                add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? 0 : entry->line,
                          lower_key + " must be less than upper.");
            }
            else
            {
                policy.interval = BoundedInterval{lower, upper};
            }
        }
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return policy;
}

void validate_cqmom_coordinate_realizability_ranks(ParsedFile& parsed, const std::string& source_name, std::size_t dimensions);

std::vector<std::optional<RealizabilityPolicy>> parse_coordinate_realizability_policies_from_parsed(ParsedFile& parsed,
                                                                                                    const std::string& source_name)
{
    // A configuration that declares its phase space bounds each coordinate
    // index by the rank, as the case loaders do, so an index beyond it is
    // reported without sizing the policy vector by it.
    std::optional<std::size_t> rank;
    std::vector<std::size_t> node_counts;
    if (const Entry* node_counts_entry = find_entry(parsed, "phase_space.node_counts");
        node_counts_entry != nullptr && parse_size_list(*node_counts_entry, node_counts))
    {
        rank = node_counts.size();
        validate_cqmom_coordinate_realizability_ranks(parsed, source_name, *rank);
    }

    const std::vector<std::pair<std::size_t, std::string>> sections = coordinate_realizability_sections(parsed);
    std::vector<std::optional<RealizabilityPolicy>> policies;
    for (const auto& item : sections)
    {
        if (rank.has_value() && item.first >= *rank)
        {
            continue;
        }
        const std::optional<std::size_t> required_size = internal::checked_size_add(item.first, 1);
        if (!required_size.has_value())
        {
            const auto section = parsed.sections.find(item.second);
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, section == parsed.sections.end() ? 0 : section->second.line,
                      "Coordinate realizability policy index cannot be represented.");
            continue;
        }
        const std::optional<RealizabilityPolicy> policy = parse_realizability_policy_from_parsed(parsed, item.second, source_name);
        if (policy.has_value())
        {
            if (policies.size() < *required_size)
            {
                // The policy vector is indexed by coordinate, so an index far
                // beyond any phase-space rank would need an impossible allocation.
                try
                {
                    policies.resize(*required_size);
                }
                catch (const std::exception&)
                {
                    const auto section = parsed.sections.find(item.second);
                    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, section == parsed.sections.end() ? 0 : section->second.line,
                              "Coordinate realizability policy index is too large.");
                    continue;
                }
            }
            policies[item.first] = *policy;
        }
    }
    return policies;
}

std::optional<CqmomDiagnosticLevel> parse_cqmom_diagnostic_level_value(ParsedFile& parsed, const Entry& entry, const std::string& source_name)
{
    if (entry.value == "detailed")
    {
        return CqmomDiagnosticLevel::detailed;
    }
    if (entry.value == "basic")
    {
        return CqmomDiagnosticLevel::basic;
    }
    if (entry.value == "none")
    {
        return CqmomDiagnosticLevel::none;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line, "cqmom.diagnostic_level must be detailed, basic, or none.");
    return std::nullopt;
}

std::optional<CqmomDiagnosticLevel> parse_cqmom_diagnostic_level(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, "cqmom.diagnostic_level");
    if (entry == nullptr)
    {
        return CqmomDiagnosticLevel::none;
    }
    return parse_cqmom_diagnostic_level_value(parsed, *entry, source_name);
}

std::optional<VandermondeSolver> parse_vandermonde_solver_value(ParsedFile& parsed, const Entry& entry, const std::string& source_name)
{
    if (entry.value == "automatic")
    {
        return VandermondeSolver::automatic;
    }
    if (entry.value == "dense")
    {
        return VandermondeSolver::dense;
    }
    if (entry.value == "pivoted_qr")
    {
        return VandermondeSolver::pivoted_qr;
    }
    if (entry.value == "svd")
    {
        return VandermondeSolver::svd;
    }
    if (entry.value == "moore_penrose")
    {
        return VandermondeSolver::moore_penrose;
    }
    if (entry.value == "bjorck_pereyra")
    {
        return VandermondeSolver::bjorck_pereyra;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line,
              "vandermonde.solver must be automatic, dense, pivoted_qr, svd, moore_penrose, or bjorck_pereyra.");
    return std::nullopt;
}

std::optional<double> parse_vandermonde_rcond(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, "vandermonde.rcond");
    if (entry == nullptr)
    {
        return 1.0e-14;
    }

    double rcond = 0.0;
    if (!parse_double(entry->value, rcond) || rcond < 0.0 || rcond >= 1.0)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line,
                  "vandermonde.rcond must be finite and satisfy 0 <= rcond < 1.");
        return std::nullopt;
    }
    return rcond;
}

std::optional<VandermondeSolver> parse_vandermonde_solver(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, "vandermonde.solver");
    if (entry == nullptr)
    {
        return VandermondeSolver::automatic;
    }
    return parse_vandermonde_solver_value(parsed, *entry, source_name);
}

bool has_output_section(const ParsedFile& parsed)
{
    if (parsed.sections.find(std::string(output_section)) != parsed.sections.end())
    {
        return true;
    }
    return std::any_of(parsed.scalars.begin(), parsed.scalars.end(),
                       [](const auto& item) { return item.first.rfind(std::string(output_key_prefix), 0) == 0; });
}

bool is_output_key(const std::string& key)
{
    return key.rfind(std::string(output_key_prefix), 0) == 0;
}

void validate_optional_phase_space_evolution_scheme(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, std::string(phase_space_evolution_scheme_key));
    if (entry != nullptr)
    {
        static_cast<void>(parse_phase_space_evolution_scheme_value(parsed, *entry, source_name));
    }
}

void validate_optional_node_advection_options(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_node_advection_options_from_parsed(parsed, source_name));
}

void validate_optional_realizability_policy(ParsedFile& parsed, const RealizabilitySettings realizability, const std::string& source_name)
{
    if (accepts_unnumbered_realizability(realizability) && has_unnumbered_realizability_section(parsed))
    {
        static_cast<void>(parse_realizability_policy_from_parsed(parsed, std::string(realizability_section), source_name));
    }
    if (!accepts_coordinate_realizability(realizability))
    {
        return;
    }
    for (const auto& item : coordinate_realizability_sections(parsed))
    {
        static_cast<void>(parse_realizability_policy_from_parsed(parsed, item.second, source_name));
    }
}

void validate_cqmom_coordinate_realizability_ranks(ParsedFile& parsed, const std::string& source_name, const std::size_t dimensions)
{
    for (const auto& item : coordinate_realizability_sections(parsed))
    {
        if (item.first >= dimensions)
        {
            const auto section = parsed.sections.find(item.second);
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, section == parsed.sections.end() ? 0 : section->second.line,
                      "CQMOM realizability coordinate index " + std::to_string(item.first) + " is outside the configured phase-space rank " +
                          std::to_string(dimensions) + ".");
        }
    }
}

void reject_unknown_phase_space_evolution_sections(ParsedFile& parsed, const RealizabilitySettings realizability, const std::string& source_name)
{
    for (const auto& item : parsed.sections)
    {
        std::size_t coordinate = 0;
        const bool accepted_realizability_section =
            (accepts_unnumbered_realizability(realizability) && item.first == realizability_section) ||
            (accepts_coordinate_realizability(realizability) && parse_coordinate_realizability_section(item.first, coordinate));
        if (item.first.rfind("time_integration.", 0) == 0 && !accepted_realizability_section && item.first != node_advection_section &&
            item.first != adaptive_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
}

void reject_unknown_phase_space_evolution_keys(ParsedFile& parsed, const RealizabilitySettings realizability, const std::string& source_name)
{
    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind("time_integration.", 0) != 0)
        {
            continue;
        }
        const bool accepted_realizability_key =
            is_realizability_key(item.first) &&
            (is_coordinate_realizability_key(item.first)
                 ? accepts_coordinate_realizability(realizability) && is_configured_coordinate_realizability_key(parsed, item.first)
                 : accepts_unnumbered_realizability(realizability));
        if (item.first == std::string(phase_space_evolution_scheme_key) || accepted_realizability_key || is_node_advection_key(item.first))
        {
            continue;
        }
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
    }
}

std::optional<PhaseSpaceEvolutionScheme> parse_phase_space_evolution_scheme_from_parsed(ParsedFile& parsed, const std::string& source_name)
{
    PhaseSpaceEvolutionScheme scheme = PhaseSpaceEvolutionScheme::strang;
    const Entry* entry = find_entry(parsed, std::string(phase_space_evolution_scheme_key));
    if (entry != nullptr)
    {
        const std::optional<PhaseSpaceEvolutionScheme> parsed_scheme = parse_phase_space_evolution_scheme_value(parsed, *entry, source_name);
        if (parsed_scheme.has_value())
        {
            scheme = *parsed_scheme;
        }
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return scheme;
}

void validate_optional_cqmom_diagnostic_level(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_cqmom_diagnostic_level(parsed, source_name));
}

void validate_optional_vandermonde_solver(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_vandermonde_solver(parsed, source_name));
    static_cast<void>(parse_vandermonde_rcond(parsed, source_name));
}

std::optional<CqmomAbscissaBounds> parse_coordinate_support(ParsedFile& parsed, const std::string& source_name, const std::size_t dimensions)
{
    const Entry* lower_entry = find_entry(parsed, "phase_space.support_lower");
    const Entry* upper_entry = find_entry(parsed, "phase_space.support_upper");
    if (lower_entry == nullptr && upper_entry == nullptr)
    {
        return CqmomAbscissaBounds{};
    }
    if (lower_entry == nullptr || upper_entry == nullptr)
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, lower_entry == nullptr ? upper_entry->line : lower_entry->line,
                  "phase_space.support_lower and phase_space.support_upper must be configured together.");
        return std::nullopt;
    }

    std::vector<std::optional<double>> lower_values;
    std::vector<std::optional<double>> upper_values;
    if (!parse_optional_double_list(*lower_entry, lower_values))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, lower_entry->line,
                  "phase_space.support_lower must contain finite numbers or none entries.");
    }
    if (!parse_optional_double_list(*upper_entry, upper_values))
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, upper_entry->line,
                  "phase_space.support_upper must contain finite numbers or none entries.");
    }
    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    if (lower_values.size() != dimensions || upper_values.size() != dimensions)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name,
                  lower_values.size() != dimensions ? lower_entry->line : upper_entry->line,
                  "phase_space.support_lower and phase_space.support_upper must contain one entry per phase-space dimension.");
        return std::nullopt;
    }

    CqmomAbscissaBounds support(dimensions);
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        if (lower_values[dimension].has_value() != upper_values[dimension].has_value())
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, lower_entry->line,
                      "Each coordinate support must configure both lower and upper or use none for both.");
            return std::nullopt;
        }
        if (!lower_values[dimension].has_value())
        {
            continue;
        }
        if (!(*lower_values[dimension] < *upper_values[dimension]))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, lower_entry->line,
                      "Bounded coordinate support requires lower < upper.");
            return std::nullopt;
        }
        support[dimension] = AbscissaBounds{*lower_values[dimension], *upper_values[dimension]};
    }
    return support;
}

ConfigurationResult<ParsedCqmomCaseSettings> parse_cqmom_case_settings(ParsedFile parsed, std::string source_name)
{
    const bool explicit_moments = has_explicit_cqmom_moments(parsed);
    reject_unknown_cqmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_cqmom_keys(parsed, explicit_moments, source_name);
    validate_optional_phase_space_evolution_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, RealizabilitySettings::cqmom, source_name);
    std::optional<CqmomDiagnosticLevel> diagnostic_level = parse_cqmom_diagnostic_level(parsed, source_name);
    std::optional<VandermondeSolver> vandermonde_solver = parse_vandermonde_solver(parsed, source_name);
    std::optional<double> moore_penrose_rcond = parse_vandermonde_rcond(parsed, source_name);
    std::optional<QmomOptions> qmom_options = parse_qmom_options_from_parsed(parsed, source_name);
    validate_format_version(parsed, source_name);
    const Entry* node_counts_entry = require_entry(parsed, "phase_space.node_counts", source_name);

    std::vector<std::size_t> node_counts;
    if (node_counts_entry != nullptr)
    {
        std::optional<std::vector<std::size_t>> parsed_node_counts =
            parse_phase_space_node_counts(parsed, *node_counts_entry, source_name, "node_counts must contain positive integers.");
        if (parsed_node_counts.has_value())
        {
            node_counts = std::move(*parsed_node_counts);
            validate_cqmom_phase_space(parsed, source_name, *node_counts_entry, node_counts);
        }
    }

    std::vector<std::size_t> conditioning_order(node_counts.size());
    // Omitted conditioning_order means identity order: dimensions are
    // conditioned in their declared phase-space order.
    for (std::size_t dimension = 0; dimension < conditioning_order.size(); ++dimension)
    {
        conditioning_order[dimension] = dimension;
    }
    const auto conditioning_entry = parsed.scalars.find("phase_space.conditioning_order");
    if (conditioning_entry != parsed.scalars.end())
    {
        conditioning_order.clear();
        if (!parse_size_list(conditioning_entry->second, conditioning_order))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, conditioning_entry->second.line,
                      "conditioning_order must contain integer dimension indices.");
        }
    }
    std::vector<bool> seen(node_counts.size(), false);
    bool permutation = conditioning_order.size() == node_counts.size();
    for (const std::size_t dimension : conditioning_order)
    {
        if (dimension >= seen.size() || seen[dimension])
        {
            permutation = false;
            break;
        }
        seen[dimension] = true;
    }
    if (!permutation)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name,
                  conditioning_entry == parsed.scalars.end() ? 0 : conditioning_entry->second.line,
                  "conditioning_order must be a permutation of the phase-space dimensions.");
    }
    std::optional<CqmomAbscissaBounds> coordinate_support = parse_coordinate_support(parsed, source_name, node_counts.size());

    if (!node_counts.empty())
    {
        validate_cqmom_coordinate_realizability_ranks(parsed, source_name, node_counts.size());
    }
    if (!parsed.errors.empty())
    {
        return failure_result<ParsedCqmomCaseSettings>(std::move(parsed));
    }

    try
    {
        ParsedCqmomCaseSettings settings{PhaseSpaceHierarchy{std::move(node_counts), std::move(conditioning_order)},
                                         qmom_options.value_or(QmomOptions{}),
                                         diagnostic_level.value_or(CqmomDiagnosticLevel::none),
                                         vandermonde_solver.value_or(VandermondeSolver::automatic),
                                         moore_penrose_rcond.value_or(1.0e-14),
                                         coordinate_support.value_or(CqmomAbscissaBounds{})};
        return ConfigurationResult<ParsedCqmomCaseSettings>{std::move(settings)};
    }
    catch (const std::exception& error)
    {
        return failure_result<ParsedCqmomCaseSettings>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0, error.what()});
    }
}

// Parses and validates the output settings. A file without output settings
// yields the defaults and no errors.
OutputOptions parse_output_options_from_parsed(ParsedFile& parsed, const std::string& source_name)
{
    OutputOptions options;
    if (!has_output_section(parsed))
    {
        return options;
    }

    static const std::vector<std::string> allowed = {"output.directory", "output.write_moments", "output.write_quadrature", "output.interval",
                                                     "output.generalized_moments"};
    for (const auto& item : parsed.scalars)
    {
        if (!is_output_key(item.first))
        {
            continue;
        }
        if (std::find(allowed.begin(), allowed.end(), item.first) == allowed.end())
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }
    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("output.", 0) == 0)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }

    const Entry* directory = find_entry(parsed, "output.directory");
    if (directory != nullptr)
    {
        const std::string cleaned_directory = trim(directory->value);
        if (cleaned_directory.empty())
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, directory->line, "output.directory must be non-empty.");
        }
        else
        {
            options.directory = cleaned_directory;
        }
    }

    parse_optional_bool(parsed, "output.write_moments", source_name, options.write_moments);
    parse_optional_bool(parsed, "output.write_quadrature", source_name, options.write_quadrature);

    const Entry* interval = find_entry(parsed, "output.interval");
    if (interval != nullptr)
    {
        double value = 0.0;
        if (!parse_double(interval->value, value))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, interval->line,
                      "output.interval must be a finite positive number.");
        }
        else if (value <= 0.0)
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, interval->line, "output.interval must be positive.");
        }
        else
        {
            options.interval = value;
        }
    }

    const Entry* generalized_moments = find_entry(parsed, "output.generalized_moments");
    if (generalized_moments != nullptr && generalized_moments->value.back() == '|')
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, generalized_moments->line,
                  "output.generalized_moments must contain non-negative finite comma-separated exponent vectors separated by '|'.");
    }
    else if (generalized_moments != nullptr)
    {
        std::istringstream requests(generalized_moments->value);
        std::string request;
        while (std::getline(requests, request, '|'))
        {
            Entry request_entry{trim(request), generalized_moments->line};
            std::vector<double> orders;
            if (!parse_double_list(request_entry, orders) ||
                std::any_of(orders.begin(), orders.end(), [](const double order) { return order < 0.0; }))
            {
                add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, generalized_moments->line,
                          "output.generalized_moments must contain non-negative finite comma-separated exponent vectors separated by '|'.");
                break;
            }
            options.generalized_moment_orders.push_back(std::move(orders));
        }
        if (options.generalized_moment_orders.empty())
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, generalized_moments->line,
                      "output.generalized_moments must contain at least one exponent vector.");
        }
    }

    return options;
}

ConfigurationResult<InitialConditionSpecification> parse_qmom_initial_condition_from_parsed(ParsedFile parsed, std::string source_name)
{
    const bool explicit_moments = has_explicit_qmom_moments(parsed);
    reject_unknown_qmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_qmom_keys(parsed, explicit_moments, source_name);
    validate_optional_phase_space_evolution_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, RealizabilitySettings::qmom, source_name);
    validate_format_version(parsed, source_name);
    static_cast<void>(parse_output_options_from_parsed(parsed, source_name));
    validate_optional_qmom_phase_space(parsed, source_name);
    if (!parsed.moments.empty())
    {
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, parsed.moments.begin()->second.line,
                  "Indexed moments are not supported for QMOM configurations.");
    }

    const Entry* moments_entry = find_entry(parsed, "moments.values");
    if (moments_entry != nullptr)
    {
        std::vector<double> moments;
        if (!parse_double_list(*moments_entry, moments))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, moments_entry->line,
                      "QMOM moments must be a comma-separated list of finite numbers.");
        }
        if (!parsed.errors.empty())
        {
            return failure_result<InitialConditionSpecification>(std::move(parsed));
        }
        return ConfigurationResult<InitialConditionSpecification>{
            InitialConditionSpecification{ExplicitMomentInitialization{MomentSet(std::move(moments))}}};
    }

    const Entry* type_entry = find_entry(parsed, "initial.type");
    if (type_entry == nullptr)
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "QMOM configuration requires either '[moments] values', '[initial] type = distribution', or '[initial] type = mixture'.");
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }
    if (type_entry->value != "distribution" && type_entry->value != "mixture")
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, type_entry->line,
                  "Only '[initial] type = distribution' and '[initial] type = mixture' are supported for QMOM initialization.");
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }

    std::size_t maximum_order = 0;
    if (find_entry(parsed, "initial.maximum_order") != nullptr)
    {
        parse_optional_size(parsed, "initial.maximum_order", source_name, maximum_order);
    }
    else if (parsed.errors.empty())
    {
        const std::optional<std::size_t> inferred_maximum_order = infer_qmom_initial_maximum_order(parsed, source_name);
        if (inferred_maximum_order.has_value())
        {
            maximum_order = *inferred_maximum_order;
        }
    }

    if (type_entry->value == "distribution")
    {
        if (has_mixture_component_sections(parsed))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, type_entry->line,
                      "Mixture component sections are only valid when '[initial] type = mixture'.");
        }
        std::optional<DistributionSpecification> distribution = parse_distribution(parsed, source_name);
        reject_unconsumed_indexed_entries(parsed, std::string(mixture_component_section_prefix), source_name);
        double number_density = 1.0;
        parse_optional_double(parsed, "initial.number_density", source_name, number_density);
        if (!std::isfinite(number_density) || number_density < 0.0)
        {
            const Entry* entry = find_entry(parsed, "initial.number_density");
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? 0 : entry->line,
                      "initial.number_density must be finite and non-negative.");
        }

        if (!parsed.errors.empty() || !distribution.has_value())
        {
            return failure_result<InitialConditionSpecification>(std::move(parsed));
        }

        return ConfigurationResult<InitialConditionSpecification>{
            InitialConditionSpecification{DistributionInitialization{*distribution, maximum_order, number_density}}};
    }

    std::optional<DistributionMixture> mixture = parse_mixture(parsed, source_name);
    reject_unconsumed_indexed_entries(parsed, std::string(mixture_component_section_prefix), source_name);
    if (!parsed.errors.empty() || !mixture.has_value())
    {
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }
    return ConfigurationResult<InitialConditionSpecification>{InitialConditionSpecification{MixtureInitialization{*mixture, maximum_order}}};
}

ConfigurationResult<InitialConditionSpecification> parse_cqmom_initial_condition_from_parsed(ParsedFile parsed, std::string source_name)
{
    const bool explicit_moments = has_explicit_cqmom_moments(parsed);
    reject_unknown_cqmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_cqmom_keys(parsed, explicit_moments, source_name);
    validate_optional_phase_space_evolution_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, RealizabilitySettings::cqmom, source_name);
    validate_optional_cqmom_diagnostic_level(parsed, source_name);
    validate_optional_vandermonde_solver(parsed, source_name);
    validate_format_version(parsed, source_name);
    static_cast<void>(parse_output_options_from_parsed(parsed, source_name));
    validate_optional_cqmom_phase_space(parsed, source_name);

    if (explicit_moments)
    {
        std::optional<RawMomentTensor> tensor = parse_explicit_tensor(parsed, source_name);
        if (!parsed.errors.empty() || !tensor.has_value())
        {
            return failure_result<InitialConditionSpecification>(std::move(parsed));
        }
        return ConfigurationResult<InitialConditionSpecification>{InitialConditionSpecification{ExplicitTensorInitialization{*tensor}}};
    }

    const Entry* type_entry = find_entry(parsed, "initial.type");
    if (type_entry == nullptr)
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "CQMOM configuration requires either explicit '[moments]' tensor entries or '[initial] type = independent_product'.");
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }
    if (type_entry->value != "independent_product")
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, type_entry->line,
                  "Only '[initial] type = independent_product' is supported for CQMOM analytical initialization.");
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }

    std::optional<std::vector<std::size_t>> maximum_orders;
    const Entry* maximum_orders_entry = find_entry(parsed, "initial.maximum_orders");
    if (maximum_orders_entry != nullptr)
    {
        std::vector<std::size_t> parsed_maximum_orders;
        if (!parse_size_list(*maximum_orders_entry, parsed_maximum_orders))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, maximum_orders_entry->line,
                      "initial.maximum_orders must contain non-negative integers.");
        }
        else
        {
            maximum_orders = std::move(parsed_maximum_orders);
        }
    }
    if (maximum_orders.has_value())
    {
        for (const std::size_t maximum_order : *maximum_orders)
        {
            if (maximum_order == std::numeric_limits<std::size_t>::max())
            {
                add_error(parsed, ConfigurationErrorCode::invalid_value, source_name,
                          maximum_orders_entry == nullptr ? 0 : maximum_orders_entry->line,
                          "initial.maximum_orders contains an order that is too large.");
                break;
            }
        }
    }
    if (maximum_orders.has_value())
    {
        const Entry* node_counts_entry = find_entry(parsed, "phase_space.node_counts");
        std::vector<std::size_t> node_counts;
        if (node_counts_entry != nullptr && parse_size_list(*node_counts_entry, node_counts) && maximum_orders->size() != node_counts.size())
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name,
                      maximum_orders_entry == nullptr ? 0 : maximum_orders_entry->line,
                      "initial.maximum_orders must contain one entry per phase-space dimension.");
        }

        std::vector<std::size_t> extents;
        extents.reserve(maximum_orders->size());
        for (const std::size_t maximum_order : *maximum_orders)
        {
            if (maximum_order != std::numeric_limits<std::size_t>::max())
            {
                extents.push_back(maximum_order + 1);
            }
        }
        std::size_t tensor_size = 0;
        if (extents.size() == maximum_orders->size() && !checked_product(extents, tensor_size))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name,
                      maximum_orders_entry == nullptr ? 0 : maximum_orders_entry->line,
                      "The full tensor of raw moments given by initial.maximum_orders is too large: the product of its extents overflows.");
        }
    }

    std::optional<IndependentProductDistribution> distribution = parse_independent_product(parsed, source_name);
    reject_unconsumed_indexed_entries(parsed, std::string(independent_product_dimension_section_prefix), source_name);
    const Entry* phase_space_node_counts = find_entry(parsed, "phase_space.node_counts");
    std::vector<std::size_t> declared_node_counts;
    if (distribution.has_value() && phase_space_node_counts != nullptr && parse_size_list(*phase_space_node_counts, declared_node_counts) &&
        distribution->marginals.size() != declared_node_counts.size())
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, phase_space_node_counts->line,
                  "Independent product initialization must define one '[initial.dimension.N]' section per phase-space dimension.");
    }
    if (distribution.has_value() && maximum_orders.has_value() && distribution->marginals.size() != maximum_orders->size())
    {
        add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name,
                  maximum_orders_entry == nullptr ? 0 : maximum_orders_entry->line,
                  "initial.maximum_orders must contain one entry per independent-product dimension.");
    }

    if (!parsed.errors.empty() || !distribution.has_value())
    {
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }

    return ConfigurationResult<InitialConditionSpecification>{
        InitialConditionSpecification{IndependentProductInitialization{*distribution, std::move(maximum_orders)}}};
}

} // namespace

ConfigurationResult<MomentSet> initialize_qmom_moments(const InitialConditionSpecification& specification, std::string source_name)
{
    if (const auto* explicit_moments = std::get_if<ExplicitMomentInitialization>(&specification))
    {
        return ConfigurationResult<MomentSet>{explicit_moments->moments};
    }
    if (std::holds_alternative<ExplicitTensorInitialization>(specification) ||
        std::holds_alternative<IndependentProductInitialization>(specification))
    {
        return failure_result<MomentSet>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                                          "One-dimensional initialization cannot use a multidimensional tensor specification."});
    }

    RawMomentResult moments;
    try
    {
        if (const auto* distribution = std::get_if<DistributionInitialization>(&specification))
        {
            moments = compute_raw_moments(distribution->distribution, distribution->maximum_order, distribution->number_density);
        }
        else
        {
            const auto& mixture = std::get<MixtureInitialization>(specification);
            moments = compute_raw_moments(mixture.mixture, mixture.maximum_order);
        }
    }
    catch (const std::exception& error)
    {
        // A maximum order far beyond any node count can exhaust memory.
        return failure_result<MomentSet>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                                          std::string("Analytical initialization moment generation failed: ") + error.what()});
    }
    if (!moments.ok())
    {
        std::string message = "Analytical initialization moment generation failed: " + moments.message;
        if (moments.failed_component != invalid_component_index)
        {
            message += " Component index: " + std::to_string(moments.failed_component) + ".";
        }
        return failure_result<MomentSet>({ConfigurationErrorCode::invalid_value, std::move(source_name), 0, std::move(message)});
    }
    return ConfigurationResult<MomentSet>{moments.moments};
}

ConfigurationResult<RawMomentTensor> initialize_cqmom_moment_tensor(const InitialConditionSpecification& specification, std::string source_name)
{
    if (const auto* explicit_tensor = std::get_if<ExplicitTensorInitialization>(&specification))
    {
        return ConfigurationResult<RawMomentTensor>{explicit_tensor->moments};
    }
    if (const auto* independent_product = std::get_if<IndependentProductInitialization>(&specification))
    {
        if (!independent_product->maximum_orders.has_value())
        {
            return failure_result<RawMomentTensor>({ConfigurationErrorCode::missing_entry, std::move(source_name), 0,
                                                    "initial.maximum_orders is required to initialize the full tensor of raw moments."});
        }
        TensorMomentResult moments;
        try
        {
            moments = compute_raw_moments(independent_product->distribution, *independent_product->maximum_orders);
        }
        catch (const std::exception& error)
        {
            // Maximum orders far beyond any node count can exhaust memory.
            return failure_result<RawMomentTensor>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                                                    std::string("Independent product initialization moment generation failed: ") + error.what()});
        }
        if (!moments.ok())
        {
            std::string message = "Independent product initialization moment generation failed: " + moments.message;
            if (moments.failed_dimension != invalid_dimension_index)
            {
                message += " Dimension index: " + std::to_string(moments.failed_dimension) + ".";
            }
            return failure_result<RawMomentTensor>({ConfigurationErrorCode::invalid_value, std::move(source_name), 0, std::move(message)});
        }
        return ConfigurationResult<RawMomentTensor>{moments.moments};
    }
    return failure_result<RawMomentTensor>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                                            "CQMOM initialization requires explicit tensor moments or an independent product specification."});
}

ConfigurationResult<CqmomMomentData> initialize_cqmom_moments(const InitialConditionSpecification& specification,
                                                              const PhaseSpaceHierarchy& hierarchy, std::string source_name)
{
    try
    {
        const CqmomMomentLayout layout = create_cqmom_moment_layout(hierarchy);
        if (const auto* explicit_tensor = std::get_if<ExplicitTensorInitialization>(&specification))
        {
            return ConfigurationResult<CqmomMomentData>{extract_cqmom_moments(explicit_tensor->moments, layout)};
        }
        if (const auto* independent_product = std::get_if<IndependentProductInitialization>(&specification))
        {
            // Requesting only the recursive CQMOM moments required by the layout
            // avoids computing tensor entries the case never reads.
            const RawMomentRequest request{layout.required_indices.to_multi_indices()};
            const SelectedRawMomentResult selected_moments = compute_raw_moments(independent_product->distribution, request);
            if (!selected_moments.ok())
            {
                std::string message = "Independent product selected moment generation failed: " + selected_moments.message;
                if (selected_moments.failed_dimension != invalid_dimension_index)
                {
                    message += " Dimension index: " + std::to_string(selected_moments.failed_dimension) + ".";
                }
                if (!selected_moments.failed_order.empty())
                {
                    message += " Moment index:";
                    for (std::size_t component = 0; component < selected_moments.failed_order.size(); ++component)
                    {
                        message += component == 0 ? " [" : ", ";
                        message += std::to_string(selected_moments.failed_order[component]);
                    }
                    message += "].";
                }
                return failure_result<CqmomMomentData>({ConfigurationErrorCode::invalid_value, std::move(source_name), 0, std::move(message)});
            }
            return ConfigurationResult<CqmomMomentData>{create_cqmom_moment_data(layout, selected_moments)};
        }
    }
    catch (const std::exception& error)
    {
        return failure_result<CqmomMomentData>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0, error.what()});
    }

    return failure_result<CqmomMomentData>({ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                                            "CQMOM initialization requires explicit tensor moments or an independent product specification."});
}

ConfigurationResult<InitialConditionSpecification> parse_qmom_initial_condition(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    return parse_qmom_initial_condition_from_parsed(std::move(parsed), std::move(source_name));
}

ConfigurationResult<InitialConditionSpecification> parse_cqmom_initial_condition(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    return parse_cqmom_initial_condition_from_parsed(std::move(parsed), std::move(source_name));
}

ConfigurationResult<QmomCase> create_qmom_case(const InitialConditionSpecification& initial_condition, const std::size_t node_count,
                                               const QmomOptions& options, std::string source_name, std::optional<AbscissaBounds> abscissa_bounds)
{
    if (node_count == 0)
    {
        return failure_result<QmomCase>({ConfigurationErrorCode::invalid_value, std::move(source_name), 0, "QMOM node_count must be positive."});
    }
    const std::optional<std::size_t> required_moment_count = internal::checked_size_multiply(2, node_count);
    if (!required_moment_count.has_value())
    {
        return failure_result<QmomCase>(
            {ConfigurationErrorCode::invalid_value, std::move(source_name), 0, "QMOM node_count is too large to derive the required moment count."});
    }

    ConfigurationResult<MomentSet> initialized = initialize_qmom_moments(initial_condition, source_name);
    if (!initialized)
    {
        return ConfigurationResult<QmomCase>::failure(initialized.error());
    }
    if (initialized.value().size() < *required_moment_count)
    {
        return failure_result<QmomCase>(
            {ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0, "QMOM requires at least 2 * node_count moments."});
    }
    return ConfigurationResult<QmomCase>{QmomCase{std::move(initialized).value(), node_count, options, std::move(abscissa_bounds)}};
}

ConfigurationResult<CqmomCase> create_cqmom_case(const InitialConditionSpecification& initial_condition, const PhaseSpaceHierarchy& hierarchy,
                                                 const QmomOptions& qmom_options, std::string source_name, CqmomAbscissaBounds coordinate_support)
{
    ConfigurationResult<CqmomMomentData> initialized = initialize_cqmom_moments(initial_condition, hierarchy, std::move(source_name));
    if (!initialized)
    {
        return ConfigurationResult<CqmomCase>::failure(initialized.error());
    }
    return ConfigurationResult<CqmomCase>{CqmomCase{
        std::move(initialized).value(), {}, qmom_options, CqmomDiagnosticLevel::none, VandermondeSolver::automatic, std::move(coordinate_support)}};
}

ConfigurationResult<CqmomCase> with_initial_moments(const CqmomCase& base_cqmom_case, RawMomentTensor explicit_moments)
{
    try
    {
        const PhaseSpaceHierarchy& hierarchy = base_cqmom_case.moments.layout.hierarchy;
        const std::optional<std::string> structure_error = validate_replacement_tensor_structure(explicit_moments, hierarchy);
        if (structure_error.has_value())
        {
            return ConfigurationResult<CqmomCase>::failure(*structure_error);
        }

        const CqmomMomentLayout layout = create_cqmom_moment_layout(hierarchy);
        for (const MultiIndexView index : layout.required_indices)
        {
            if (!internal::multi_index_within_extents(index, explicit_moments.extents))
            {
                return ConfigurationResult<CqmomCase>::failure("The full tensor of raw moments does not contain required CQMOM moment index " +
                                                               format_moment_index(index) + ".");
            }
        }

        return ConfigurationResult<CqmomCase>{qmomkit::internal::with_moments(base_cqmom_case, extract_cqmom_moments(explicit_moments, layout))};
    }
    catch (const std::exception& error)
    {
        return ConfigurationResult<CqmomCase>::failure("Unable to replace CQMOM initial moments: " + std::string(error.what()));
    }
}

ConfigurationResult<CqmomCase> with_initial_moments(const CqmomCase& base_cqmom_case, CqmomMomentData explicit_moments)
{
    try
    {
        internal::validate_cqmom_moment_data_matches_hierarchy(
            explicit_moments, base_cqmom_case.moments.layout.hierarchy,
            {"CQMOM moment data layout must match the required CQMOM moment layout of the base hierarchy.",
             "CQMOM moment data value count does not match its required CQMOM moment layout.", "CQMOM moment data values must be finite."});

        return ConfigurationResult<CqmomCase>{qmomkit::internal::with_moments(base_cqmom_case, std::move(explicit_moments))};
    }
    catch (const std::exception& error)
    {
        return ConfigurationResult<CqmomCase>::failure("Unable to replace CQMOM initial moments: " + std::string(error.what()));
    }
}

namespace
{
// The prepared overloads validate through the CqmomCase overloads, so both
// report the same errors. The replaced case keeps the base hierarchy, so the
// base case's prepared metadata stays valid for it and is shared.
ConfigurationResult<CqmomPreparedCase> with_prepared_metadata_of(ConfigurationResult<CqmomCase> replaced_cqmom_case_result,
                                                                 const CqmomPreparedCase& base_cqmom_prepared_case)
{
    if (!replaced_cqmom_case_result)
    {
        return ConfigurationResult<CqmomPreparedCase>::failure(replaced_cqmom_case_result.error());
    }
    try
    {
        return ConfigurationResult<CqmomPreparedCase>{qmomkit::internal::CqmomConstructionAccess::create_prepared_case(
            std::move(replaced_cqmom_case_result).value(),
            qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(base_cqmom_prepared_case))};
    }
    catch (const std::exception& error)
    {
        return ConfigurationResult<CqmomPreparedCase>::failure("Unable to replace CQMOM initial moments: " + std::string(error.what()));
    }
}
} // namespace

ConfigurationResult<CqmomPreparedCase> with_initial_moments(const CqmomPreparedCase& base_cqmom_prepared_case, RawMomentTensor explicit_moments)
{
    return with_prepared_metadata_of(
        with_initial_moments(qmomkit::internal::CqmomConstructionAccess::case_definition(base_cqmom_prepared_case), std::move(explicit_moments)),
        base_cqmom_prepared_case);
}

ConfigurationResult<CqmomPreparedCase> with_initial_moments(const CqmomPreparedCase& base_cqmom_prepared_case, CqmomMomentData explicit_moments)
{
    return with_prepared_metadata_of(
        with_initial_moments(qmomkit::internal::CqmomConstructionAccess::case_definition(base_cqmom_prepared_case), std::move(explicit_moments)),
        base_cqmom_prepared_case);
}

ConfigurationResult<QmomCase> parse_qmom(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    const bool explicit_moments = has_explicit_qmom_moments(parsed);
    reject_unknown_qmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_qmom_keys(parsed, explicit_moments, source_name);
    validate_optional_phase_space_evolution_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, RealizabilitySettings::qmom, source_name);
    validate_format_version(parsed, source_name);
    std::optional<QmomOptions> qmom_options = parse_qmom_options_from_parsed(parsed, source_name);
    const Entry* node_counts_entry = require_entry(parsed, "phase_space.node_counts", source_name);
    if (!parsed.moments.empty())
    {
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, parsed.moments.begin()->second.line,
                  "Indexed moments are not supported for QMOM configurations.");
    }

    std::vector<std::size_t> node_counts;
    if (node_counts_entry != nullptr)
    {
        std::optional<std::vector<std::size_t>> parsed_node_counts =
            parse_phase_space_node_counts(parsed, *node_counts_entry, source_name, "QMOM node_counts must contain exactly one positive integer.");
        if (parsed_node_counts.has_value())
        {
            node_counts = std::move(*parsed_node_counts);
            validate_qmom_phase_space(parsed, source_name, *node_counts_entry, node_counts);
        }
    }
    if (!parsed.errors.empty())
    {
        return failure_result<QmomCase>(std::move(parsed));
    }

    std::optional<CqmomAbscissaBounds> coordinate_support = parse_coordinate_support(parsed, source_name, node_counts.size());
    if (!coordinate_support.has_value() || !parsed.errors.empty())
    {
        return failure_result<QmomCase>(std::move(parsed));
    }

    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_qmom_initial_condition_from_parsed(std::move(parsed), source_name);
    if (!initial_condition)
    {
        return ConfigurationResult<QmomCase>::failure(initial_condition.error());
    }
    const CqmomAbscissaBounds& parsed_support = coordinate_support.value();
    const std::optional<AbscissaBounds> abscissa_bounds = parsed_support.empty() ? std::nullopt : parsed_support.front();
    return create_qmom_case(initial_condition.value(), node_counts[0], qmom_options.value_or(QmomOptions{}), std::move(source_name), abscissa_bounds);
}

ConfigurationResult<CqmomCase> parse_cqmom(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    ConfigurationResult<ParsedCqmomCaseSettings> settings = parse_cqmom_case_settings(parsed, source_name);
    if (!settings)
    {
        return ConfigurationResult<CqmomCase>::failure(settings.error());
    }
    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_cqmom_initial_condition_from_parsed(std::move(parsed), source_name);
    if (!initial_condition)
    {
        return ConfigurationResult<CqmomCase>::failure(initial_condition.error());
    }
    ConfigurationResult<CqmomCase> cqmom_case_result =
        create_cqmom_case(initial_condition.value(), settings.value().hierarchy, settings.value().qmom_options, std::move(source_name),
                          settings.value().coordinate_support);
    if (cqmom_case_result)
    {
        cqmom_case_result.value().diagnostic_level = settings.value().diagnostic_level;
        cqmom_case_result.value().vandermonde_solver = settings.value().vandermonde_solver;
        cqmom_case_result.value().conditional_options.moore_penrose_rcond = settings.value().moore_penrose_rcond;
    }
    return cqmom_case_result;
}

ConfigurationResult<CqmomPreparedCase> parse_cqmom_prepared(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    ConfigurationResult<ParsedCqmomCaseSettings> settings = parse_cqmom_case_settings(parsed, source_name);
    if (!settings)
    {
        return ConfigurationResult<CqmomPreparedCase>::failure(settings.error());
    }

    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_cqmom_initial_condition_from_parsed(std::move(parsed), source_name);
    if (!initial_condition)
    {
        return ConfigurationResult<CqmomPreparedCase>::failure(initial_condition.error());
    }

    ConfigurationResult<CqmomCase> cqmom_case_result = create_cqmom_case(
        initial_condition.value(), settings.value().hierarchy, settings.value().qmom_options, source_name, settings.value().coordinate_support);
    if (!cqmom_case_result)
    {
        return ConfigurationResult<CqmomPreparedCase>::failure(cqmom_case_result.error());
    }
    cqmom_case_result.value().diagnostic_level = settings.value().diagnostic_level;
    cqmom_case_result.value().vandermonde_solver = settings.value().vandermonde_solver;
    cqmom_case_result.value().conditional_options.moore_penrose_rcond = settings.value().moore_penrose_rcond;
    return ConfigurationResult<CqmomPreparedCase>{
        qmomkit::internal::CqmomConstructionAccess::create_prepared_case(std::move(cqmom_case_result).value())};
}

ConfigurationResult<RealizabilityPolicy> parse_realizability_policy(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<RealizabilityPolicy>(std::move(parsed));
    }
    if (!has_realizability_section(parsed))
    {
        return ConfigurationResult<RealizabilityPolicy>{RealizabilityPolicy{}};
    }

    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("time_integration.", 0) == 0 && item.first != realizability_section && item.first != node_advection_section &&
            item.first != adaptive_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }

    // Per-coordinate settings belong to CQMOM; this parser reads the QMOM policy only.
    for (const auto& item : parsed.scalars)
    {
        if (is_coordinate_realizability_key(item.first))
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }

    const std::optional<RealizabilityPolicy> policy = parse_realizability_policy_from_parsed(parsed, std::string(realizability_section), source_name);
    if (!parsed.errors.empty() || !policy.has_value())
    {
        return failure_result<RealizabilityPolicy>(std::move(parsed));
    }
    return ConfigurationResult<RealizabilityPolicy>{*policy};
}

ConfigurationResult<std::vector<std::optional<RealizabilityPolicy>>> parse_coordinate_realizability_policies(const std::string_view text,
                                                                                                             std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<std::vector<std::optional<RealizabilityPolicy>>>(std::move(parsed));
    }

    reject_unknown_phase_space_evolution_sections(parsed, RealizabilitySettings::cqmom, source_name);
    reject_unknown_phase_space_evolution_keys(parsed, RealizabilitySettings::cqmom, source_name);
    std::vector<std::optional<RealizabilityPolicy>> policies = parse_coordinate_realizability_policies_from_parsed(parsed, source_name);

    if (!parsed.errors.empty())
    {
        return failure_result<std::vector<std::optional<RealizabilityPolicy>>>(std::move(parsed));
    }
    return ConfigurationResult<std::vector<std::optional<RealizabilityPolicy>>>{std::move(policies)};
}

ConfigurationResult<PhaseSpaceEvolutionScheme> parse_phase_space_evolution_scheme(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<PhaseSpaceEvolutionScheme>(std::move(parsed));
    }

    reject_unknown_phase_space_evolution_sections(parsed, RealizabilitySettings::any, source_name);
    reject_unknown_phase_space_evolution_keys(parsed, RealizabilitySettings::any, source_name);
    const std::optional<PhaseSpaceEvolutionScheme> scheme = parse_phase_space_evolution_scheme_from_parsed(parsed, source_name);

    if (!parsed.errors.empty() || !scheme.has_value())
    {
        return failure_result<PhaseSpaceEvolutionScheme>(std::move(parsed));
    }
    return ConfigurationResult<PhaseSpaceEvolutionScheme>{*scheme};
}

ConfigurationResult<NodeAdvectionOptions> parse_node_advection_options(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<NodeAdvectionOptions>(std::move(parsed));
    }

    reject_unknown_phase_space_evolution_sections(parsed, RealizabilitySettings::any, source_name);
    reject_unknown_phase_space_evolution_keys(parsed, RealizabilitySettings::any, source_name);
    const std::optional<NodeAdvectionOptions> options = parse_node_advection_options_from_parsed(parsed, source_name);

    if (!parsed.errors.empty() || !options.has_value())
    {
        return failure_result<NodeAdvectionOptions>(std::move(parsed));
    }
    return ConfigurationResult<NodeAdvectionOptions>{*options};
}

ConfigurationResult<QmomPhaseSpaceEvolutionOptions> parse_qmom_phase_space_evolution_options(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<QmomPhaseSpaceEvolutionOptions>(std::move(parsed));
    }

    reject_unknown_phase_space_evolution_sections(parsed, RealizabilitySettings::qmom, source_name);
    reject_unknown_phase_space_evolution_keys(parsed, RealizabilitySettings::qmom, source_name);
    reject_unknown_qmom_option_keys(parsed, source_name);
    const std::optional<PhaseSpaceEvolutionScheme> scheme = parse_phase_space_evolution_scheme_from_parsed(parsed, source_name);
    const std::optional<NodeAdvectionOptions> node_options = parse_node_advection_options_from_parsed(parsed, source_name);
    (void)parse_qmom_options_from_parsed(parsed, source_name);
    const std::optional<RealizabilityPolicy> realizability_policy =
        parse_realizability_policy_from_parsed(parsed, std::string(realizability_section), source_name);

    if (!parsed.errors.empty() || !scheme.has_value() || !node_options.has_value() || !realizability_policy.has_value())
    {
        return failure_result<QmomPhaseSpaceEvolutionOptions>(std::move(parsed));
    }

    QmomPhaseSpaceEvolutionOptions options;
    options.scheme = *scheme;
    options.node_options = *node_options;
    options.realizability_policy = *realizability_policy;
    return ConfigurationResult<QmomPhaseSpaceEvolutionOptions>{std::move(options)};
}

ConfigurationResult<CqmomPhaseSpaceEvolutionOptions> parse_cqmom_phase_space_evolution_options(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<CqmomPhaseSpaceEvolutionOptions>(std::move(parsed));
    }

    reject_unknown_phase_space_evolution_sections(parsed, RealizabilitySettings::cqmom, source_name);
    reject_unknown_phase_space_evolution_keys(parsed, RealizabilitySettings::cqmom, source_name);
    const std::optional<PhaseSpaceEvolutionScheme> scheme = parse_phase_space_evolution_scheme_from_parsed(parsed, source_name);
    const std::optional<NodeAdvectionOptions> node_options = parse_node_advection_options_from_parsed(parsed, source_name);
    std::vector<std::optional<RealizabilityPolicy>> coordinate_policies = parse_coordinate_realizability_policies_from_parsed(parsed, source_name);

    if (!parsed.errors.empty() || !scheme.has_value() || !node_options.has_value())
    {
        return failure_result<CqmomPhaseSpaceEvolutionOptions>(std::move(parsed));
    }

    CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = *scheme;
    options.node_options = *node_options;
    options.source_options.coordinate_realizability_policies = std::move(coordinate_policies);
    return ConfigurationResult<CqmomPhaseSpaceEvolutionOptions>{std::move(options)};
}

ConfigurationResult<OutputOptions> parse_output_options(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<OutputOptions>(std::move(parsed));
    }
    OutputOptions options = parse_output_options_from_parsed(parsed, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<OutputOptions>(std::move(parsed));
    }
    return ConfigurationResult<OutputOptions>{options};
}

namespace
{

template <typename T, typename ParseFunction>
[[nodiscard]] ConfigurationResult<T> load_configuration(const std::filesystem::path& filename, ParseFunction parse)
{
    std::optional<std::string> contents = read_configuration_file(filename);
    if (!contents.has_value())
    {
        return file_failure<T>(filename);
    }
    return parse(*contents, filename.string());
}

} // namespace

ConfigurationResult<QmomCase> load_qmom(const std::filesystem::path& filename)
{
    return load_configuration<QmomCase>(filename, parse_qmom);
}

ConfigurationResult<CqmomCase> load_cqmom(const std::filesystem::path& filename)
{
    return load_configuration<CqmomCase>(filename, parse_cqmom);
}

ConfigurationResult<CqmomPreparedCase> load_cqmom_prepared(const std::filesystem::path& filename)
{
    return load_configuration<CqmomPreparedCase>(filename, parse_cqmom_prepared);
}

ConfigurationResult<OutputOptions> load_output_options(const std::filesystem::path& filename)
{
    return load_configuration<OutputOptions>(filename, parse_output_options);
}

ConfigurationResult<NodeAdvectionOptions> load_node_advection_options(const std::filesystem::path& filename)
{
    return load_configuration<NodeAdvectionOptions>(filename, parse_node_advection_options);
}

ConfigurationResult<RealizabilityPolicy> load_realizability_policy(const std::filesystem::path& filename)
{
    return load_configuration<RealizabilityPolicy>(filename, parse_realizability_policy);
}

ConfigurationResult<QmomPhaseSpaceEvolutionOptions> load_qmom_phase_space_evolution_options(const std::filesystem::path& filename)
{
    return load_configuration<QmomPhaseSpaceEvolutionOptions>(filename, parse_qmom_phase_space_evolution_options);
}

ConfigurationResult<CqmomPhaseSpaceEvolutionOptions> load_cqmom_phase_space_evolution_options(const std::filesystem::path& filename)
{
    return load_configuration<CqmomPhaseSpaceEvolutionOptions>(filename, parse_cqmom_phase_space_evolution_options);
}

} // namespace qmomkit::config
