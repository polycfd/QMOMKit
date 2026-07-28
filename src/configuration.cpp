#include "qmomkit/configuration.hpp"

#include "internal/cqmom_case.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
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
};

constexpr std::string_view mixture_component_section_prefix = "initial.mixture.component.";
constexpr std::string_view independent_product_dimension_section_prefix = "initial.dimension.";
constexpr std::string_view realizability_section = "time_integration.realizability";
constexpr std::string_view realizability_key_prefix = "time_integration.realizability.";
constexpr std::string_view transport_scheme_key = "time_integration.scheme";
constexpr std::string_view node_advection_section = "time_integration.node_advection";
constexpr std::string_view node_advection_key_prefix = "time_integration.node_advection.";
constexpr std::string_view node_advection_scheme_key = "time_integration.node_advection.scheme";
constexpr std::string_view rk23_section = "time_integration.node_advection.rk23";
constexpr std::string_view rk23_relative_tolerance_key = "time_integration.node_advection.rk23.relative_tolerance";
constexpr std::string_view rk23_absolute_tolerance_key = "time_integration.node_advection.rk23.absolute_tolerance";
constexpr std::string_view rk23_initial_step_key = "time_integration.node_advection.rk23.initial_step";
constexpr std::string_view rk23_max_step_key = "time_integration.node_advection.rk23.max_step";
constexpr std::string_view rk23_min_step_key = "time_integration.node_advection.rk23.min_step";
constexpr std::string_view rk23_max_steps_key = "time_integration.node_advection.rk23.max_steps";
constexpr std::string_view cqmom_time_integration_section = "time_integration.cqmom";
constexpr std::string_view cqmom_source_node_count_evolution_key = "time_integration.cqmom.source_node_count_evolution";
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

bool parse_double(const std::string& text, double& value)
{
    try
    {
        std::size_t consumed = 0;
        value = std::stod(text, &consumed);
        return consumed == text.size() && std::isfinite(value);
    }
    catch (const std::exception&)
    {
        return false;
    }
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

bool is_mixture_component_section(const std::string& section)
{
    if (section.rfind(std::string(mixture_component_section_prefix), 0) != 0)
    {
        return false;
    }
    std::size_t component = 0;
    return parse_size(section.substr(mixture_component_section_prefix.size()), component);
}

bool is_independent_product_dimension_section(const std::string& section)
{
    if (section.rfind(std::string(independent_product_dimension_section_prefix), 0) != 0)
    {
        return false;
    }
    std::size_t dimension = 0;
    return parse_size(section.substr(independent_product_dimension_section_prefix.size()), dimension);
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
                add_error(parsed, ConfigurationErrorCode::syntax_error, source_name, line_number, "Malformed section header.");
                continue;
            }
            section = trim(line_text.substr(1, line_text.size() - 2));
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
                         const std::string& source_name)
{
    for (const auto& item : parsed.scalars)
    {
        bool allowed_key = std::find(allowed.begin(), allowed.end(), item.first) != allowed.end();
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
                             const std::string& source_name)
{
    for (const auto& item : parsed.sections)
    {
        bool allowed_section = std::find(allowed.begin(), allowed.end(), item.first) != allowed.end();
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

bool validate_format_version(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* version = find_entry(parsed, "format_version");
    if (version != nullptr && version->value != "1")
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, version->line, "format_version must be 1.");
    }
    return parsed.errors.empty();
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

bool parse_required_size(ParsedFile& parsed, const std::string& key, const std::string& source_name, std::size_t& value)
{
    const Entry* entry = require_entry(parsed, key, source_name);
    if (entry == nullptr)
    {
        return false;
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

bool parse_optional_positive_double(ParsedFile& parsed, const std::string& key, const std::string& source_name, std::optional<double>& value)
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
    static const std::vector<std::string> fields{"family", "mean",  "standard_deviation", "mu_log", "sigma_log", "alpha", "beta",
                                                 "shape",  "scale", "degrees_of_freedom"};
    if (key.rfind(prefix + ".", 0) != 0)
    {
        return false;
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

void append_time_integration_allowed_keys(std::vector<std::string>& allowed)
{
    allowed.push_back(std::string(transport_scheme_key));
    allowed.push_back(std::string(node_advection_scheme_key));
    allowed.push_back(std::string(rk23_relative_tolerance_key));
    allowed.push_back(std::string(rk23_absolute_tolerance_key));
    allowed.push_back(std::string(rk23_initial_step_key));
    allowed.push_back(std::string(rk23_max_step_key));
    allowed.push_back(std::string(rk23_min_step_key));
    allowed.push_back(std::string(rk23_max_steps_key));
    allowed.push_back(std::string(cqmom_source_node_count_evolution_key));
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
}

std::vector<std::string> with_runtime_allowed_keys(std::vector<std::string> allowed)
{
    append_time_integration_allowed_keys(allowed);
    append_output_allowed_keys(allowed);
    return allowed;
}

std::vector<std::string> time_integration_sections()
{
    return {"time_integration", "time_integration.realizability", "time_integration.node_advection", "time_integration.node_advection.rk23",
            "time_integration.cqmom"};
}

void reject_unknown_qmom_sections(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    std::vector<std::string> allowed = {"phase_space", "moments", "qmom", "initial", "output"};
    const std::vector<std::string> time_sections = time_integration_sections();
    allowed.insert(allowed.end(), time_sections.begin(), time_sections.end());
    if (explicit_moments)
    {
        reject_unknown_sections(parsed, allowed, {"initial."}, source_name);
        return;
    }
    reject_unknown_sections(parsed, allowed, {std::string(mixture_component_section_prefix)}, source_name);
}

void reject_unknown_qmom_keys(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    if (explicit_moments)
    {
        reject_unknown_keys(
            parsed, with_runtime_allowed_keys({"format_version", "phase_space.names", "phase_space.node_counts", "qmom.adaptive", "moments.values"}),
            {"initial."}, source_name);
        return;
    }
    reject_unknown_keys(parsed, with_runtime_allowed_keys(qmom_allowed_keys()), {std::string(mixture_component_section_prefix)}, source_name);
}

std::optional<DistributionMixture> parse_mixture(ParsedFile& parsed, const std::string& source_name)
{
    std::map<std::size_t, Entry> component_sections;
    for (const auto& item : parsed.sections)
    {
        if (!is_mixture_component_section(item.first))
        {
            continue;
        }
        std::size_t component = 0;
        parse_size(item.first.substr(mixture_component_section_prefix.size()), component);
        component_sections.emplace(component, item.second);
    }

    if (component_sections.empty())
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "Mixture initialization requires at least one '[initial.mixture.component.N]' section.");
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

        const std::string prefix = std::string(mixture_component_section_prefix) + std::to_string(item.first);
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

bool has_mixture_component_sections(const ParsedFile& parsed)
{
    for (const auto& item : parsed.sections)
    {
        if (is_mixture_component_section(item.first))
        {
            return true;
        }
    }
    return false;
}

std::vector<std::string> cqmom_explicit_allowed_keys()
{
    return with_runtime_allowed_keys({"format_version", "phase_space.names", "phase_space.node_counts", "phase_space.conditioning_order",
                                      "phase_space.adaptive_branch_counts", "qmom.adaptive", "cqmom.diagnostic_level", "cqmom.moment_storage",
                                      "vandermonde.solver", "moments.extents"});
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
    const std::vector<std::string> time_sections = time_integration_sections();
    allowed.insert(allowed.end(), time_sections.begin(), time_sections.end());
    if (explicit_moments)
    {
        reject_unknown_sections(parsed, allowed, {"initial."}, source_name);
        return;
    }
    reject_unknown_sections(parsed, allowed, {std::string(independent_product_dimension_section_prefix)}, source_name);
}

void reject_unknown_cqmom_keys(ParsedFile& parsed, const bool explicit_moments, const std::string& source_name)
{
    if (explicit_moments)
    {
        reject_unknown_keys(parsed, cqmom_explicit_allowed_keys(), {"initial."}, source_name);
        return;
    }
    reject_unknown_keys(parsed, cqmom_initial_allowed_keys(), {std::string(independent_product_dimension_section_prefix)}, source_name);
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
        const auto lines = parsed.section_lines.find(item.first);
        if (lines != parsed.section_lines.end() && lines->second.size() > 1)
        {
            add_error(parsed, ConfigurationErrorCode::duplicate_entry, source_name, lines->second[1], "Duplicate section '" + item.first + "'.");
        }
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

    std::vector<DistributionSpecification> marginals;
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
        for (const auto& scalar : parsed.scalars)
        {
            if (scalar.first.rfind(prefix + ".", 0) == 0 && !marginal_key_allowed(scalar.first, prefix))
            {
                add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, scalar.second.line, "Unknown key '" + scalar.first + "'.");
            }
        }

        std::optional<DistributionSpecification> distribution = parse_distribution(parsed, prefix, source_name);
        if (distribution.has_value())
        {
            marginals.push_back(*distribution);
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
    product = 1;
    for (const std::size_t extent : extents)
    {
        if (extent == 0 || product > std::numeric_limits<std::size_t>::max() / extent)
        {
            return false;
        }
        product *= extent;
    }
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

    std::vector<double> values(moment_count, 0.0);
    std::vector<bool> supplied(moment_count, false);
    for (const auto& item : parsed.moments)
    {
        if (item.first.size() != extents.size())
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, item.second.line,
                      "Moment index rank does not match the number of phase-space dimensions.");
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
        else if (!values.empty())
        {
            const std::size_t flat = flatten_multi_index(item.first, extents);
            values[flat] = value;
            supplied[flat] = true;
        }
    }

    const std::size_t missing_count = static_cast<std::size_t>(std::count(supplied.begin(), supplied.end(), false));
    if (missing_count != 0)
    {
        add_error(parsed, ConfigurationErrorCode::missing_entry, source_name, 0,
                  "CQMOM moment tensor is incomplete; " + std::to_string(missing_count) + " entries are missing.");
    }
    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return RawMomentTensor{std::move(extents), std::move(values)};
}

template <class T> ConfigurationResult<T> failure_result(ParsedFile&& parsed)
{
    return {std::nullopt, std::move(parsed.errors)};
}

template <class T> ConfigurationResult<T> file_failure(const std::filesystem::path& filename)
{
    return {std::nullopt, {{ConfigurationErrorCode::file_error, filename.string(), 0, "Unable to open configuration file."}}};
}

struct ParsedCqmomCaseSettings
{
    PhaseSpaceHierarchy hierarchy;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
    CqmomMomentStorage moment_storage = CqmomMomentStorage::lean;
    bool moment_storage_configured = false;
    std::size_t moment_storage_line = 0;
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

bool is_realizability_key(const std::string& key)
{
    return key.rfind(std::string(realizability_key_prefix), 0) == 0;
}

const std::vector<std::string>& realizability_allowed_keys()
{
    static const std::vector<std::string> allowed{
        "time_integration.realizability.check_mode", "time_integration.realizability.accept_boundary",    "time_integration.realizability.support",
        "time_integration.realizability.tolerance",  "time_integration.realizability.relative_tolerance", "time_integration.realizability.lower",
        "time_integration.realizability.upper"};
    return allowed;
}

bool is_node_advection_key(const std::string& key)
{
    return key.rfind(std::string(node_advection_key_prefix), 0) == 0;
}

const std::vector<std::string>& node_advection_allowed_keys()
{
    static const std::vector<std::string> allowed{std::string(node_advection_scheme_key),   std::string(rk23_relative_tolerance_key),
                                                  std::string(rk23_absolute_tolerance_key), std::string(rk23_initial_step_key),
                                                  std::string(rk23_max_step_key),           std::string(rk23_min_step_key),
                                                  std::string(rk23_max_steps_key)};
    return allowed;
}

std::optional<CqmomNodeCountEvolution> parse_cqmom_source_node_count_evolution_value(ParsedFile& parsed, const Entry& entry,
                                                                                     const std::string& source_name)
{
    if (entry.value == "preserve_active_counts")
    {
        return CqmomNodeCountEvolution::preserve_active_counts;
    }
    if (entry.value == "allow_growth_to_configured_maximum")
    {
        return CqmomNodeCountEvolution::allow_growth_to_configured_maximum;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line,
              "time_integration.cqmom.source_node_count_evolution must be preserve_active_counts or allow_growth_to_configured_maximum.");
    return std::nullopt;
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

std::optional<RealizabilityCheckMode> parse_realizability_check_mode(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* check_mode = find_entry(parsed, "time_integration.realizability.check_mode");
    if (check_mode != nullptr)
    {
        return parse_realizability_check_mode_value(parsed, *check_mode, "time_integration.realizability.check_mode", source_name);
    }

    return RealizabilityCheckMode::None;
}

std::optional<MomentSupport> parse_moment_support(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = require_entry(parsed, "time_integration.realizability.support", source_name);
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
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry->line,
              "time_integration.realizability.support must be real_line, nonnegative, or bounded_interval.");
    return std::nullopt;
}

std::optional<TransportScheme> parse_transport_scheme_value(ParsedFile& parsed, const Entry& entry, const std::string& source_name)
{
    if (entry.value == "strang")
    {
        return TransportScheme::strang;
    }
    if (entry.value == "cd")
    {
        return TransportScheme::cd;
    }
    if (entry.value == "dc")
    {
        return TransportScheme::dc;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line, "time_integration.scheme must be strang, cd, or dc.");
    return std::nullopt;
}

std::optional<NodeAdvectionScheme> parse_node_advection_scheme_value(ParsedFile& parsed, const Entry& entry, const std::string& source_name)
{
    if (entry.value == "ssp_rk3")
    {
        return NodeAdvectionScheme::ssp_rk3;
    }
    if (entry.value == "adaptive_rk23")
    {
        return NodeAdvectionScheme::adaptive_rk23;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line,
              "time_integration.node_advection.scheme must be ssp_rk3 or adaptive_rk23.");
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

    parse_optional_positive_double(parsed, std::string(rk23_relative_tolerance_key), source_name, options.rk23_options.relative_tolerance);
    parse_optional_positive_double(parsed, std::string(rk23_absolute_tolerance_key), source_name, options.rk23_options.absolute_tolerance);
    parse_optional_positive_double(parsed, std::string(rk23_initial_step_key), source_name, options.rk23_options.initial_step);
    parse_optional_positive_double(parsed, std::string(rk23_max_step_key), source_name, options.rk23_options.max_step);
    parse_optional_positive_double(parsed, std::string(rk23_min_step_key), source_name, options.rk23_options.min_step);
    parse_optional_positive_size(parsed, std::string(rk23_max_steps_key), source_name, options.rk23_options.max_steps);

    const Entry* min_step = find_entry(parsed, std::string(rk23_min_step_key));
    const Entry* max_step = find_entry(parsed, std::string(rk23_max_step_key));
    if (min_step != nullptr && max_step != nullptr && parsed.errors.empty() && options.rk23_options.min_step > options.rk23_options.max_step)
    {
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, min_step->line,
                  "time_integration.node_advection.rk23.min_step must not exceed max_step.");
    }

    if (!parsed.errors.empty())
    {
        return std::nullopt;
    }
    return options;
}

std::optional<RealizabilityPolicy> parse_realizability_policy_from_parsed(ParsedFile& parsed, const std::string& source_name)
{
    if (!has_realizability_section(parsed))
    {
        return RealizabilityPolicy{};
    }

    for (const auto& item : parsed.scalars)
    {
        if (!is_realizability_key(item.first))
        {
            continue;
        }
        const std::vector<std::string>& allowed = realizability_allowed_keys();
        if (std::find(allowed.begin(), allowed.end(), item.first) == allowed.end())
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
        }
    }

    std::optional<RealizabilityCheckMode> check_mode = parse_realizability_check_mode(parsed, source_name);
    if (!check_mode.has_value())
    {
        return std::nullopt;
    }

    RealizabilityPolicy policy;
    policy.check_mode = *check_mode;
    parse_optional_bool(parsed, "time_integration.realizability.accept_boundary", source_name, policy.accept_boundary);
    parse_optional_double(parsed, "time_integration.realizability.tolerance", source_name, policy.options.tolerance);
    parse_optional_bool(parsed, "time_integration.realizability.relative_tolerance", source_name, policy.options.relative_tolerance);
    if (policy.options.tolerance <= 0.0)
    {
        const Entry* entry = find_entry(parsed, "time_integration.realizability.tolerance");
        add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? 0 : entry->line,
                  "time_integration.realizability.tolerance must be positive.");
    }

    if (policy.check_mode != RealizabilityCheckMode::None || find_entry(parsed, "time_integration.realizability.support") != nullptr)
    {
        std::optional<MomentSupport> support = parse_moment_support(parsed, source_name);
        if (support.has_value())
        {
            policy.support = *support;
        }
    }

    if (policy.support == MomentSupport::bounded_interval)
    {
        double lower = 0.0;
        double upper = 0.0;
        const bool parsed_lower = parse_required_double(parsed, "time_integration.realizability.lower", source_name, lower);
        const bool parsed_upper = parse_required_double(parsed, "time_integration.realizability.upper", source_name, upper);
        if (parsed_lower && parsed_upper)
        {
            if (lower >= upper)
            {
                const Entry* entry = find_entry(parsed, "time_integration.realizability.lower");
                add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry == nullptr ? 0 : entry->line,
                          "time_integration.realizability.lower must be less than upper.");
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
        return CqmomDiagnosticLevel::detailed;
    }
    return parse_cqmom_diagnostic_level_value(parsed, *entry, source_name);
}

std::optional<CqmomMomentStorage> parse_cqmom_moment_storage_value(ParsedFile& parsed, const Entry& entry, const std::string& source_name)
{
    if (entry.value == "lean")
    {
        return CqmomMomentStorage::lean;
    }
    if (entry.value == "dense")
    {
        return CqmomMomentStorage::dense;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line, "cqmom.moment_storage must be lean or dense.");
    return std::nullopt;
}

std::optional<CqmomMomentStorage> parse_cqmom_moment_storage(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, "cqmom.moment_storage");
    if (entry == nullptr)
    {
        return CqmomMomentStorage::lean;
    }
    return parse_cqmom_moment_storage_value(parsed, *entry, source_name);
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
    if (entry.value == "bjorck_pereyra")
    {
        return VandermondeSolver::bjorck_pereyra;
    }
    if (entry.value == "bjorck_pereyra_with_dense_fallback")
    {
        return VandermondeSolver::bjorck_pereyra_with_dense_fallback;
    }
    add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, entry.line,
              "vandermonde.solver must be automatic, dense, bjorck_pereyra, or bjorck_pereyra_with_dense_fallback.");
    return std::nullopt;
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

void validate_optional_transport_scheme(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, std::string(transport_scheme_key));
    if (entry != nullptr)
    {
        static_cast<void>(parse_transport_scheme_value(parsed, *entry, source_name));
    }
}

void validate_optional_node_advection_options(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_node_advection_options_from_parsed(parsed, source_name));
}

void validate_optional_realizability_policy(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_realizability_policy_from_parsed(parsed, source_name));
}

void validate_optional_cqmom_source_node_count_evolution(ParsedFile& parsed, const std::string& source_name)
{
    const Entry* entry = find_entry(parsed, std::string(cqmom_source_node_count_evolution_key));
    if (entry != nullptr)
    {
        static_cast<void>(parse_cqmom_source_node_count_evolution_value(parsed, *entry, source_name));
    }
}

void validate_optional_cqmom_diagnostic_level(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_cqmom_diagnostic_level(parsed, source_name));
}

void validate_optional_cqmom_moment_storage(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_cqmom_moment_storage(parsed, source_name));
}

void validate_optional_vandermonde_solver(ParsedFile& parsed, const std::string& source_name)
{
    static_cast<void>(parse_vandermonde_solver(parsed, source_name));
}

ConfigurationResult<ParsedCqmomCaseSettings> parse_cqmom_case_settings(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    const bool explicit_moments = has_explicit_cqmom_moments(parsed);
    reject_unknown_cqmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_cqmom_keys(parsed, explicit_moments, source_name);
    validate_optional_transport_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, source_name);
    validate_optional_cqmom_source_node_count_evolution(parsed, source_name);
    std::optional<CqmomDiagnosticLevel> diagnostic_level = parse_cqmom_diagnostic_level(parsed, source_name);
    const Entry* moment_storage_entry = find_entry(parsed, "cqmom.moment_storage");
    std::optional<CqmomMomentStorage> moment_storage = parse_cqmom_moment_storage(parsed, source_name);
    std::optional<VandermondeSolver> vandermonde_solver = parse_vandermonde_solver(parsed, source_name);
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

    bool adaptive_branch_counts = false;
    const auto adaptive_branch_counts_entry = parsed.scalars.find("phase_space.adaptive_branch_counts");
    if (adaptive_branch_counts_entry != parsed.scalars.end())
    {
        if (adaptive_branch_counts_entry->second.value == "true")
        {
            adaptive_branch_counts = true;
        }
        else if (adaptive_branch_counts_entry->second.value != "false")
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, adaptive_branch_counts_entry->second.line,
                      "adaptive_branch_counts must be true or false.");
        }
    }
    bool adaptive_wheeler = false;
    const auto adaptive_wheeler_entry = parsed.scalars.find("qmom.adaptive");
    if (adaptive_wheeler_entry != parsed.scalars.end())
    {
        if (adaptive_wheeler_entry->second.value == "true")
        {
            adaptive_wheeler = true;
        }
        else if (adaptive_wheeler_entry->second.value == "false")
        {
            adaptive_wheeler = false;
        }
        else
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, adaptive_wheeler_entry->second.line,
                      "qmom.adaptive must be true or false.");
        }
    }

    if (!parsed.errors.empty())
    {
        return failure_result<ParsedCqmomCaseSettings>(std::move(parsed));
    }

    try
    {
        ParsedCqmomCaseSettings settings{PhaseSpaceHierarchy{std::move(node_counts), std::move(conditioning_order)},
                                         {adaptive_wheeler},
                                         adaptive_branch_counts,
                                         diagnostic_level.value_or(CqmomDiagnosticLevel::detailed),
                                         vandermonde_solver.value_or(VandermondeSolver::automatic),
                                         moment_storage.value_or(CqmomMomentStorage::lean),
                                         moment_storage_entry != nullptr,
                                         moment_storage_entry == nullptr ? 0 : moment_storage_entry->line};
        return {std::move(settings), {}};
    }
    catch (const std::exception& error)
    {
        return {std::nullopt, {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0, error.what()}}};
    }
}

} // namespace

ConfigurationResult<MomentSet> initialize_moments(const InitialConditionSpecification& specification, std::string source_name)
{
    if (const auto* explicit_moments = std::get_if<ExplicitMomentInitialization>(&specification))
    {
        return {explicit_moments->moments, {}};
    }
    if (std::holds_alternative<ExplicitTensorInitialization>(specification) ||
        std::holds_alternative<IndependentProductInitialization>(specification))
    {
        return {std::nullopt,
                {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                  "One-dimensional initialization cannot use a multidimensional tensor specification."}}};
    }

    RawMomentResult moments;
    if (const auto* distribution = std::get_if<DistributionInitialization>(&specification))
    {
        moments = compute_raw_moments(distribution->distribution, distribution->maximum_order, distribution->number_density);
    }
    else
    {
        const auto& mixture = std::get<MixtureInitialization>(specification);
        moments = compute_raw_moments(mixture.mixture, mixture.maximum_order);
    }
    if (!moments.ok())
    {
        std::string message = "Analytical initialization moment generation failed: " + moments.message;
        if (moments.failed_component != invalid_component_index)
        {
            message += " Component index: " + std::to_string(moments.failed_component) + ".";
        }
        return {std::nullopt, {{ConfigurationErrorCode::invalid_value, std::move(source_name), 0, std::move(message)}}};
    }
    return {moments.moments, {}};
}

ConfigurationResult<RawMomentTensor> initialize_cqmom_dense_moments(const InitialConditionSpecification& specification, std::string source_name)
{
    if (const auto* explicit_tensor = std::get_if<ExplicitTensorInitialization>(&specification))
    {
        return {explicit_tensor->moments, {}};
    }
    if (const auto* independent_product = std::get_if<IndependentProductInitialization>(&specification))
    {
        if (!independent_product->maximum_orders.has_value())
        {
            return {std::nullopt,
                    {{ConfigurationErrorCode::missing_entry, std::move(source_name), 0,
                      "initial.maximum_orders is required for dense/full tensor independent-product initialization."}}};
        }
        TensorMomentResult moments = compute_raw_moments(independent_product->distribution, *independent_product->maximum_orders);
        if (!moments.ok())
        {
            std::string message = "Independent product initialization moment generation failed: " + moments.message;
            if (moments.failed_dimension != invalid_dimension_index)
            {
                message += " Dimension index: " + std::to_string(moments.failed_dimension) + ".";
            }
            return {std::nullopt, {{ConfigurationErrorCode::invalid_value, std::move(source_name), 0, std::move(message)}}};
        }
        return {moments.moments, {}};
    }
    return {std::nullopt,
            {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
              "Multidimensional initialization requires explicit tensor moments or an independent product specification."}}};
}

ConfigurationResult<CqmomMomentData> initialize_cqmom_moment_data(const InitialConditionSpecification& specification,
                                                                  const PhaseSpaceHierarchy& hierarchy, std::string source_name)
{
    try
    {
        const CqmomMomentLayout layout = create_cqmom_moment_layout(hierarchy);
        if (const auto* explicit_tensor = std::get_if<ExplicitTensorInitialization>(&specification))
        {
            return {extract_cqmom_moments(explicit_tensor->moments, layout), {}};
        }
        if (const auto* independent_product = std::get_if<IndependentProductInitialization>(&specification))
        {
            // Lean initialization requests only the recursive CQMOM moments
            // required by the layout, avoiding unnecessary dense tensor values.
            const RawMomentRequest request{layout.required_indices};
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
                return {std::nullopt, {{ConfigurationErrorCode::invalid_value, std::move(source_name), 0, std::move(message)}}};
            }
            return {create_cqmom_moment_data(layout, selected_moments), {}};
        }
    }
    catch (const std::exception& error)
    {
        return {std::nullopt, {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0, error.what()}}};
    }

    return {std::nullopt,
            {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
              "Lean CQMOM initialization requires explicit tensor moments or an independent product specification."}}};
}

ConfigurationResult<InitialConditionSpecification> parse_qmom_initial_condition(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    const bool explicit_moments = has_explicit_qmom_moments(parsed);
    reject_unknown_qmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_qmom_keys(parsed, explicit_moments, source_name);
    validate_optional_transport_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, source_name);
    validate_optional_cqmom_source_node_count_evolution(parsed, source_name);
    validate_format_version(parsed, source_name);
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
        return {InitialConditionSpecification{ExplicitMomentInitialization{MomentSet(std::move(moments))}}, {}};
    }

    const Entry* type_entry = require_entry(parsed, "initial.type", source_name);
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
    parse_required_size(parsed, "initial.maximum_order", source_name, maximum_order);

    if (type_entry->value == "distribution")
    {
        if (has_mixture_component_sections(parsed))
        {
            add_error(parsed, ConfigurationErrorCode::invalid_configuration, source_name, type_entry->line,
                      "Mixture component sections are only valid when '[initial] type = mixture'.");
        }
        std::optional<DistributionSpecification> distribution = parse_distribution(parsed, source_name);
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

        return {InitialConditionSpecification{DistributionInitialization{*distribution, maximum_order, number_density}}, {}};
    }

    std::optional<DistributionMixture> mixture = parse_mixture(parsed, source_name);
    if (!parsed.errors.empty() || !mixture.has_value())
    {
        return failure_result<InitialConditionSpecification>(std::move(parsed));
    }
    return {InitialConditionSpecification{MixtureInitialization{*mixture, maximum_order}}, {}};
}

ConfigurationResult<InitialConditionSpecification> parse_cqmom_initial_condition(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    const bool explicit_moments = has_explicit_cqmom_moments(parsed);
    reject_unknown_cqmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_cqmom_keys(parsed, explicit_moments, source_name);
    validate_optional_transport_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, source_name);
    validate_optional_cqmom_source_node_count_evolution(parsed, source_name);
    validate_optional_cqmom_diagnostic_level(parsed, source_name);
    validate_optional_cqmom_moment_storage(parsed, source_name);
    validate_optional_vandermonde_solver(parsed, source_name);
    validate_format_version(parsed, source_name);
    validate_optional_cqmom_phase_space(parsed, source_name);

    if (explicit_moments)
    {
        std::optional<RawMomentTensor> tensor = parse_explicit_tensor(parsed, source_name);
        if (!parsed.errors.empty() || !tensor.has_value())
        {
            return failure_result<InitialConditionSpecification>(std::move(parsed));
        }
        return {InitialConditionSpecification{ExplicitTensorInitialization{*tensor}}, {}};
    }

    const Entry* type_entry = require_entry(parsed, "initial.type", source_name);
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
                      "initial.maximum_orders produce tensor extents whose product overflows.");
        }
    }

    std::optional<IndependentProductDistribution> distribution = parse_independent_product(parsed, source_name);
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

    return {InitialConditionSpecification{IndependentProductInitialization{*distribution, std::move(maximum_orders)}}, {}};
}

ConfigurationResult<QmomCase> create_qmom_case(const InitialConditionSpecification& initial_condition, const std::size_t node_count,
                                               const QmomOptions& options, std::string source_name)
{
    if (node_count == 0)
    {
        return {std::nullopt, {{ConfigurationErrorCode::invalid_value, std::move(source_name), 0, "QMOM node_count must be positive."}}};
    }

    ConfigurationResult<MomentSet> initialized = initialize_moments(initial_condition, source_name);
    if (!initialized.ok())
    {
        return {std::nullopt, std::move(initialized.errors)};
    }
    if (node_count <= std::numeric_limits<std::size_t>::max() / 2 && initialized.value->size() < 2 * node_count)
    {
        return {std::nullopt,
                {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0, "QMOM requires at least 2 * node_count moments."}}};
    }
    return {QmomCase{std::move(*initialized.value), node_count, options}, {}};
}

ConfigurationResult<CqmomDenseCase> create_cqmom_dense_case(const InitialConditionSpecification& initial_condition,
                                                            const PhaseSpaceHierarchy& hierarchy, const QmomOptions& qmom_options,
                                                            const bool adaptive_branch_counts, std::string source_name)
{
    ConfigurationResult<RawMomentTensor> initialized = initialize_cqmom_dense_moments(initial_condition, source_name);
    if (!initialized.ok())
    {
        return {std::nullopt, std::move(initialized.errors)};
    }
    if (initialized.value->extents.size() != hierarchy.number_of_dimensions())
    {
        return {std::nullopt,
                {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                  "CQMOM moment tensor rank must match the phase-space dimension count."}}};
    }
    const std::vector<std::size_t>& node_counts = hierarchy.node_counts();
    for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
    {
        if (node_counts[dimension] > std::numeric_limits<std::size_t>::max() / 2 ||
            initialized.value->extents[dimension] < 2 * node_counts[dimension])
        {
            return {std::nullopt,
                    {{ConfigurationErrorCode::invalid_configuration, std::move(source_name), 0,
                      "Each moment extent must be at least 2 * node_count for that phase-space dimension."}}};
        }
    }

    return {CqmomDenseCase{hierarchy, std::move(*initialized.value), {}, qmom_options, adaptive_branch_counts}, {}};
}

ConfigurationResult<CqmomLeanCase> create_cqmom_lean_case(const InitialConditionSpecification& initial_condition,
                                                          const PhaseSpaceHierarchy& hierarchy, const QmomOptions& qmom_options,
                                                          const bool adaptive_branch_counts, std::string source_name)
{
    ConfigurationResult<CqmomMomentData> initialized = initialize_cqmom_moment_data(initial_condition, hierarchy, std::move(source_name));
    if (!initialized.ok())
    {
        return {std::nullopt, std::move(initialized.errors)};
    }
    return {CqmomLeanCase{std::move(*initialized.value), {}, qmom_options, adaptive_branch_counts}, {}};
}

ConfigurationResult<QmomCase> parse_qmom(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    const bool explicit_moments = has_explicit_qmom_moments(parsed);
    reject_unknown_qmom_sections(parsed, explicit_moments, source_name);
    reject_unknown_qmom_keys(parsed, explicit_moments, source_name);
    validate_optional_transport_scheme(parsed, source_name);
    validate_optional_node_advection_options(parsed, source_name);
    validate_optional_realizability_policy(parsed, source_name);
    validate_optional_cqmom_source_node_count_evolution(parsed, source_name);
    validate_format_version(parsed, source_name);
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
    bool adaptive = false;
    const auto adaptive_entry = parsed.scalars.find("qmom.adaptive");
    if (adaptive_entry != parsed.scalars.end())
    {
        if (adaptive_entry->second.value == "true")
        {
            adaptive = true;
        }
        else if (adaptive_entry->second.value != "false")
        {
            add_error(parsed, ConfigurationErrorCode::invalid_value, source_name, adaptive_entry->second.line,
                      "qmom.adaptive must be true or false.");
        }
    }
    if (!parsed.errors.empty())
    {
        return failure_result<QmomCase>(std::move(parsed));
    }

    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_qmom_initial_condition(text, source_name);
    if (!initial_condition.ok())
    {
        return {std::nullopt, std::move(initial_condition.errors)};
    }
    return create_qmom_case(*initial_condition.value, node_counts[0], {adaptive}, std::move(source_name));
}

ConfigurationResult<CqmomDenseCase> parse_cqmom_dense(const std::string_view text, std::string source_name)
{
    ConfigurationResult<ParsedCqmomCaseSettings> settings = parse_cqmom_case_settings(text, source_name);
    if (!settings.ok())
    {
        return {std::nullopt, std::move(settings.errors)};
    }
    if (settings.value->moment_storage_configured && settings.value->moment_storage != CqmomMomentStorage::dense)
    {
        return {std::nullopt,
                {{ConfigurationErrorCode::invalid_value, std::move(source_name), settings.value->moment_storage_line,
                  "parse_cqmom_dense requires cqmom.moment_storage = dense when moment_storage is configured."}}};
    }

    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_cqmom_initial_condition(text, source_name);
    if (!initial_condition.ok())
    {
        return {std::nullopt, std::move(initial_condition.errors)};
    }
    ConfigurationResult<CqmomDenseCase> dense_case =
        create_cqmom_dense_case(*initial_condition.value, settings.value->hierarchy, settings.value->qmom_options,
                                settings.value->adaptive_branch_counts, std::move(source_name));
    if (dense_case.ok())
    {
        dense_case.value->diagnostic_level = settings.value->diagnostic_level;
        dense_case.value->vandermonde_solver = settings.value->vandermonde_solver;
    }
    return dense_case;
}

ConfigurationResult<CqmomLeanCase> parse_cqmom_lean(const std::string_view text, std::string source_name)
{
    ConfigurationResult<ParsedCqmomCaseSettings> settings = parse_cqmom_case_settings(text, source_name);
    if (!settings.ok())
    {
        return {std::nullopt, std::move(settings.errors)};
    }
    if (settings.value->moment_storage_configured && settings.value->moment_storage != CqmomMomentStorage::lean)
    {
        return {std::nullopt,
                {{ConfigurationErrorCode::invalid_value, std::move(source_name), settings.value->moment_storage_line,
                  "parse_cqmom_lean requires cqmom.moment_storage = lean when moment_storage is configured."}}};
    }

    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_cqmom_initial_condition(text, source_name);
    if (!initial_condition.ok())
    {
        return {std::nullopt, std::move(initial_condition.errors)};
    }
    ConfigurationResult<CqmomLeanCase> lean_case =
        create_cqmom_lean_case(*initial_condition.value, settings.value->hierarchy, settings.value->qmom_options,
                               settings.value->adaptive_branch_counts, std::move(source_name));
    if (lean_case.ok())
    {
        lean_case.value->diagnostic_level = settings.value->diagnostic_level;
        lean_case.value->vandermonde_solver = settings.value->vandermonde_solver;
    }
    return lean_case;
}

ConfigurationResult<CqmomCase> parse_cqmom(const std::string_view text, std::string source_name)
{
    ConfigurationResult<ParsedCqmomCaseSettings> settings = parse_cqmom_case_settings(text, source_name);
    if (!settings.ok())
    {
        return {std::nullopt, std::move(settings.errors)};
    }

    ConfigurationResult<InitialConditionSpecification> initial_condition = parse_cqmom_initial_condition(text, source_name);
    if (!initial_condition.ok())
    {
        return {std::nullopt, std::move(initial_condition.errors)};
    }

    if (settings.value->moment_storage == CqmomMomentStorage::dense)
    {
        ConfigurationResult<CqmomDenseCase> dense_case = create_cqmom_dense_case(
            *initial_condition.value, settings.value->hierarchy, settings.value->qmom_options, settings.value->adaptive_branch_counts, source_name);
        if (!dense_case.ok())
        {
            return {std::nullopt, std::move(dense_case.errors)};
        }
        dense_case.value->diagnostic_level = settings.value->diagnostic_level;
        dense_case.value->vandermonde_solver = settings.value->vandermonde_solver;
        return {qmomkit::internal::create_cqmom_case(std::move(*dense_case.value)), {}};
    }

    ConfigurationResult<CqmomLeanCase> lean_case = create_cqmom_lean_case(
        *initial_condition.value, settings.value->hierarchy, settings.value->qmom_options, settings.value->adaptive_branch_counts, source_name);
    if (!lean_case.ok())
    {
        return {std::nullopt, std::move(lean_case.errors)};
    }
    lean_case.value->diagnostic_level = settings.value->diagnostic_level;
    lean_case.value->vandermonde_solver = settings.value->vandermonde_solver;
    return {qmomkit::internal::create_cqmom_case(std::move(*lean_case.value)), {}};
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
        return {RealizabilityPolicy{}, {}};
    }

    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("time_integration.", 0) == 0 && item.first != realizability_section && item.first != node_advection_section &&
            item.first != rk23_section && item.first != cqmom_time_integration_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }

    const std::optional<RealizabilityPolicy> policy = parse_realizability_policy_from_parsed(parsed, source_name);
    if (!parsed.errors.empty() || !policy.has_value())
    {
        return failure_result<RealizabilityPolicy>(std::move(parsed));
    }
    return {*policy, {}};
}

ConfigurationResult<TransportScheme> parse_transport_scheme(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<TransportScheme>(std::move(parsed));
    }

    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("time_integration.", 0) == 0 && item.first != realizability_section && item.first != node_advection_section &&
            item.first != rk23_section && item.first != cqmom_time_integration_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind("time_integration.", 0) != 0)
        {
            continue;
        }
        if (item.first == std::string(transport_scheme_key) || is_realizability_key(item.first) || is_node_advection_key(item.first) ||
            item.first == std::string(cqmom_source_node_count_evolution_key))
        {
            continue;
        }
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
    }

    TransportScheme scheme = TransportScheme::strang;
    const Entry* entry = find_entry(parsed, std::string(transport_scheme_key));
    if (entry != nullptr)
    {
        const std::optional<TransportScheme> parsed_scheme = parse_transport_scheme_value(parsed, *entry, source_name);
        if (parsed_scheme.has_value())
        {
            scheme = *parsed_scheme;
        }
    }

    if (!parsed.errors.empty())
    {
        return failure_result<TransportScheme>(std::move(parsed));
    }
    return {scheme, {}};
}

ConfigurationResult<NodeAdvectionOptions> parse_node_advection_options(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<NodeAdvectionOptions>(std::move(parsed));
    }

    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("time_integration.", 0) == 0 && item.first != realizability_section && item.first != node_advection_section &&
            item.first != rk23_section && item.first != cqmom_time_integration_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind("time_integration.", 0) != 0)
        {
            continue;
        }
        if (item.first == std::string(transport_scheme_key) || is_node_advection_key(item.first) || is_realizability_key(item.first) ||
            item.first == std::string(cqmom_source_node_count_evolution_key))
        {
            continue;
        }
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
    }

    const std::optional<NodeAdvectionOptions> options = parse_node_advection_options_from_parsed(parsed, source_name);

    if (!parsed.errors.empty() || !options.has_value())
    {
        return failure_result<NodeAdvectionOptions>(std::move(parsed));
    }
    return {*options, {}};
}

ConfigurationResult<NodeAdvectionScheme> parse_node_advection_scheme(const std::string_view text, std::string source_name)
{
    ConfigurationResult<NodeAdvectionOptions> options = parse_node_advection_options(text, std::move(source_name));
    if (!options.ok())
    {
        return {std::nullopt, std::move(options.errors)};
    }
    return {options.value->scheme, {}};
}

ConfigurationResult<CqmomNodeCountEvolution> parse_cqmom_source_node_count_evolution(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<CqmomNodeCountEvolution>(std::move(parsed));
    }

    for (const auto& item : parsed.sections)
    {
        if (item.first.rfind("time_integration.", 0) == 0 && item.first != realizability_section && item.first != node_advection_section &&
            item.first != rk23_section && item.first != cqmom_time_integration_section)
        {
            add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown section '" + item.first + "'.");
        }
    }
    for (const auto& item : parsed.scalars)
    {
        if (item.first.rfind("time_integration.", 0) != 0)
        {
            continue;
        }
        if (item.first == std::string(transport_scheme_key) || is_node_advection_key(item.first) || is_realizability_key(item.first))
        {
            continue;
        }
        if (item.first == std::string(cqmom_source_node_count_evolution_key))
        {
            continue;
        }
        add_error(parsed, ConfigurationErrorCode::unknown_entry, source_name, item.second.line, "Unknown key '" + item.first + "'.");
    }

    CqmomNodeCountEvolution policy = CqmomNodeCountEvolution::allow_growth_to_configured_maximum;
    const Entry* entry = find_entry(parsed, std::string(cqmom_source_node_count_evolution_key));
    if (entry != nullptr)
    {
        const std::optional<CqmomNodeCountEvolution> parsed_policy = parse_cqmom_source_node_count_evolution_value(parsed, *entry, source_name);
        if (parsed_policy.has_value())
        {
            policy = *parsed_policy;
        }
    }

    if (!parsed.errors.empty())
    {
        return failure_result<CqmomNodeCountEvolution>(std::move(parsed));
    }
    return {policy, {}};
}

ConfigurationResult<OutputOptions> parse_output_options(const std::string_view text, std::string source_name)
{
    ParsedFile parsed = parse_file(text, source_name);
    if (!parsed.errors.empty())
    {
        return failure_result<OutputOptions>(std::move(parsed));
    }
    if (!has_output_section(parsed))
    {
        return {OutputOptions{}, {}};
    }

    static const std::vector<std::string> allowed = {"output.directory", "output.write_moments", "output.write_quadrature", "output.interval"};
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

    OutputOptions options;
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

    if (!parsed.errors.empty())
    {
        return failure_result<OutputOptions>(std::move(parsed));
    }
    return {options, {}};
}

ConfigurationResult<QmomCase> load_qmom(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<QmomCase>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_qmom(contents.str(), filename.string());
}

ConfigurationResult<CqmomDenseCase> load_cqmom_dense(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<CqmomDenseCase>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_cqmom_dense(contents.str(), filename.string());
}

ConfigurationResult<CqmomLeanCase> load_cqmom_lean(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<CqmomLeanCase>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_cqmom_lean(contents.str(), filename.string());
}

ConfigurationResult<CqmomCase> load_cqmom(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<CqmomCase>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_cqmom(contents.str(), filename.string());
}

ConfigurationResult<OutputOptions> load_output_options(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<OutputOptions>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_output_options(contents.str(), filename.string());
}

ConfigurationResult<NodeAdvectionOptions> load_node_advection_options(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<NodeAdvectionOptions>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_node_advection_options(contents.str(), filename.string());
}

ConfigurationResult<NodeAdvectionScheme> load_node_advection_scheme(const std::filesystem::path& filename)
{
    ConfigurationResult<NodeAdvectionOptions> options = load_node_advection_options(filename);
    if (!options.ok())
    {
        return {std::nullopt, std::move(options.errors)};
    }
    return {options.value->scheme, {}};
}

ConfigurationResult<CqmomNodeCountEvolution> load_cqmom_source_node_count_evolution(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        return file_failure<CqmomNodeCountEvolution>(filename);
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return parse_cqmom_source_node_count_evolution(contents.str(), filename.string());
}

} // namespace qmomkit::config
