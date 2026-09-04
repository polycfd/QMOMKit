#include "conditioning_report.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cqmom_2d_borderline_conditioning
{

namespace
{

constexpr std::size_t solver_count = 4;

std::string optional_value(const std::optional<double>& value)
{
    if (!value.has_value() || !std::isfinite(*value))
    {
        return "-";
    }

    std::ostringstream stream;
    stream << std::scientific << std::setprecision(6) << *value;
    return stream.str();
}

std::string order_name(const std::vector<std::size_t>& conditioning_order)
{
    return conditioning_order == std::vector<std::size_t>{0, 1} ? "x_then_y" : "y_then_x";
}

} // namespace

void write_header(std::ostream& stream)
{
    stream << "# mode delta conditioning_order requested_solver selected_solver success status active_counts diagnostic_level parent_path "
              "condition_estimate_inf sigma_min_over_sigma_max local_svd_rank condition_number_2 residual_norm moment_error failed_dimension "
              "message\n";
}

void write_row(std::ostream& stream, const RunRecord& record)
{
    stream << record.mode << '\t' << std::scientific << std::setprecision(6) << record.delta << '\t' << record.conditioning_order << '\t'
           << record.requested_solver << '\t' << record.selected_solver << '\t' << (record.success ? "true" : "false") << '\t' << record.status
           << '\t' << record.active_counts << '\t' << record.diagnostic_level << '\t' << record.parent_path << '\t'
           << optional_value(record.condition_estimate_inf) << '\t' << optional_value(std::optional<double>{record.matrix.ratio}) << '\t'
           << record.matrix.rank << '\t' << optional_value(record.matrix.condition_number) << '\t' << optional_value(record.residual_norm) << '\t'
           << optional_value(record.moment_error) << '\t' << record.failed_dimension << '\t' << record.message << '\n';
}

std::string report_number(const double value)
{
    if (std::isinf(value))
    {
        return "inf";
    }
    if (!std::isfinite(value))
    {
        return "-";
    }

    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

std::string report_number(const std::optional<double>& value)
{
    return value.has_value() ? report_number(*value) : "-";
}

std::string outcome_code(const RunRecord& record)
{
    if (record.success)
    {
        return "OK";
    }
    if (record.status == "conditional_moment_failed")
    {
        return "FAIL-C";
    }
    if (record.status == "qmom_failed")
    {
        return "FAIL-W";
    }
    return "FAIL";
}

std::string failure_stage_label(const RunRecord& record)
{
    if (record.success)
    {
        return "-";
    }
    if (record.status == "conditional_moment_failed")
    {
        return "conditional";
    }
    if (record.status == "qmom_failed")
    {
        return "Wheeler/QMOM";
    }
    return "other";
}

std::string concise_solver_name(const std::string& label)
{
    if (label == "dense")
    {
        return "Dense";
    }
    if (label == "bjorck_pereyra")
    {
        return "BP";
    }
    if (label == "pivoted_qr")
    {
        return "QR";
    }
    if (label == "svd")
    {
        return "SVD";
    }
    return label;
}

struct FailureGroup
{
    double delta = 0.0;
    std::string conditioning_order;
    std::vector<std::string> affected_policies;
    std::string stage;
    std::string message;
};

std::vector<FailureGroup> group_failure_details(const std::vector<RunRecord>& records)
{
    std::vector<FailureGroup> groups;
    for (const RunRecord& record : records)
    {
        const bool fixed_failure = record.mode == "fixed" && !record.success;
        if (!fixed_failure)
        {
            continue;
        }

        const std::string stage = failure_stage_label(record);
        auto group = std::find_if(groups.begin(), groups.end(),
                                  [&](const FailureGroup& candidate)
                                  {
                                      return candidate.delta == record.delta && candidate.conditioning_order == record.conditioning_order &&
                                             candidate.stage == stage && candidate.message == record.message;
                                  });
        if (group == groups.end())
        {
            groups.push_back({record.delta, record.conditioning_order, {}, stage, record.message});
            group = std::prev(groups.end());
        }
        group->affected_policies.push_back(concise_solver_name(record.requested_solver));
    }
    return groups;
}

std::string join_policies(const std::vector<std::string>& policies)
{
    std::ostringstream stream;
    for (std::size_t index = 0; index < policies.size(); ++index)
    {
        if (index > 0)
        {
            stream << ", ";
        }
        stream << policies[index];
    }
    return stream.str();
}

const RunRecord* find_record(const std::vector<RunRecord>& records, const std::string& mode, const double delta,
                             const std::string& conditioning_order, const std::string& requested_solver)
{
    for (const RunRecord& record : records)
    {
        if (record.mode == mode && record.delta == delta && record.conditioning_order == conditioning_order &&
            record.requested_solver == requested_solver)
        {
            return &record;
        }
    }
    return nullptr;
}

const RunRecord& require_record(const std::vector<RunRecord>& records, const std::string& mode, const double delta,
                                const std::string& conditioning_order, const std::string& requested_solver)
{
    const RunRecord* record = find_record(records, mode, delta, conditioning_order, requested_solver);
    if (record == nullptr)
    {
        throw std::runtime_error("The report input is missing an expected reconstruction result.");
    }
    return *record;
}

const char* solver_column_name(const SolverChoice& choice)
{
    if (choice.solver == qmomkit::VandermondeSolver::dense)
    {
        return "Dense";
    }
    if (choice.solver == qmomkit::VandermondeSolver::bjorck_pereyra)
    {
        return "BP";
    }
    if (choice.solver == qmomkit::VandermondeSolver::pivoted_qr)
    {
        return "QR";
    }
    return "SVD";
}

std::optional<double> maximum_successful_moment_error(const std::vector<RunRecord>& records, const double delta,
                                                      const std::string& conditioning_order)
{
    std::optional<double> maximum_error;
    for (const SolverChoice& choice : solver_choices)
    {
        const RunRecord& record = require_record(records, "fixed", delta, conditioning_order, choice.label);
        if (record.success && record.moment_error.has_value())
        {
            maximum_error = maximum_error.has_value() ? std::max(*maximum_error, *record.moment_error) : record.moment_error;
        }
    }
    return maximum_error;
}

std::string report_section(const std::string& report, const std::string& heading)
{
    const std::size_t begin = report.find(heading);
    if (begin == std::string::npos)
    {
        throw std::runtime_error("The Markdown report is missing the section: " + heading);
    }
    const std::size_t next_heading = report.find("\n#", begin + heading.size());
    return report.substr(begin, next_heading == std::string::npos ? std::string::npos : next_heading - begin);
}

void validate_report_file(const std::filesystem::path& filename)
{
    std::ifstream input(filename);
    if (!input)
    {
        throw std::runtime_error("Could not reopen output/report.md for validation.");
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.good() && !input.eof())
    {
        throw std::runtime_error("Could not read output/report.md for validation.");
    }
    const std::string report = contents.str();

    for (const std::string& heading : {"# 2D CQMOM Borderline-Conditioning Report", "## Purpose Of The Example", "## Test Population And Moment Set",
                                       "## Why The Reconstruction Is Numerically Difficult", "## Quantities Used To Assess The Reconstruction",
                                       "## Results For Fixed Node Counts", "## Effect Of The Vandermonde Solver", "## Effect Of Conditioning Order",
                                       "## Failed Reconstructions", "## Results With Adaptive Node Counts", "## Main Observations"})
    {
        if (report.find(heading) == std::string::npos)
        {
            throw std::runtime_error("The Markdown report is missing the heading: " + heading);
        }
    }
    if (report.find("[diagnostics.txt](diagnostics.txt)") == std::string::npos)
    {
        throw std::runtime_error("The Markdown report is missing its diagnostics link.");
    }
    for (const std::string& expression : {"$\\delta$", "$r_\\sigma$", "$\\kappa_2$"})
    {
        if (report.find(expression) == std::string::npos)
        {
            throw std::runtime_error("The Markdown report is missing a metric expression: " + expression);
        }
    }
    for (const std::string& obsolete_heading : {"| sigma_min/sigma_max |", "| kappa_2 |", "| max moment error |"})
    {
        if (report.find(obsolete_heading) != std::string::npos)
        {
            throw std::runtime_error("The Markdown report contains an obsolete metric heading: " + obsolete_heading);
        }
    }
    if (report.find("four direct solvers") == std::string::npos || report.find("does not establish a universal solver ranking") == std::string::npos)
    {
        throw std::runtime_error("The Markdown report is missing its solver-comparison conclusion.");
    }
    for (const std::string& conditioning_order : {"(0,1)", "(1,0)"})
    {
        const std::string section = report_section(report, conditioning_order == "(0,1)" ? "### x_then_y" : "### y_then_x");
        for (const double delta : delta_values)
        {
            const std::string row_prefix = "| " + report_number(delta) + " |";
            if (section.find(row_prefix) == std::string::npos)
            {
                throw std::runtime_error("The Markdown report is missing a fixed-count delta row.");
            }
        }
        for (const char* column : {"Dense", "BP", "QR", "SVD"})
        {
            if (section.find(std::string("| ") + column + " |") == std::string::npos)
            {
                throw std::runtime_error("The Markdown report is missing a solver column.");
            }
        }
    }
}

void write_order_report(std::ostream& stream, const std::vector<RunRecord>& records, const std::string& conditioning_order, const std::string& title)
{
    stream << "### " << title << "\n\n";
    stream << "| delta | relative singular value | local rank | condition number |";
    for (const std::size_t solver_index : report_solver_indices)
    {
        stream << ' ' << solver_column_name(solver_choices[solver_index]) << " |";
    }
    stream << " worst successful moment error |\n";
    stream << "|---:|---:|---:|---:|---|---|---|---|---:|\n";
    for (const double delta : delta_values)
    {
        const RunRecord& first = require_record(records, "fixed", delta, conditioning_order, solver_choices.front().label);
        stream << "| " << report_number(delta) << " | " << report_number(first.matrix.ratio) << " | " << first.matrix.rank << "/2 | "
               << report_number(first.matrix.condition_number) << " |";
        for (const std::size_t solver_index : report_solver_indices)
        {
            stream << ' ' << outcome_code(require_record(records, "fixed", delta, conditioning_order, solver_choices[solver_index].label)) << " |";
        }
        stream << ' ' << report_number(maximum_successful_moment_error(records, delta, conditioning_order)) << " |\n";
    }
    stream << '\n';
}

void write_diagnostics_file(const std::vector<RunRecord>& records, const std::filesystem::path& filename)
{
    std::ofstream output(filename);
    if (!output)
    {
        throw std::runtime_error("Could not open output/diagnostics.txt.");
    }
    write_header(output);
    for (const RunRecord& record : records)
    {
        write_row(output, record);
    }
    output.flush();
    if (!output)
    {
        throw std::runtime_error("Could not finish writing output/diagnostics.txt.");
    }
}

void write_markdown_report(const std::vector<RunRecord>& records, const std::filesystem::path& filename)
{
    std::ofstream report(filename);
    if (!report)
    {
        throw std::runtime_error("Could not open output/report.md.");
    }

    std::size_t conditional_failures = 0;
    std::size_t qmom_failures = 0;
    std::size_t successful_fixed = 0;
    for (const RunRecord& record : records)
    {
        if (record.success)
        {
            successful_fixed += record.mode == "fixed" ? 1 : 0;
        }
        else if (record.status == "conditional_moment_failed")
        {
            ++conditional_failures;
        }
        else if (record.status == "qmom_failed")
        {
            ++qmom_failures;
        }
    }

    std::optional<double> first_x_rank_deficient_delta;
    bool all_direct_x_rank_failure_conditional = true;
    bool y_parent_remained_full_rank = true;
    for (const double delta : delta_values)
    {
        const RunRecord& x_record = require_record(records, "fixed", delta, "(0,1)", solver_choices.front().label);
        const RunRecord& y_record = require_record(records, "fixed", delta, "(1,0)", solver_choices.front().label);
        if (!first_x_rank_deficient_delta.has_value() && x_record.matrix.rank < 2)
        {
            first_x_rank_deficient_delta = delta;
            for (const std::size_t solver_index : report_solver_indices)
            {
                all_direct_x_rank_failure_conditional =
                    all_direct_x_rank_failure_conditional &&
                    require_record(records, "fixed", delta, "(0,1)", solver_choices[solver_index].label).status == "conditional_moment_failed";
            }
        }
        y_parent_remained_full_rank = y_parent_remained_full_rank && y_record.matrix.rank == 2;
    }

    const std::vector<FailureGroup> failure_groups = group_failure_details(records);
    report << "# 2D CQMOM Borderline-Conditioning Report\n\n"
              "## Purpose Of The Example\n\n"
              "This example follows a fixed two-dimensional population while the two x-nodes move together. It shows how the conditioning order "
              "changes the stage at which the reconstruction becomes difficult, how four Vandermonde solvers compare on the same input, and how "
              "allowing fewer nodes changes the result at exact degeneracy.\n\n";
    report << "## Test Population And Moment Set\n\n"
              "For separation $\\delta$, the four atoms are\n\n"
              "$$\n"
              "\\left(-\\frac{\\delta}{2},-1\\right),\\;w=0.10,\\qquad "
              "\\left(-\\frac{\\delta}{2},+1\\right),\\;w=0.30,\n"
              "$$\n\n"
              "$$\n"
              "\\left(+\\frac{\\delta}{2},-1\\right),\\;w=0.15,\\qquad "
              "\\left(+\\frac{\\delta}{2},+1\\right),\\;w=0.45.\n"
              "$$\n\n"
              "The raw-moment tensor contains $M_{p,q}$ for $0\\le p,q\\le3$, computed by direct summation:\n\n"
              "$$\n"
              "M_{p,q}=\\sum_j w_j x_j^p y_j^q.\n"
              "$$\n\n"
              "Useful checks are $M_{0,0}=1$, $M_{1,0}=0.1\\delta$, $M_{0,1}=0.5$, and $M_{2,0}=0.25\\delta^2$.\n\n";
    report
        << "## Why The Reconstruction Is Numerically Difficult\n\n"
           "In x_then_y, the weighted parent system is\n\n"
           "$$\n"
           "A_x=\\begin{bmatrix}0.4&0.6\\\\-0.2\\delta&0.3\\delta\\end{bmatrix}.\n"
           "$$\n\n"
           "Its infinity-norm condition number grows approximately as\n\n"
           "$$\n"
           "\\kappa_\\infty(A_x)\\sim\\frac{2.5}{\\delta}.\n"
           "$$\n\n"
           "The y-first parent system remains well conditioned. Reversing the order therefore moves the weak x-reconstruction later in the recursive "
           "calculation; it does not remove the near-degeneracy.\n\n";
    report << "## Quantities Used To Assess The Reconstruction\n\n"
              "The example reports four different quantities. The relative singular value measures column independence:\n\n"
              "$$\n"
              "r_\\sigma=\\frac{\\sigma_{\\min}}{\\sigma_{\\max}}.\n"
              "$$\n\n"
              "The local 2-norm condition number measures sensitivity:\n\n"
              "$$\n"
              "\\kappa_2=\\frac{\\sigma_{\\max}}{\\sigma_{\\min}}=\\frac{1}{r_\\sigma}.\n"
              "$$\n\n"
              "For a successful quadrature, the worst scaled moment error is\n\n"
              "$$\n"
              "E_M=\\max_{p,q}\\frac{|\\widehat M_{p,q}-M_{p,q}|}{\\max(1,|M_{p,q}|)}.\n"
              "$$\n\n"
              "A residual checks how closely the computed weights satisfy the represented linear system. A small residual does not imply that the "
              "system is insensitive. The local rank indicator uses the relative threshold $1\\times10^{-14}$; these matrix quantities are computed "
              "from the prescribed population, not copied from an internal QR or SVD factorization.\n\n";

    report << "## Results For Fixed Node Counts\n\n"
              "The fixed sweep contains 64 reconstructions. "
           << successful_fixed << " succeeded, " << conditional_failures << " failed in a conditional Vandermonde solve, and " << qmom_failures
           << " failed in a Wheeler or broader QMOM stage.\n\n"
              "| outcome | fixed-count rows |\n"
              "|---|---:|\n"
              "| successful | "
           << successful_fixed
           << " |\n"
              "| conditional Vandermonde failure | "
           << conditional_failures
           << " |\n"
              "| Wheeler/QMOM failure | "
           << qmom_failures << " |\n\n";
    if (first_x_rank_deficient_delta.has_value())
    {
        report << "The first rank-deficient x-parent occurs at $\\delta = " << report_number(*first_x_rank_deficient_delta)
               << "$. At that separation, "
               << (all_direct_x_rank_failure_conditional ? "all four direct solvers fail in the conditional Vandermonde solve."
                                                         : "the four direct solvers do not all fail at the same stage.")
               << "\n\n";
    }
    report << (y_parent_remained_full_rank
                   ? "In y_then_x, the y-parent remains full rank throughout the sweep. The difficult x-reconstruction "
                     "therefore appears later as a Wheeler/QMOM failure at exact coincidence.\n\n"
                   : "In y_then_x, the y-parent is not full rank for every separation; interpret the order comparison together "
                     "with the local metrics.\n\n");

    report << "## Effect Of The Vandermonde Solver\n\n"
              "The four direct solvers are compared to determine whether a failed reconstruction depends on the linear solver. If their outcomes "
              "disagree, the case may be sensitive to solver-specific rank or pivot decisions. If they fail at the same stage, the failure is not "
              "fixed by changing the solver.\n\n"
              "The outcome codes used in the tables are:\n\n"
              "| Code | Meaning |\n"
              "|---|---|\n"
              "| OK | Successful reconstruction |\n"
              "| FAIL-C | Conditional-moment / conditional Vandermonde failure |\n"
              "| FAIL-W | Wheeler or broader QMOM reconstruction failure |\n"
              "| FAIL | Other failed reconstruction |\n\n";

    report << "## Effect Of Conditioning Order\n\n"
              "The tables below keep the same separation values and solver columns while changing only the conditioning order. The x-first order "
              "reaches the weak conditional solve early. The y-first order starts from the well-separated parent system and reaches the weak "
              "x-direction "
              "later in the recursive reconstruction.\n\n";
    write_order_report(report, records, "(0,1)", "x_then_y — Clustered Coordinate First");
    write_order_report(report, records, "(1,0)", "y_then_x — Well-Separated Coordinate First");

    report << "## Failed Reconstructions\n\n"
              "The reported stage distinguishes a conditional Vandermonde failure from a Wheeler/QMOM failure. The latter covers a failed Wheeler "
              "recurrence, either in a marginal reconstruction or later in the recursive calculation. Identical failures are combined into one row "
              "so the table shows which solver policies share the same observed result.\n\n"
              "| delta | order | solvers with this result | stage | message |\n"
              "|---:|---|---|---|---|\n";
    for (const FailureGroup& group : failure_groups)
    {
        report << "| " << report_number(group.delta) << " | " << group.conditioning_order << " | " << join_policies(group.affected_policies) << " | "
               << group.stage << " | " << group.message << " |\n";
    }
    report << '\n';

    report << "## Results With Adaptive Node Counts\n\n"
              "| order | outcome | active/requested branch counts | interpretation |\n"
              "|---|---|---|---|\n";
    for (const std::string& conditioning_order : {"(0,1)", "(1,0)"})
    {
        const RunRecord& record = require_record(records, "adaptive_reference", 0.0, conditioning_order, "automatic");
        const std::string interpretation = conditioning_order == "(0,1)"
                                               ? "The root x-quadrature reduces to one active node; the remaining y-branch retains two."
                                               : "Both parent y-nodes remain active; each conditional x-branch reduces to one node.";
        report << "| " << order_name(conditioning_order == "(0,1)" ? std::vector<std::size_t>{0, 1} : std::vector<std::size_t>{1, 0}) << " | "
               << outcome_code(record) << " | " << record.active_counts << " | " << interpretation << " |\n";
    }
    report << "\nAt $\\delta=0$, the adaptive cases succeed by constructing lower-node quadratures supported by the moments. They do not make the "
              "singular "
              "fixed two-node request well conditioned. Node reduction occurs in Wheeler inversion before the conditional Vandermonde solve; QR and "
              "SVD do not trigger it.\n\n";

    report
        << "## Main Observations\n\n"
           "1. As $\\delta$ decreases, $r_\\sigma$ decreases and $\\kappa_2$ increases for the x-parent.\n"
           "2. At $\\delta = 2.00e-14$, all four direct solvers fail in the same conditional Vandermonde solve.\n"
           "3. At $\\delta=0$, x_then_y fails after the x-parent collapses, while y_then_x reaches a Wheeler failure after its well-conditioned "
           "y-parent step.\n"
           "4. A small residual or small moment error describes the current computed result; neither proves robustness as the nodes coalesce.\n"
           "5. The adaptive cases succeed because they use fewer active nodes. They are lower-node reconstructions, not solver fixes.\n\n"
           "These results describe this prescribed population. This example does not establish a universal solver ranking. The complete dataset is "
           "available in [diagnostics.txt](diagnostics.txt).\n";

    report.flush();
    if (!report)
    {
        throw std::runtime_error("Could not finish writing output/report.md.");
    }
}

void validate_successful_rows(const std::vector<RunRecord>& records)
{
    const std::size_t principal_rows = delta_values.size() * 2 * solver_choices.size();
    if (records.size() != principal_rows + 2)
    {
        throw std::runtime_error("The diagnostic sweep did not produce the expected number of rows.");
    }

    bool saw_x_then_y = false;
    bool saw_y_then_x = false;
    std::array<bool, solver_count> saw_solver{};
    bool saw_healthy_success = false;
    bool saw_non_success = false;
    for (const RunRecord& record : records)
    {
        saw_x_then_y = saw_x_then_y || record.conditioning_order == "(0,1)";
        saw_y_then_x = saw_y_then_x || record.conditioning_order == "(1,0)";
        for (std::size_t solver = 0; solver < solver_choices.size(); ++solver)
        {
            saw_solver[solver] = saw_solver[solver] || record.requested_solver == solver_choices[solver].label;
        }
        if (record.mode == "fixed" && record.delta == delta_values.front())
        {
            saw_healthy_success = saw_healthy_success || record.success;
        }
        saw_non_success = saw_non_success || !record.success;
    }

    if (!saw_x_then_y || !saw_y_then_x || !saw_healthy_success || !saw_non_success)
    {
        throw std::runtime_error("The diagnostic sweep did not exercise the required order and status cases.");
    }
    for (const bool exercised : saw_solver)
    {
        if (!exercised)
        {
            throw std::runtime_error("The diagnostic sweep did not exercise every requested solver.");
        }
    }
}

void validate_report_inputs(const std::vector<RunRecord>& records)
{
    for (const double delta : delta_values)
    {
        for (const std::string& conditioning_order : {"(0,1)", "(1,0)"})
        {
            for (const SolverChoice& choice : solver_choices)
            {
                (void)require_record(records, "fixed", delta, conditioning_order, choice.label);
            }
        }
    }
    for (const std::string& conditioning_order : {"(0,1)", "(1,0)"})
    {
        (void)require_record(records, "adaptive_reference", 0.0, conditioning_order, "automatic");
    }
}

void print_terminal_summary(const std::vector<RunRecord>& records)
{
    std::size_t successful = 0;
    std::size_t conditional_failures = 0;
    std::size_t qmom_failures = 0;
    std::size_t other_failures = 0;
    std::size_t fixed_successful = 0;
    std::size_t fixed_conditional_failures = 0;
    std::size_t fixed_qmom_failures = 0;
    std::size_t fixed_other_failures = 0;
    for (const RunRecord& record : records)
    {
        if (record.success)
        {
            ++successful;
            fixed_successful += record.mode == "fixed" ? 1 : 0;
        }
        else if (record.status == "conditional_moment_failed")
        {
            ++conditional_failures;
            fixed_conditional_failures += record.mode == "fixed" ? 1 : 0;
        }
        else if (record.status == "qmom_failed")
        {
            ++qmom_failures;
            fixed_qmom_failures += record.mode == "fixed" ? 1 : 0;
        }
        else
        {
            ++other_failures;
            fixed_other_failures += record.mode == "fixed" ? 1 : 0;
        }
    }

    const RunRecord& adaptive_x_then_y = require_record(records, "adaptive_reference", 0.0, "(0,1)", "automatic");
    const RunRecord& adaptive_y_then_x = require_record(records, "adaptive_reference", 0.0, "(1,0)", "automatic");
    std::cout << "\n2D CQMOM borderline-conditioning sweep complete\n\n"
              << "All results (" << records.size() << " runs):\n"
              << "  successful reconstructions   : " << successful << "\n"
              << "  conditional failures         : " << conditional_failures << "\n"
              << "  Wheeler/QMOM failures        : " << qmom_failures << "\n"
              << "  other failures               : " << other_failures << "\n\n"
              << "Fixed-count results (" << records.size() - 2 << " runs):\n"
              << "  successful reconstructions   : " << fixed_successful << "\n"
              << "  conditional failures         : " << fixed_conditional_failures << "\n"
              << "  Wheeler/QMOM failures        : " << fixed_qmom_failures << "\n"
              << "  other failures               : " << fixed_other_failures << "\n\n"
              << "Adaptive delta=0 references:\n"
              << "  x_then_y : " << outcome_code(adaptive_x_then_y) << ", " << adaptive_x_then_y.active_counts << "\n"
              << "  y_then_x : " << outcome_code(adaptive_y_then_x) << ", " << adaptive_y_then_x.active_counts << "\n\n"
              << "Human-readable report : output/report.md\n"
              << "Complete diagnostics   : output/diagnostics.txt\n";
}

} // namespace cqmom_2d_borderline_conditioning
