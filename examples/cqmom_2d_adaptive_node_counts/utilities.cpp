#include "utilities.hpp"

#include <cstddef>
#include <iomanip>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>

namespace cqmom_2d_adaptive_node_counts
{
namespace
{

// A count that does not exist before a branch is created or after it is
// removed is shown as "-".
std::string format_count(const std::optional<std::size_t>& count)
{
    return count.has_value() ? std::to_string(*count) : "-";
}

const char* kind_name(const qmomkit::CqmomNodeCountTransitionKind kind)
{
    switch (kind)
    {
    case qmomkit::CqmomNodeCountTransitionKind::reduced:
        return "reduced";
    case qmomkit::CqmomNodeCountTransitionKind::grew:
        return "grew";
    case qmomkit::CqmomNodeCountTransitionKind::created:
        return "created";
    case qmomkit::CqmomNodeCountTransitionKind::removed:
        return "removed";
    }
    return "";
}

} // namespace

void print_transition(std::ostream& stream, const double time, const qmomkit::CqmomNodeCountTransition& transition)
{
    std::ostringstream line;
    line << "  t = " << std::fixed << std::setprecision(3) << time << "  branch [";
    for (std::size_t level = 0; level < transition.parent_path.size(); ++level)
    {
        line << (level > 0 ? "," : "") << transition.parent_path.node_index_at_level(level);
    }
    line << "]  " << format_count(transition.active_node_count_before) << " -> " << format_count(transition.active_node_count_after) << " nodes ("
         << kind_name(transition.kind) << ")\n";
    stream << line.str();
}

} // namespace cqmom_2d_adaptive_node_counts
