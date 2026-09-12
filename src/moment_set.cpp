#include "qmomkit/moment_set.hpp"

#include <stdexcept>
#include <utility>

namespace qmomkit
{

MomentSet::MomentSet(std::vector<double> moments) : moments_(std::move(moments)) {}

MomentSet::MomentSet(std::initializer_list<double> moments) : moments_(moments) {}

std::size_t MomentSet::size() const noexcept
{
    return moments_.size();
}

bool MomentSet::empty() const noexcept
{
    return moments_.empty();
}

const std::vector<double>& MomentSet::values() const noexcept
{
    return moments_;
}

double MomentSet::operator[](const std::size_t moment_order) const
{
    if (moment_order >= moments_.size())
    {
        throw std::out_of_range("Moment order is outside the stored sequence.");
    }
    return moments_[moment_order];
}

} // namespace qmomkit
