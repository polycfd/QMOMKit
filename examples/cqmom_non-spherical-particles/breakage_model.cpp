#include "breakage_model.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cqmom_non_spherical_particles
{

qmomkit::CqmomMomentSource binary_breakage_source()
{
    return [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> source(moments.size(), 0.0);
        for (std::size_t index = 0; index < moments.size(); ++index)
        {
            if (moments.index(index).size() != 5)
            {
                throw std::invalid_argument("The 5D binary breakage source received a moment index with the wrong rank.");
            }
        }
        return qmomkit::CqmomMomentRate{std::move(source)};
    };
}

void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, double, std::mt19937_64&)
{
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        if (sample.size() != 5)
        {
            throw std::runtime_error("The 5D MC particle example expects samples ordered as (d, chi, u, v, w).");
        }
    }
}

} // namespace cqmom_non_spherical_particles
