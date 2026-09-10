#include "qmomkit/moment_reconstruction.hpp"

#include "internal/cqmom_validation.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/numeric_utils.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{

[[nodiscard]] double integer_power(const double base, const std::size_t exponent)
{
    double result = 1.0;
    double factor = base;
    std::size_t remaining = exponent;
    while (remaining > 0)
    {
        if ((remaining % 2) == 1)
        {
            result *= factor;
        }
        remaining /= 2;
        if (remaining > 0)
        {
            factor *= factor;
        }
    }
    return result;
}

void validate_reconstructed_value(const double value)
{
    if (!std::isfinite(value))
    {
        throw std::overflow_error("Raw-moment reconstruction produced a non-finite value.");
    }
}

[[nodiscard]] bool is_integer_exponent(const double exponent)
{
    return std::trunc(exponent) == exponent;
}

void validate_generalized_value(const double value)
{
    if (!std::isfinite(value))
    {
        throw std::overflow_error("Generalized moment evaluation produced a non-finite value.");
    }
}

void validate_generalized_exponent(const double exponent)
{
    if (!std::isfinite(exponent))
    {
        throw std::invalid_argument("Generalized moment exponent must be finite.");
    }
}

[[nodiscard]] double generalized_power(const double base, const double exponent)
{
    if (base < 0.0 && !is_integer_exponent(exponent))
    {
        throw std::domain_error("Generalized moment evaluation requires non-negative abscissas for non-integer exponents.");
    }

    const double value = std::pow(base, exponent);
    validate_generalized_value(value);
    return value;
}

[[nodiscard]] internal::CqmomQuadratureStorageValidationMessages cqmom_reconstruction_quadrature_messages()
{
    return {"CQMOM quadrature must contain at least one terminal node.", "CQMOM quadrature abscissas must contain one row per phase-space dimension.",
            "Every CQMOM abscissa row must contain one value per terminal node.", "CQMOM quadrature weights must be finite.",
            "CQMOM quadrature abscissas must be finite."};
}

[[nodiscard]] double reconstruct_cqmom_moment_value(const CqmomQuadrature& quadrature, const MultiIndex& moment_index)
{
    // CQMOM raw moment reconstruction uses terminal nodes:
    // M_alpha = sum_j W_j * product_d x_{d,j}^{alpha_d}.
    double value = 0.0;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        double contribution = quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            contribution *= integer_power(quadrature.abscissas()[dimension][node], moment_index[dimension]);
            validate_reconstructed_value(contribution);
        }
        value += contribution;
        validate_reconstructed_value(value);
    }
    return value;
}

[[nodiscard]] CqmomMomentData reconstruct_cqmom_moments_impl(const CqmomQuadrature& quadrature, const CqmomMomentLayout& layout,
                                                             const bool trusted_metadata)
{
    if (!trusted_metadata && (layout.hierarchy.node_counts() != quadrature.hierarchy().node_counts() ||
                              layout.hierarchy.conditioning_order() != quadrature.hierarchy().conditioning_order()))
    {
        throw std::invalid_argument("CQMOM moment layout hierarchy must match the quadrature hierarchy.");
    }
    internal::validate_cqmom_quadrature_storage(quadrature, cqmom_reconstruction_quadrature_messages());
    if (quadrature.has_rectangular_storage() && quadrature.weights().size() != layout.hierarchy.total_nodes())
    {
        throw std::invalid_argument("Fixed CQMOM quadrature node count must match the moment layout hierarchy.");
    }
    if (!quadrature.has_rectangular_storage() && quadrature.terminal_paths().size() != quadrature.weights().size())
    {
        throw std::invalid_argument("Adaptive CQMOM quadrature must provide one terminal path per terminal node.");
    }
    if (!trusted_metadata)
    {
        internal::validate_fixed_cqmom_moment_layout(layout, "CQMOM moment layout must be the fixed required layout for its hierarchy.");
    }

    const std::vector<std::size_t> moment_extents = layout.hierarchy.moment_extents();
    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const MultiIndex& index : layout.required_indices)
    {
        if (!trusted_metadata && !internal::multi_index_within_extents(index, moment_extents))
        {
            throw std::out_of_range("CQMOM moment layout contains an index outside the fixed required moment extents.");
        }
        values.push_back(reconstruct_cqmom_moment_value(quadrature, index));
    }

    return CqmomMomentData{layout, std::move(values)};
}

} // namespace

std::vector<double> reconstruct_raw_moments(const Quadrature1D& quadrature, const std::size_t moment_count)
{
    if (!quadrature.ok())
    {
        throw std::invalid_argument("Cannot reconstruct moments from an unsuccessful QMOM quadrature.");
    }
    if (moment_count == 0)
    {
        throw std::invalid_argument("Moment reconstruction requires at least one moment.");
    }
    if (quadrature.abscissas.size() != quadrature.weights.size() || quadrature.weights.empty())
    {
        throw std::invalid_argument("QMOM quadrature has inconsistent node storage.");
    }
    if (!internal::all_finite(quadrature.abscissas) || !internal::all_finite(quadrature.weights))
    {
        throw std::invalid_argument("QMOM quadrature contains non-finite values.");
    }

    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
        {
            const double contribution = quadrature.weights[node] * integer_power(quadrature.abscissas[node], order);
            validate_reconstructed_value(contribution);
            moments[order] += contribution;
            validate_reconstructed_value(moments[order]);
        }
    }
    return moments;
}

double evaluate_generalized_moment(const Quadrature1D& quadrature, const double exponent)
{
    validate_generalized_exponent(exponent);
    if (!quadrature.ok())
    {
        throw std::invalid_argument("Cannot evaluate generalized moment from an unsuccessful QMOM quadrature.");
    }
    if (quadrature.abscissas.size() != quadrature.weights.size() || quadrature.weights.empty())
    {
        throw std::invalid_argument("QMOM quadrature has inconsistent node storage.");
    }
    if (!internal::all_finite(quadrature.abscissas) || !internal::all_finite(quadrature.weights))
    {
        throw std::invalid_argument("QMOM quadrature contains non-finite values.");
    }

    double value = 0.0;
    for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
    {
        const double contribution = quadrature.weights[node] * generalized_power(quadrature.abscissas[node], exponent);
        validate_generalized_value(contribution);
        value += contribution;
        validate_generalized_value(value);
    }
    return value;
}

RawMomentTensor reconstruct_raw_moments(const CqmomQuadrature& quadrature, const std::vector<std::size_t>& extents)
{
    if (extents.size() != quadrature.number_of_dimensions())
    {
        throw std::invalid_argument("Moment reconstruction extents must match the CQMOM phase-space dimension count.");
    }
    internal::validate_cqmom_quadrature_storage(quadrature, cqmom_reconstruction_quadrature_messages());

    RawMomentTensor moments{extents, std::vector<double>(product_extents(extents), 0.0)};
    for (std::size_t flat_moment = 0; flat_moment < moments.values.size(); ++flat_moment)
    {
        const MultiIndex moment_index = unflatten_multi_index(flat_moment, extents);
        moments.values[flat_moment] = reconstruct_cqmom_moment_value(quadrature, moment_index);
    }
    return moments;
}

SelectedRawMomentResult reconstruct_raw_moments(const CqmomQuadrature& quadrature, const RawMomentRequest& request)
{
    internal::validate_cqmom_quadrature_storage(quadrature, cqmom_reconstruction_quadrature_messages());

    std::vector<double> values;
    values.reserve(request.indices.size());
    for (const MultiIndex& index : request.indices)
    {
        if (index.size() != quadrature.number_of_dimensions())
        {
            throw std::invalid_argument("Requested raw moment rank must match the CQMOM phase-space dimension count.");
        }
        values.push_back(reconstruct_cqmom_moment_value(quadrature, index));
    }

    return SelectedRawMomentResult{
        request.indices, std::move(values), RawMomentStatus::success, "", invalid_component_index, invalid_dimension_index, {}};
}

double evaluate_generalized_moment(const CqmomQuadrature& quadrature, const std::vector<double>& exponents)
{
    if (exponents.size() != quadrature.number_of_dimensions())
    {
        throw std::invalid_argument("Generalized moment exponent count must match the CQMOM phase-space dimension count.");
    }
    for (const double exponent : exponents)
    {
        validate_generalized_exponent(exponent);
    }
    internal::validate_cqmom_quadrature_storage(quadrature, cqmom_reconstruction_quadrature_messages());

    double value = 0.0;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        double contribution = quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            contribution *= generalized_power(quadrature.abscissas()[dimension][node], exponents[dimension]);
            validate_generalized_value(contribution);
        }
        value += contribution;
        validate_generalized_value(value);
    }
    return value;
}

CqmomMomentData reconstruct_cqmom_moments(const CqmomQuadrature& quadrature, const CqmomMomentLayout& layout)
{
    return reconstruct_cqmom_moments_impl(quadrature, layout, false);
}

namespace internal
{

CqmomMomentData reconstruct_trusted_cqmom_moments(const CqmomQuadrature& quadrature, const CqmomMomentLayout& layout)
{
    return reconstruct_cqmom_moments_impl(quadrature, layout, true);
}

} // namespace internal

} // namespace qmomkit
