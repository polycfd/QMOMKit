#pragma once

#include "numeric_utils.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/multi_index.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace qmomkit::internal
{

struct CqmomMomentDataValidationMessages
{
    std::string fixed_layout;
    std::string value_count;
    std::string finite_values;
};

struct CqmomQuadratureStorageValidationMessages
{
    std::string empty_weights;
    std::string abscissa_rows;
    std::string abscissa_node_count;
    std::string finite_weights;
    std::string finite_abscissas;
};

struct CqmomMomentRateValidationMessages
{
    std::string value_count;
    std::string finite_values;
};

[[nodiscard]] inline bool same_cqmom_layout(const CqmomMomentLayout& lhs, const CqmomMomentLayout& rhs)
{
    return lhs.hierarchy.node_counts() == rhs.hierarchy.node_counts() && lhs.hierarchy.conditioning_order() == rhs.hierarchy.conditioning_order() &&
           lhs.required_indices == rhs.required_indices;
}

[[nodiscard]] inline bool multi_index_within_extents(const MultiIndex& index, const std::vector<std::size_t>& extents)
{
    if (index.size() != extents.size())
    {
        return false;
    }
    for (std::size_t axis = 0; axis < extents.size(); ++axis)
    {
        if (index[axis] >= extents[axis])
        {
            return false;
        }
    }
    return true;
}

inline void validate_cqmom_moment_layout_matches_expected(const CqmomMomentLayout& layout, const CqmomMomentLayout& expected_layout,
                                                          const std::string& message)
{
    if (layout.hierarchy.node_counts() != expected_layout.hierarchy.node_counts() ||
        layout.hierarchy.conditioning_order() != expected_layout.hierarchy.conditioning_order() ||
        layout.required_indices != expected_layout.required_indices)
    {
        throw std::invalid_argument(message);
    }
}

inline void validate_fixed_cqmom_moment_layout(const CqmomMomentLayout& layout, const std::string& message)
{
    const CqmomMomentLayout expected_layout = create_cqmom_moment_layout(layout.hierarchy);
    validate_cqmom_moment_layout_matches_expected(layout, expected_layout, message);
}

inline void validate_cqmom_moment_layout_matches_hierarchy(const CqmomMomentLayout& layout, const PhaseSpaceHierarchy& hierarchy,
                                                           const std::string& message)
{
    const CqmomMomentLayout expected_layout = create_cqmom_moment_layout(hierarchy);
    validate_cqmom_moment_layout_matches_expected(layout, expected_layout, message);
}

inline void validate_cqmom_moment_values_match_layout(const CqmomMomentData& moments, const CqmomMomentDataValidationMessages& messages)
{
    if (moments.values.size() != moments.layout.required_indices.size())
    {
        throw std::invalid_argument(messages.value_count);
    }
    if (!all_finite(moments.values))
    {
        throw std::invalid_argument(messages.finite_values);
    }
}

inline void validate_fixed_cqmom_moment_data(const CqmomMomentData& moments, const CqmomMomentDataValidationMessages& messages)
{
    validate_fixed_cqmom_moment_layout(moments.layout, messages.fixed_layout);
    validate_cqmom_moment_values_match_layout(moments, messages);
}

inline void validate_cqmom_moment_data_matches_hierarchy(const CqmomMomentData& moments, const PhaseSpaceHierarchy& hierarchy,
                                                         const CqmomMomentDataValidationMessages& messages)
{
    const CqmomMomentLayout expected_layout = create_cqmom_moment_layout(hierarchy);
    validate_cqmom_moment_layout_matches_expected(moments.layout, expected_layout, messages.fixed_layout);
    validate_cqmom_moment_values_match_layout(moments, messages);
}

inline void validate_dense_cqmom_moment_tensor(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const char* rank_message,
                                               const char* value_count_message, const char* finite_values_message)
{
    if (moments.extents.size() != hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument(rank_message);
    }
    if (moments.values.size() != product_extents(moments.extents))
    {
        throw std::invalid_argument(value_count_message);
    }
    if (!all_finite(moments.values))
    {
        throw std::invalid_argument(finite_values_message);
    }
}

inline void validate_dense_cqmom_reconstruction_tensor(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const char* rank_message,
                                                       const char* value_count_message, const char* finite_values_message,
                                                       const char* minimum_extents_message)
{
    validate_dense_cqmom_moment_tensor(moments, hierarchy, rank_message, value_count_message, finite_values_message);
    for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
    {
        if (moments.extents[dimension] < 2 * hierarchy.node_count(dimension))
        {
            throw std::invalid_argument(minimum_extents_message);
        }
    }
}

inline void validate_cqmom_quadrature_storage(const CqmomQuadrature& quadrature, const CqmomQuadratureStorageValidationMessages& messages)
{
    if (quadrature.weights().empty())
    {
        throw std::invalid_argument(messages.empty_weights);
    }
    if (quadrature.abscissas().size() != quadrature.number_of_dimensions())
    {
        throw std::invalid_argument(messages.abscissa_rows);
    }
    for (const std::vector<double>& dimension_abscissas : quadrature.abscissas())
    {
        if (dimension_abscissas.size() != quadrature.weights().size())
        {
            throw std::invalid_argument(messages.abscissa_node_count);
        }
    }
    if (!all_finite(quadrature.weights()))
    {
        throw std::invalid_argument(messages.finite_weights);
    }
    for (const std::vector<double>& dimension_abscissas : quadrature.abscissas())
    {
        if (!all_finite(dimension_abscissas))
        {
            throw std::invalid_argument(messages.finite_abscissas);
        }
    }
}

inline void validate_cqmom_moment_rate(const CqmomMomentRate& rate, const std::size_t expected_size,
                                       const CqmomMomentRateValidationMessages& messages)
{
    if (rate.values.size() != expected_size)
    {
        throw std::invalid_argument(messages.value_count);
    }
    if (!all_finite(rate.values))
    {
        throw std::invalid_argument(messages.finite_values);
    }
}

} // namespace qmomkit::internal
