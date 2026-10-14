#include "monte_carlo/moments.hpp"

#include "ensemble_checks.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit::monte_carlo
{
namespace
{

[[nodiscard]] bool valid_number_density(const double number_density)
{
    return std::isfinite(number_density) && number_density >= 0.0;
}

[[nodiscard]] EnsembleMomentStatus validate_ensemble(const Ensemble& ensemble, std::size_t& dimensions, std::string& message,
                                                     const bool dimensions_first = false)
{
    const internal::EnsembleCheck check = internal::check_ensemble(ensemble, dimensions_first);
    dimensions = check.dimensions;
    if (check.problem == internal::EnsembleProblem::none)
    {
        return EnsembleMomentStatus::success;
    }
    message = internal::ensemble_problem_message(check.problem);
    return check.problem == internal::EnsembleProblem::non_finite_coordinate ? EnsembleMomentStatus::non_finite_sample
                                                                             : EnsembleMomentStatus::invalid_input;
}

[[nodiscard]] double integer_power(const double base, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        result *= base;
        if (!std::isfinite(result))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    return result;
}

[[nodiscard]] double monomial_value(const Sample& sample, const qmomkit::MultiIndex& order)
{
    double value = 1.0;
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        value *= integer_power(sample[dimension], order[dimension]);
        if (!std::isfinite(value))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    return value;
}

[[nodiscard]] double standard_error(const double sum, const double sum_squares, const std::size_t sample_count, const double number_density)
{
    if (sample_count < 2)
    {
        return 0.0;
    }
    // Estimate the standard error of n0 * mean(Y), where Y is the sampled
    // monomial. The unbiased sample variance is formed before the n0 scaling.
    const double count = static_cast<double>(sample_count);
    const double raw_variance = (sum_squares - (sum * sum) / count) / (count - 1.0);
    const double variance = std::max(0.0, raw_variance);
    return number_density * std::sqrt(variance / count);
}

[[nodiscard]] std::vector<std::size_t> extents_from_maximum_orders(const std::vector<std::size_t>& maximum_orders)
{
    std::vector<std::size_t> extents(maximum_orders.size(), 0);
    for (std::size_t dimension = 0; dimension < maximum_orders.size(); ++dimension)
    {
        extents[dimension] = maximum_orders[dimension] + 1;
    }
    return extents;
}

[[nodiscard]] bool has_extent_overflow(const std::vector<std::size_t>& maximum_orders)
{
    return std::any_of(maximum_orders.begin(), maximum_orders.end(),
                       [](const std::size_t maximum_order) { return maximum_order == std::numeric_limits<std::size_t>::max(); });
}

[[nodiscard]] SelectedEnsembleMomentEstimate selected_failure(const EnsembleMomentStatus status, std::string message)
{
    return SelectedEnsembleMomentEstimate{{}, {}, {}, status, std::move(message)};
}

[[nodiscard]] GeneralizedEnsembleMomentEstimate generalized_failure(const std::vector<std::vector<double>>& orders, const EnsembleMomentStatus status,
                                                                    std::string message)
{
    return GeneralizedEnsembleMomentEstimate{orders, {}, status, std::move(message)};
}

[[nodiscard]] bool is_integer(const double value)
{
    return std::trunc(value) == value;
}

[[nodiscard]] double generalized_monomial_value(const Sample& sample, const std::vector<double>& order, bool& domain_error)
{
    double value = 1.0;
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        if (sample[dimension] < 0.0 && !is_integer(order[dimension]))
        {
            domain_error = true;
            return std::numeric_limits<double>::quiet_NaN();
        }
        value *= std::pow(sample[dimension], order[dimension]);
        if (!std::isfinite(value))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    return value;
}

} // namespace

OneDimensionalMomentEstimate compute_ensemble_moments(const Ensemble& ensemble, const double number_density, const std::size_t maximum_order)
{
    std::size_t dimensions = 0;
    std::string message;
    const EnsembleMomentStatus ensemble_status = validate_ensemble(ensemble, dimensions, message);
    if (ensemble_status != EnsembleMomentStatus::success)
    {
        return OneDimensionalMomentEstimate{{}, {}, ensemble_status, message};
    }
    if (dimensions != 1)
    {
        return OneDimensionalMomentEstimate{
            {}, {}, EnsembleMomentStatus::invalid_input, "One-dimensional moment extraction requires one-dimensional samples."};
    }
    if (!valid_number_density(number_density))
    {
        return OneDimensionalMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "Number density must be finite and non-negative."};
    }
    if (maximum_order == std::numeric_limits<std::size_t>::max())
    {
        return OneDimensionalMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "Maximum order is too large."};
    }

    std::vector<double> moments(maximum_order + 1, 0.0);
    std::vector<double> errors(maximum_order + 1, 0.0);
    const double inverse_count = 1.0 / static_cast<double>(ensemble.size());
    // Accumulate each raw power over the same finite realization:
    // m_k = n0 / N * sum_p x_p^k.
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        double sum = 0.0;
        double sum_squares = 0.0;
        const qmomkit::MultiIndex multi_index{order};
        for (const Sample& sample : ensemble)
        {
            const double value = monomial_value(sample, multi_index);
            sum += value;
            sum_squares += value * value;
            if (!std::isfinite(sum) || !std::isfinite(sum_squares))
            {
                return OneDimensionalMomentEstimate{{}, {}, EnsembleMomentStatus::non_finite_moment, "Moment accumulation became non-finite."};
            }
        }
        moments[order] = number_density * sum * inverse_count;
        errors[order] = standard_error(sum, sum_squares, ensemble.size(), number_density);
        if (!std::isfinite(moments[order]) || !std::isfinite(errors[order]))
        {
            return OneDimensionalMomentEstimate{{}, {}, EnsembleMomentStatus::non_finite_moment, "Moment or standard error became non-finite."};
        }
    }

    return OneDimensionalMomentEstimate{qmomkit::MomentSet(std::move(moments)), std::move(errors), EnsembleMomentStatus::success, ""};
}

TensorMomentEstimate compute_ensemble_moments(const Ensemble& ensemble, const double number_density, const std::vector<std::size_t>& maximum_orders)
{
    std::size_t dimensions = 0;
    std::string message;
    const EnsembleMomentStatus ensemble_status = validate_ensemble(ensemble, dimensions, message);
    if (ensemble_status != EnsembleMomentStatus::success)
    {
        return TensorMomentEstimate{{}, {}, ensemble_status, message};
    }
    if (maximum_orders.size() != dimensions)
    {
        return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "Maximum-order rank must match the sample phase-space dimension."};
    }
    if (!valid_number_density(number_density))
    {
        return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "Number density must be finite and non-negative."};
    }
    if (has_extent_overflow(maximum_orders))
    {
        return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, "At least one maximum order is too large."};
    }

    const std::vector<std::size_t> extents = extents_from_maximum_orders(maximum_orders);
    std::size_t moment_count = 0;
    try
    {
        moment_count = qmomkit::product_extents(extents);
    }
    catch (const std::exception& error)
    {
        return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::invalid_input, std::string("Tensor moment extents are invalid: ") + error.what()};
    }
    std::vector<double> values(moment_count, 0.0);
    std::vector<double> errors(moment_count, 0.0);
    const double inverse_count = 1.0 / static_cast<double>(ensemble.size());
    // Multi-index enumeration follows RawMomentTensor phase-space order; this
    // keeps sampled moments directly comparable with CQMOM tensor moments.
    const std::vector<qmomkit::MultiIndex> orders = qmomkit::enumerate_multi_indices(extents);

    // Each tensor entry estimates M_alpha = n0 / N * sum_p prod_d xi[p,d]^alpha[d].
    for (std::size_t flat = 0; flat < moment_count; ++flat)
    {
        double sum = 0.0;
        double sum_squares = 0.0;
        for (const Sample& sample : ensemble)
        {
            const double value = monomial_value(sample, orders[flat]);
            sum += value;
            sum_squares += value * value;
            if (!std::isfinite(sum) || !std::isfinite(sum_squares))
            {
                return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::non_finite_moment, "Tensor moment accumulation became non-finite."};
            }
        }
        values[flat] = number_density * sum * inverse_count;
        errors[flat] = standard_error(sum, sum_squares, ensemble.size(), number_density);
        if (!std::isfinite(values[flat]) || !std::isfinite(errors[flat]))
        {
            return TensorMomentEstimate{{}, {}, EnsembleMomentStatus::non_finite_moment, "Tensor moment or standard error became non-finite."};
        }
    }

    return TensorMomentEstimate{qmomkit::RawMomentTensor{extents, std::move(values)}, std::move(errors), EnsembleMomentStatus::success, ""};
}

SelectedEnsembleMomentEstimate compute_requested_ensemble_moments(const Ensemble& ensemble, const double number_density,
                                                                  const qmomkit::RawMomentRequest& request)
{
    std::size_t dimensions = 0;
    std::string message;
    const EnsembleMomentStatus ensemble_status = validate_ensemble(ensemble, dimensions, message, true);
    if (ensemble_status != EnsembleMomentStatus::success)
    {
        return selected_failure(ensemble_status, std::move(message));
    }

    for (const qmomkit::MultiIndex& index : request.indices)
    {
        if (index.size() != dimensions)
        {
            return selected_failure(EnsembleMomentStatus::invalid_input, "Requested raw moment rank must match the sample phase-space dimension.");
        }
    }

    if (!valid_number_density(number_density))
    {
        return selected_failure(EnsembleMomentStatus::invalid_input, "Number density must be finite and non-negative.");
    }

    for (const qmomkit::MultiIndex& index : request.indices)
    {
        for (const std::size_t exponent : index)
        {
            if (exponent == std::numeric_limits<std::size_t>::max())
            {
                return selected_failure(EnsembleMomentStatus::invalid_input, "Requested raw moment exponent is too large.");
            }
        }
    }

    if (request.indices.empty())
    {
        return SelectedEnsembleMomentEstimate{{}, {}, {}, EnsembleMomentStatus::success, ""};
    }

    std::vector<double> values(request.indices.size(), 0.0);
    std::vector<double> errors(request.indices.size(), 0.0);
    const double inverse_count = 1.0 / static_cast<double>(ensemble.size());
    for (std::size_t requested = 0; requested < request.indices.size(); ++requested)
    {
        double sum = 0.0;
        double sum_squares = 0.0;
        for (const Sample& sample : ensemble)
        {
            const double value = monomial_value(sample, request.indices[requested]);
            if (!std::isfinite(value))
            {
                return selected_failure(EnsembleMomentStatus::non_finite_moment, "Selected moment accumulation became non-finite.");
            }
            sum += value;
            sum_squares += value * value;
            if (!std::isfinite(sum) || !std::isfinite(sum_squares))
            {
                return selected_failure(EnsembleMomentStatus::non_finite_moment, "Selected moment accumulation became non-finite.");
            }
        }

        values[requested] = number_density * sum * inverse_count;
        const double count = static_cast<double>(ensemble.size());
        const double raw_variance = ensemble.size() < 2 ? 0.0 : (sum_squares - (sum * sum) / count) / (count - 1.0);
        errors[requested] = standard_error(sum, sum_squares, ensemble.size(), number_density);
        if (!std::isfinite(values[requested]) || !std::isfinite(raw_variance) || !std::isfinite(errors[requested]))
        {
            return selected_failure(EnsembleMomentStatus::non_finite_moment, "Selected moment or standard error became non-finite.");
        }
    }

    return SelectedEnsembleMomentEstimate{request.indices, std::move(values), std::move(errors), EnsembleMomentStatus::success, ""};
}

GeneralizedEnsembleMomentEstimate compute_generalized_ensemble_moments(const Ensemble& ensemble, const double number_density,
                                                                       const std::vector<std::vector<double>>& orders)
{
    std::size_t dimensions = 0;
    std::string message;
    const EnsembleMomentStatus ensemble_status = validate_ensemble(ensemble, dimensions, message);
    if (ensemble_status != EnsembleMomentStatus::success)
    {
        return generalized_failure(orders, ensemble_status, std::move(message));
    }
    if (!valid_number_density(number_density))
    {
        return generalized_failure(orders, EnsembleMomentStatus::invalid_input, "Number density must be finite and non-negative.");
    }
    for (const std::vector<double>& order : orders)
    {
        if (order.size() != dimensions)
        {
            return generalized_failure(orders, EnsembleMomentStatus::invalid_input,
                                       "Generalized moment rank must match the sample phase-space dimension.");
        }
        if (!std::all_of(order.begin(), order.end(), [](const double exponent) { return std::isfinite(exponent) && exponent >= 0.0; }))
        {
            return generalized_failure(orders, EnsembleMomentStatus::invalid_input, "Generalized moment exponents must be finite and non-negative.");
        }
    }

    std::vector<double> values(orders.size(), 0.0);
    const double inverse_count = 1.0 / static_cast<double>(ensemble.size());
    for (std::size_t requested = 0; requested < orders.size(); ++requested)
    {
        double sum = 0.0;
        for (const Sample& sample : ensemble)
        {
            bool domain_error = false;
            const double value = generalized_monomial_value(sample, orders[requested], domain_error);
            if (domain_error)
            {
                return generalized_failure(orders, EnsembleMomentStatus::invalid_input,
                                           "Generalized moment evaluation requires non-negative coordinates for non-integer exponents.");
            }
            if (!std::isfinite(value))
            {
                return generalized_failure(orders, EnsembleMomentStatus::non_finite_moment, "Generalized moment evaluation became non-finite.");
            }
            sum += value;
            if (!std::isfinite(sum))
            {
                return generalized_failure(orders, EnsembleMomentStatus::non_finite_moment, "Generalized moment accumulation became non-finite.");
            }
        }
        values[requested] = number_density * sum * inverse_count;
        if (!std::isfinite(values[requested]))
        {
            return generalized_failure(orders, EnsembleMomentStatus::non_finite_moment, "Generalized moment evaluation became non-finite.");
        }
    }
    return GeneralizedEnsembleMomentEstimate{orders, std::move(values), EnsembleMomentStatus::success, ""};
}

} // namespace qmomkit::monte_carlo
