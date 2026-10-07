#pragma once

#include "qmomkit/cqmom.hpp"

#include <optional>
#include <vector>

namespace qmomkit::internal
{

struct RecursiveQuadratureResult;

struct QmomScalingContext
{
    std::vector<double> moments;
    std::optional<AbscissaBounds> abscissa_bounds;
    double scale = 1.0;
};

struct CqmomScalingContext
{
    CqmomMomentData moments;
    CqmomAbscissaBounds coordinate_support;
    std::vector<double> scales;
};

[[nodiscard]] std::optional<QmomScalingContext> make_qmom_scaling_context(const std::vector<double>& moments,
                                                                          std::optional<AbscissaBounds> abscissa_bounds);
// The scale of each phase-space coordinate, from the marginal moments of the
// whole distribution: sqrt(m2 / m0), else |m1| / m0, else one.
[[nodiscard]] std::vector<double> cqmom_coordinate_scales(const CqmomMomentData& moments);
[[nodiscard]] std::optional<CqmomScalingContext> make_cqmom_scaling_context(const CqmomMomentData& moments,
                                                                            const CqmomAbscissaBounds& coordinate_support);

[[nodiscard]] bool restore_abscissas(std::vector<double>& abscissas, double scale);
[[nodiscard]] bool restore_cqmom_abscissas(std::vector<std::vector<double>>& abscissas, const std::vector<double>& scales);
[[nodiscard]] bool restore_bounds(AbscissaBounds& bounds, double scale);
[[nodiscard]] bool can_restore_cqmom_result(const RecursiveQuadratureResult& result, const std::vector<double>& scales);

} // namespace qmomkit::internal
