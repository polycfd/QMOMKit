#pragma once

#include "qmomkit/moment_source.hpp"

namespace gqmom_1d_absolute_value_loss
{

/**
 * The loss of elements at the rate c |xi| as a moment source written over the
 * nodes of the source-stage quadrature. xi is dimensionless, and c = 1.
 */
[[nodiscard]] qmomkit::QmomMomentSource absolute_value_loss();

struct ExactMoments
{
    double m0 = 0.0;
    double m2 = 0.0;
};

/**
 * The moments m0 and m2 of the exact solution n(t, xi) = n(0, xi) exp(-c |xi| t)
 * at time t, for a normal n(0, xi) with the given number density, mean and
 * standard deviation.
 */
[[nodiscard]] ExactMoments exact_moments(double number_density, double mean, double standard_deviation, double time);

} // namespace gqmom_1d_absolute_value_loss
