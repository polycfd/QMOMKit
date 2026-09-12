#pragma once

namespace qmomkit
{

/**
 * Closed interval used to accept reconstructed one-dimensional quadrature
 * abscissas. This is reconstruction metadata and is independent of moment
 * realizability support checks.
 */
struct AbscissaBounds
{
    double lower = 0.0;
    double upper = 0.0;
};

} // namespace qmomkit
