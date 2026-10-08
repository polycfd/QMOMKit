#pragma once

namespace qmomkit
{

/**
 * Closed interval used to accept reconstructed one-dimensional quadrature
 * abscissas, up to a margin of 1e-12 times max(|lower|, |upper|) that absorbs
 * rounding. Abscissas are never clipped to it. This is reconstruction metadata
 * and is independent of moment realizability support checks.
 */
struct AbscissaBounds
{
    double lower = 0.0;
    double upper = 0.0;
};

} // namespace qmomkit
