#pragma once

#include "numeric_utils.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace qmomkit::internal
{

class QmomStateAccess
{
public:
    [[nodiscard]] static QmomState from_quadrature(QmomResult quadrature, const std::size_t moment_count)
    {
        if (moment_count == 0)
        {
            throw std::invalid_argument("QMOM state requires at least one stored moment.");
        }
        if (!quadrature.ok())
        {
            throw std::invalid_argument("QMOM state quadrature must be successful.");
        }
        if (quadrature.weights.empty())
        {
            throw std::invalid_argument("QMOM state quadrature must contain at least one node.");
        }
        if (quadrature.abscissas.size() != quadrature.weights.size())
        {
            throw std::invalid_argument("QMOM state quadrature has inconsistent node storage.");
        }
        if (!all_finite(quadrature.abscissas) || !all_finite(quadrature.weights))
        {
            throw std::invalid_argument("QMOM state quadrature must be finite.");
        }

        const std::vector<double> moments =
            reconstruct_raw_moments(Quadrature1D{quadrature.abscissas, quadrature.weights, QuadratureStatus::success, ""}, moment_count);
        return QmomState{std::move(moments), std::move(quadrature)};
    }
};

} // namespace qmomkit::internal
