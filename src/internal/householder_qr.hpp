#pragma once

#include "linear_algebra.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit::internal
{

enum class QrStatus
{
    success,
    invalid_input,
    rank_deficient
};

struct PivotedQrFactorization
{
    // Pivoted Householder QR for square or tall matrices (m >= n), stored as
    // AP = QR with column_permutation[pivoted_column] = original_column. The
    // upper triangle stores R; entries below each diagonal store Householder
    // vector tails with implicit leading one and one coefficient per reflector.
    DenseMatrix factors{0, 0};
    std::vector<double> householder_coefficients;
    std::vector<std::size_t> column_permutation;

    std::size_t rows = 0;
    std::size_t cols = 0;
    std::size_t numerical_rank = 0;

    double relative_rank_tolerance = 1.0e-14;

    QrStatus status = QrStatus::success;
    std::string message;
};

struct QrSolveResult
{
    DenseMatrix solution{0, 0};
    std::size_t numerical_rank = 0;
    QrStatus status = QrStatus::success;
    std::string message;
};

/**
 * Factorize a square or tall dense matrix as AP = QR.
 *
 * The numerical rank uses the threshold
 * relative_rank_tolerance * abs(R(0, 0)). A rank_deficient result is
 * inspectable, but solve_least_squares refuses rank-deficient factorizations:
 * it does not produce truncated, pseudoinverse, regularized, or minimum-norm
 * solutions.
 */
[[nodiscard]] PivotedQrFactorization factorize_pivoted_householder_qr(DenseMatrix matrix, double relative_rank_tolerance = 1.0e-14);

/**
 * Reuse a full-column-rank factorization for one or more right-hand sides.
 *
 * For tall matrices this returns the least-squares solution without forming
 * explicit Q or normal equations. Rank-deficient or internally inconsistent
 * factorizations are rejected before triangular back-substitution.
 */
[[nodiscard]] QrSolveResult solve_least_squares(const PivotedQrFactorization& factorization, DenseMatrix rhs);

[[nodiscard]] QrSolveResult solve_pivoted_householder_qr(DenseMatrix matrix, DenseMatrix rhs, double relative_rank_tolerance = 1.0e-14);

} // namespace qmomkit::internal
