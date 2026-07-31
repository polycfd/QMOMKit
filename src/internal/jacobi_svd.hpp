#pragma once

#include "linear_algebra.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit::internal
{

enum class SvdStatus
{
    success,
    invalid_input,
    rank_deficient,
    no_convergence
};

struct SvdOptions
{
    double relative_rank_tolerance = 1.0e-14;
    double jacobi_tolerance = 0.0;
    std::size_t max_sweeps = 0;
};

struct SvdFactorization
{
    DenseMatrix left_singular_vectors{0, 0};
    std::vector<double> singular_values;
    DenseMatrix right_singular_vectors{0, 0};

    std::size_t rows = 0;
    std::size_t cols = 0;
    std::size_t numerical_rank = 0;
    std::size_t sweeps = 0;

    double relative_rank_tolerance = 1.0e-14;
    double jacobi_tolerance = 0.0;

    SvdStatus status = SvdStatus::success;
    std::string message;
};

struct SvdSolveResult
{
    DenseMatrix solution{0, 0};
    std::size_t numerical_rank = 0;
    SvdStatus status = SvdStatus::success;
    std::string message;
};

/**
 * Factorize a square or tall dense matrix as A = U Sigma V^T.
 *
 * The implementation uses a dependency-free one-sided Jacobi iteration. It
 * stores thin left singular vectors U, descending singular values, and full
 * right singular vectors V. Square full-rank factorizations support ordinary
 * linear solves; tall full-column-rank factorizations support least-squares
 * solves. Rank-deficient and non-converged factorizations remain inspectable,
 * but solve_svd refuses them and does not produce pseudoinverse, truncated,
 * regularized, or minimum-norm solutions.
 */
[[nodiscard]] SvdFactorization factorize_svd(DenseMatrix matrix, const SvdOptions& options = {});

/**
 * Reuse a successful full-column-rank SVD for one or more right-hand sides.
 */
[[nodiscard]] SvdSolveResult solve_svd(const SvdFactorization& factorization, DenseMatrix rhs);

[[nodiscard]] SvdSolveResult solve_svd(DenseMatrix matrix, DenseMatrix rhs, const SvdOptions& options = {});

} // namespace qmomkit::internal
