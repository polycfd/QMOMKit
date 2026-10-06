#pragma once

#include "checked_size.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace qmomkit::internal
{

class DenseMatrix
{
public:
    DenseMatrix(const std::size_t rows, const std::size_t cols) : rows_(rows), cols_(cols), values_(storage_size(rows, cols), 0.0) {}

    [[nodiscard]] std::size_t rows() const noexcept
    {
        return rows_;
    }

    [[nodiscard]] std::size_t cols() const noexcept
    {
        return cols_;
    }

    double& operator()(const std::size_t row, const std::size_t col)
    {
        assert(row < rows_);
        assert(col < cols_);
        const std::size_t index = row * cols_ + col;
        assert(index < values_.size());
        return values_[index];
    }

    double operator()(const std::size_t row, const std::size_t col) const
    {
        assert(row < rows_);
        assert(col < cols_);
        const std::size_t index = row * cols_ + col;
        assert(index < values_.size());
        return values_[index];
    }

private:
    [[nodiscard]] static std::size_t storage_size(const std::size_t rows, const std::size_t cols)
    {
        const std::optional<std::size_t> size = checked_size_multiply(rows, cols);
        if (!size.has_value())
        {
            throw std::length_error("DenseMatrix storage dimensions overflow size_t.");
        }
        return *size;
    }

    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> values_;
};

struct LinearSolveResult
{
    DenseMatrix solution;
    bool success = true;
    std::string message;
};

struct EigenDecomposition
{
    std::vector<double> eigenvalues;
    DenseMatrix eigenvectors;
    // False when the rotation budget ran out before every off-diagonal entry
    // met the tolerance; the eigenpairs are then only approximate.
    bool converged = true;
};

// The solver first scales the matrix by the power of two that places its
// largest entry in [1, 2), so `tolerance` bounds the off-diagonal entries
// relative to the largest entry and the result does not depend on the
// matrix's overall magnitude. The default stops at machine precision, which
// the quadratic convergence of the final rotations reaches at almost no extra
// cost. `max_rotations` limits the number of single pair rotations; zero
// selects 10 n^2, five times what classical Jacobi was measured to need
// (under 2 n^2 for tridiagonal Jacobi and dense Hankel matrices).
struct JacobiOptions
{
    double tolerance = std::numeric_limits<double>::epsilon();
    std::size_t max_rotations = 0;
};

[[nodiscard]] double matrix_infinity_norm(const DenseMatrix& matrix);

[[nodiscard]] LinearSolveResult solve_dense(DenseMatrix matrix, DenseMatrix rhs, double singular_tolerance, const std::string& singular_message);

[[nodiscard]] DenseMatrix make_identity(std::size_t n);

[[nodiscard]] double residual_norm(const DenseMatrix& matrix, const DenseMatrix& solution, const DenseMatrix& rhs);

[[nodiscard]] double condition_estimate(const DenseMatrix& matrix, double singular_tolerance);

[[nodiscard]] EigenDecomposition symmetric_jacobi_eigen(DenseMatrix matrix, const JacobiOptions& options);

} // namespace qmomkit::internal
