#pragma once

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit::internal
{

class DenseMatrix
{
public:
    DenseMatrix(const std::size_t rows, const std::size_t cols) : rows_(rows), cols_(cols), values_(rows * cols, 0.0) {}

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
};

struct JacobiOptions
{
    double tolerance = 1.0e-14;
    std::size_t max_sweeps = 0;
    std::size_t initial_q = 0;
};

[[nodiscard]] double matrix_infinity_norm(const DenseMatrix& matrix);

[[nodiscard]] double off_diagonal_norm(const DenseMatrix& matrix);

[[nodiscard]] LinearSolveResult solve_dense(DenseMatrix matrix, DenseMatrix rhs, double singular_tolerance, const std::string& singular_message);

[[nodiscard]] DenseMatrix make_identity(std::size_t n);

[[nodiscard]] double residual_norm(const DenseMatrix& matrix, const DenseMatrix& solution, const DenseMatrix& rhs);

[[nodiscard]] double condition_estimate(const DenseMatrix& matrix, double singular_tolerance);

[[nodiscard]] EigenDecomposition symmetric_jacobi_eigen(DenseMatrix matrix, const JacobiOptions& options);

} // namespace qmomkit::internal
