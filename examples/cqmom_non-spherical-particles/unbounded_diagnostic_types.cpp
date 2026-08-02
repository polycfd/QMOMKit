#include "unbounded_diagnostic_types.hpp"

namespace cqmom_non_spherical_particles
{

Matrix::Matrix(const std::size_t row_count, const std::size_t col_count) : rows(row_count), cols(col_count), values(row_count * col_count, 0.0) {}

double& Matrix::operator()(const std::size_t row, const std::size_t col)
{
    return values.at(row * cols + col);
}

double Matrix::operator()(const std::size_t row, const std::size_t col) const
{
    return values.at(row * cols + col);
}

} // namespace cqmom_non_spherical_particles
