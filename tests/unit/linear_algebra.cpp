#include "internal/linear_algebra.hpp"

#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace
{

using qmomkit::internal::DenseMatrix;
using qmomkit::internal::EigenDecomposition;

constexpr double epsilon = std::numeric_limits<double>::epsilon();

DenseMatrix make_matrix(const std::vector<std::vector<double>>& rows)
{
    DenseMatrix matrix(rows.size(), rows.size());
    for (std::size_t row = 0; row < rows.size(); ++row)
    {
        assert(rows[row].size() == rows.size());
        for (std::size_t col = 0; col < rows.size(); ++col)
        {
            matrix(row, col) = rows[row][col];
        }
    }
    return matrix;
}

// The options the Wheeler inversion passes to the eigen solver.
EigenDecomposition eigen(const DenseMatrix& matrix)
{
    return qmomkit::internal::symmetric_jacobi_eigen(matrix, {1.0e-14, 50 * std::max<std::size_t>(matrix.rows(), 1)});
}

// Jacobi matrix of the Legendre polynomials: zero diagonal and off-diagonal
// sqrt(beta_k) = k / sqrt(4k^2 - 1), from the monic recurrence coefficient
// beta_k = k^2 / (4k^2 - 1). Its eigenvalues are the Gauss-Legendre nodes on
// [-1, 1], and 2 v_0i^2 are the weights, because the Legendre weight function
// has total mass 2.
DenseMatrix legendre_jacobi_matrix(const std::size_t n)
{
    DenseMatrix matrix(n, n);
    for (std::size_t k = 1; k < n; ++k)
    {
        const double kk = static_cast<double>(k);
        const double off_diagonal = kk / std::sqrt(4.0 * kk * kk - 1.0);
        matrix(k - 1, k) = off_diagonal;
        matrix(k, k - 1) = off_diagonal;
    }
    return matrix;
}

std::vector<double> legendre_weights(const EigenDecomposition& decomposition)
{
    std::vector<double> weights;
    for (std::size_t node = 0; node < decomposition.eigenvalues.size(); ++node)
    {
        weights.push_back(2.0 * decomposition.eigenvectors(0, node) * decomposition.eigenvectors(0, node));
    }
    return weights;
}

double infinity_norm(const DenseMatrix& matrix)
{
    return qmomkit::internal::matrix_infinity_norm(matrix);
}

// max_j || A v_j - lambda_j v_j ||_inf: every column must be an eigenvector
// of the eigenvalue at the same position.
double eigen_residual(const DenseMatrix& matrix, const EigenDecomposition& decomposition)
{
    const std::size_t n = matrix.rows();
    double residual = 0.0;
    for (std::size_t col = 0; col < n; ++col)
    {
        for (std::size_t row = 0; row < n; ++row)
        {
            double product = 0.0;
            for (std::size_t k = 0; k < n; ++k)
            {
                product += matrix(row, k) * decomposition.eigenvectors(k, col);
            }
            residual = std::max(residual, std::abs(product - decomposition.eigenvalues[col] * decomposition.eigenvectors(row, col)));
        }
    }
    return residual;
}

// max_ij |(V^T V - I)_ij|: the eigenvectors must be orthonormal.
double orthogonality_error(const EigenDecomposition& decomposition)
{
    const std::size_t n = decomposition.eigenvalues.size();
    double error = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            double product = 0.0;
            for (std::size_t k = 0; k < n; ++k)
            {
                product += decomposition.eigenvectors(k, i) * decomposition.eigenvectors(k, j);
            }
            error = std::max(error, std::abs(product - (i == j ? 1.0 : 0.0)));
        }
    }
    return error;
}

// A backward-stable symmetric eigen solver leaves a residual and an
// orthogonality error of order n * epsilon relative to the matrix norm; the
// factor 16 leaves room for the accumulated rotations.
void assert_valid_decomposition(const DenseMatrix& matrix, const EigenDecomposition& decomposition)
{
    const std::size_t n = matrix.rows();
    const double tolerance = 16.0 * static_cast<double>(n) * epsilon;
    assert(decomposition.eigenvalues.size() == n);
    assert(decomposition.eigenvectors.rows() == n);
    assert(decomposition.eigenvectors.cols() == n);
    assert(std::is_sorted(decomposition.eigenvalues.begin(), decomposition.eigenvalues.end()));
    assert(eigen_residual(matrix, decomposition) <= tolerance * infinity_norm(matrix));
    assert(orthogonality_error(decomposition) <= tolerance);
}

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

// A 2 x 2 symmetric matrix [[a, b], [b, c]] has the eigenvalues
// (a + c) / 2 -/+ sqrt(((a - c) / 2)^2 + b^2). The solver diagonalizes it with
// one rotation, so the eigenvalues are accurate to a few rounding errors of
// the matrix norm.
void test_two_by_two_matrices_match_the_closed_form()
{
    const std::vector<std::vector<double>> entries{{2.0, 1.0, 3.0}, {1.0, 0.5, -4.0}, {0.0, 1.0, 0.0}, {3.0, 2.0, 3.0}, {-1.0e3, 7.0, 2.5e3}};
    for (const std::vector<double>& abc : entries)
    {
        const double a = abc[0];
        const double b = abc[1];
        const double c = abc[2];
        const DenseMatrix matrix = make_matrix({{a, b}, {b, c}});
        const EigenDecomposition decomposition = eigen(matrix);

        const double centre = 0.5 * (a + c);
        const double radius = std::hypot(0.5 * (a - c), b);
        const double tolerance = 8.0 * epsilon * infinity_norm(matrix);
        assert_near(decomposition.eigenvalues[0], centre - radius, tolerance);
        assert_near(decomposition.eigenvalues[1], centre + radius, tolerance);
        assert_valid_decomposition(matrix, decomposition);
    }
}

// A diagonal matrix is already diagonal: the eigenvalues are its diagonal in
// ascending order and the eigenvectors are the matching unit vectors, exactly.
void test_diagonal_matrix_is_returned_sorted_and_exact()
{
    const std::vector<double> diagonal{3.0, -1.0, 2.0, 0.5};
    DenseMatrix matrix(diagonal.size(), diagonal.size());
    for (std::size_t i = 0; i < diagonal.size(); ++i)
    {
        matrix(i, i) = diagonal[i];
    }
    const EigenDecomposition decomposition = eigen(matrix);

    assert(decomposition.eigenvalues == (std::vector<double>{-1.0, 0.5, 2.0, 3.0}));
    const std::vector<std::size_t> original_position{1, 3, 2, 0};
    for (std::size_t col = 0; col < diagonal.size(); ++col)
    {
        for (std::size_t row = 0; row < diagonal.size(); ++row)
        {
            assert(decomposition.eigenvectors(row, col) == (row == original_position[col] ? 1.0 : 0.0));
        }
    }
}

// [[2, 1, 1], [1, 2, 1], [1, 1, 2]] = I + J has the eigenvalues 1, 1, 4. A
// repeated eigenvalue has no unique eigenvectors, so the residual and the
// orthogonality are checked instead.
void test_repeated_eigenvalue_gives_an_orthonormal_basis()
{
    const DenseMatrix matrix = make_matrix({{2.0, 1.0, 1.0}, {1.0, 2.0, 1.0}, {1.0, 1.0, 2.0}});
    const EigenDecomposition decomposition = eigen(matrix);

    const double tolerance = 16.0 * epsilon * infinity_norm(matrix);
    assert_near(decomposition.eigenvalues[0], 1.0, tolerance);
    assert_near(decomposition.eigenvalues[1], 1.0, tolerance);
    assert_near(decomposition.eigenvalues[2], 4.0, tolerance);
    assert_valid_decomposition(matrix, decomposition);
}

// Gauss-Legendre nodes for two to five nodes are the roots of P_n, which have
// closed forms up to P_5, and the weights are 2 / ((1 - x^2) P_n'(x)^2). Both
// must be recovered from the Jacobi matrix of the Legendre polynomials.
void test_gauss_legendre_nodes_and_weights_match_the_closed_form()
{
    const double s65 = std::sqrt(6.0 / 5.0);
    const double s107 = std::sqrt(10.0 / 7.0);
    const double s30 = std::sqrt(30.0);
    const double s70 = std::sqrt(70.0);
    const std::vector<std::vector<double>> nodes{{-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)},
                                                 {-std::sqrt(0.6), 0.0, std::sqrt(0.6)},
                                                 {-std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * s65), -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * s65),
                                                  std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * s65), std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * s65)},
                                                 {-std::sqrt(5.0 + 2.0 * s107) / 3.0, -std::sqrt(5.0 - 2.0 * s107) / 3.0, 0.0,
                                                  std::sqrt(5.0 - 2.0 * s107) / 3.0, std::sqrt(5.0 + 2.0 * s107) / 3.0}};
    const std::vector<std::vector<double>> weights{
        {1.0, 1.0},
        {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0},
        {(18.0 - s30) / 36.0, (18.0 + s30) / 36.0, (18.0 + s30) / 36.0, (18.0 - s30) / 36.0},
        {(322.0 - 13.0 * s70) / 900.0, (322.0 + 13.0 * s70) / 900.0, 128.0 / 225.0, (322.0 + 13.0 * s70) / 900.0, (322.0 - 13.0 * s70) / 900.0}};

    for (std::size_t index = 0; index < nodes.size(); ++index)
    {
        const std::size_t n = nodes[index].size();
        const DenseMatrix matrix = legendre_jacobi_matrix(n);
        const EigenDecomposition decomposition = eigen(matrix);
        const std::vector<double> computed_weights = legendre_weights(decomposition);

        for (std::size_t node = 0; node < n; ++node)
        {
            assert_near(decomposition.eigenvalues[node], nodes[index][node], 1.0e-14);
            assert_near(computed_weights[node], weights[index][node], 1.0e-14);
        }
        assert_valid_decomposition(matrix, decomposition);
    }
}

// An n-node Gauss-Legendre rule integrates every polynomial of degree up to
// 2n - 1 exactly: sum_i w_i x_i^k = 2 / (k + 1) for even k and 0 for odd k.
// This checks the nodes and the first eigenvector row together, for up to
// sixteen nodes, against that exact integral.
void test_gauss_legendre_rules_integrate_polynomials_exactly()
{
    for (std::size_t n = 2; n <= 16; ++n)
    {
        const DenseMatrix matrix = legendre_jacobi_matrix(n);
        const EigenDecomposition decomposition = eigen(matrix);
        const std::vector<double> weights = legendre_weights(decomposition);

        for (std::size_t k = 0; k <= 2 * n - 1; ++k)
        {
            double integral = 0.0;
            for (std::size_t node = 0; node < n; ++node)
            {
                integral += weights[node] * std::pow(decomposition.eigenvalues[node], static_cast<double>(k));
            }
            const double exact = k % 2 == 0 ? 2.0 / static_cast<double>(k + 1) : 0.0;
            // The weighted sum of 2n bounded terms carries a few rounding
            // errors per term; measured errors stay below 1e-14.
            assert_near(integral, exact, 1.0e-13);
        }
        // The nodes are symmetric about zero.
        for (std::size_t node = 0; node < n; ++node)
        {
            assert_near(decomposition.eigenvalues[node], -decomposition.eigenvalues[n - 1 - node], 1.0e-14);
        }
        assert_valid_decomposition(matrix, decomposition);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"two_by_two_matrices_match_the_closed_form", test_two_by_two_matrices_match_the_closed_form},
            {"diagonal_matrix_is_returned_sorted_and_exact", test_diagonal_matrix_is_returned_sorted_and_exact},
            {"repeated_eigenvalue_gives_an_orthonormal_basis", test_repeated_eigenvalue_gives_an_orthonormal_basis},
            {"gauss_legendre_nodes_and_weights_match_the_closed_form", test_gauss_legendre_nodes_and_weights_match_the_closed_form},
            {"gauss_legendre_rules_integrate_polynomials_exactly", test_gauss_legendre_rules_integrate_polynomials_exactly},
        });
}
