#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <qmomkit/moment_rates.hpp>
#include <qmomkit/phase_space_evolution.hpp>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr std::size_t basis_size = 3;
constexpr std::size_t quadrature_count = 3;

using Matrix3 = std::array<std::array<double, basis_size>, basis_size>;

// Keep the element deliberately small: it supplies just enough Galerkin
// machinery to show the handoff to QMOMKit without becoming an FE framework.
struct ReferenceElement
{
    std::array<double, quadrature_count> points;
    std::array<double, quadrature_count> weights;
    Matrix3 basis_values;
    Matrix3 inverse_mass;
};

[[nodiscard]] std::array<double, basis_size> quadratic_lagrange_basis(const double x)
{
    return {0.5 * x * (x - 1.0), 1.0 - x * x, 0.5 * x * (x + 1.0)};
}

[[nodiscard]] Matrix3 invert_3x3(const Matrix3& matrix)
{
    const double determinant = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
                               matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
                               matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    if (std::abs(determinant) <= 1.0e-14)
    {
        throw std::runtime_error("Reference-element mass matrix is singular.");
    }

    Matrix3 inverse{};
    inverse[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant;
    inverse[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant;
    inverse[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant;
    inverse[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant;
    inverse[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant;
    inverse[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant;
    inverse[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant;
    inverse[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant;
    inverse[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant;
    return inverse;
}

[[nodiscard]] ReferenceElement make_quadratic_reference_element()
{
    const double point = std::sqrt(3.0 / 5.0);
    ReferenceElement element{{-point, 0.0, point}, {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}, {}, {}};

    Matrix3 mass{};
    for (std::size_t q = 0; q < quadrature_count; ++q)
    {
        element.basis_values[q] = quadratic_lagrange_basis(element.points[q]);
        for (std::size_t row = 0; row < basis_size; ++row)
        {
            for (std::size_t column = 0; column < basis_size; ++column)
            {
                mass[row][column] += element.weights[q] * element.basis_values[q][row] * element.basis_values[q][column];
            }
        }
    }
    element.inverse_mass = invert_3x3(mass);
    return element;
}

[[nodiscard]] std::size_t coefficient_index(const std::size_t basis, const std::size_t moment, const std::size_t moment_count)
{
    return basis * moment_count + moment;
}

// QMOMKit works with one local moment set at a time. The Galerkin solution is
// stored as basis coefficients, so first evaluate those fields at the spatial
// quadrature points.
[[nodiscard]] std::vector<double> interpolate_coefficients_to_local_moments(const ReferenceElement& element, const std::vector<double>& coefficients,
                                                                            const std::size_t moment_count)
{
    if (coefficients.size() != basis_size * moment_count)
    {
        throw std::invalid_argument("Coefficient vector size does not match the reference element.");
    }

    std::vector<double> local_moments(quadrature_count * moment_count, 0.0);
    for (std::size_t q = 0; q < quadrature_count; ++q)
    {
        for (std::size_t basis = 0; basis < basis_size; ++basis)
        {
            for (std::size_t moment = 0; moment < moment_count; ++moment)
            {
                local_moments[q * moment_count + moment] +=
                    element.basis_values[q][basis] * coefficients[coefficient_index(basis, moment, moment_count)];
            }
        }
    }
    return local_moments;
}

// After a local QMOMKit call, the values live at the spatial Gauss points. The
// weighted projection forms the right-hand side of the element mass-matrix
// system; M^{-1} turns that into basis coefficients. The values are rates in
// Route A and evolved moments in Route B.
[[nodiscard]] std::vector<double> project_local_values_to_coefficients(const ReferenceElement& element, const std::vector<double>& local_values,
                                                                       const std::size_t moment_count)
{
    if (local_values.size() != quadrature_count * moment_count)
    {
        throw std::invalid_argument("Local value vector size does not match the reference element.");
    }

    std::vector<double> projected(basis_size * moment_count, 0.0);
    for (std::size_t q = 0; q < quadrature_count; ++q)
    {
        for (std::size_t basis = 0; basis < basis_size; ++basis)
        {
            const double weighted_basis = element.weights[q] * element.basis_values[q][basis];
            for (std::size_t moment = 0; moment < moment_count; ++moment)
            {
                projected[coefficient_index(basis, moment, moment_count)] += weighted_basis * local_values[q * moment_count + moment];
            }
        }
    }

    std::vector<double> coefficients(basis_size * moment_count, 0.0);
    for (std::size_t moment = 0; moment < moment_count; ++moment)
    {
        for (std::size_t row = 0; row < basis_size; ++row)
        {
            for (std::size_t column = 0; column < basis_size; ++column)
            {
                coefficients[coefficient_index(row, moment, moment_count)] +=
                    element.inverse_mass[row][column] * projected[coefficient_index(column, moment, moment_count)];
            }
        }
    }
    return coefficients;
}

[[nodiscard]] std::vector<double> local_moment_slice(const std::vector<double>& local_moments, const std::size_t quadrature_point,
                                                     const std::size_t moment_count)
{
    const auto first = local_moments.begin() + static_cast<std::ptrdiff_t>(quadrature_point * moment_count);
    return {first, first + static_cast<std::ptrdiff_t>(moment_count)};
}

[[nodiscard]] std::vector<double> make_spatially_varying_coefficients(const std::vector<double>& moments)
{
    // The scale field 1 + 0.15 x stays positive over the element. It therefore
    // preserves realizability while making the moment field visibly spatially
    // nonconstant.
    const std::array<double, basis_size> nodal_scales{0.85, 1.0, 1.15};
    std::vector<double> coefficients(basis_size * moments.size(), 0.0);
    for (std::size_t basis = 0; basis < basis_size; ++basis)
    {
        for (std::size_t moment = 0; moment < moments.size(); ++moment)
        {
            coefficients[coefficient_index(basis, moment, moments.size())] = nodal_scales[basis] * moments[moment];
        }
    }
    return coefficients;
}

void print_coefficients(const std::string& label, const std::vector<double>& coefficients, const std::size_t moment_count)
{
    std::cout << label << '\n';
    for (std::size_t basis = 0; basis < basis_size; ++basis)
    {
        std::cout << "  basis " << basis << ':';
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            std::cout << ' ' << coefficients[coefficient_index(basis, moment, moment_count)];
        }
        std::cout << '\n';
    }
}

} // namespace

int main()
{
    try
    {
        // One one-dimensional QMOM population on one quadratic reference
        // element is enough to show the coupling. There is no mesh or
        // physical-space operator in this example.
        constexpr std::size_t node_count = 2;
        const double initial_time = 0.0;
        const double final_time = 0.05;
        const double dt = 0.01;
        const std::size_t time_step_count = 5;
        if (std::abs(initial_time + static_cast<double>(time_step_count) * dt - final_time) > 1.0e-14)
        {
            throw std::runtime_error("The common time interval is not an integer number of time steps.");
        }

        // Keep the phase-space model and timing identical in both routes. The
        // zero source leaves a pure phase-space advection problem.
        const qmomkit::QmomOptions qmom_options;
        const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.1}; };
        const qmomkit::MomentSource zero_source = [](double, const qmomkit::QmomSourceContext& context)
        { return std::vector<double>(context.size(), 0.0); };

        qmomkit::QmomPhaseSpaceEvolutionOptions evolution_options;
        evolution_options.qmom_options = qmom_options;
        evolution_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
        evolution_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

        const ReferenceElement element = make_quadratic_reference_element();

        // QMOM nodes -0.5 and 1.2 with weights 0.4 and 0.6 define this local
        // raw-moment sequence. Positive spatial scaling keeps every local
        // sequence at the spatial quadrature points realizable.
        const std::vector<double> base_moments{1.0, 0.52, 0.964, 0.9868};

        // Conceptually, initial_coefficients[i][k] is raw moment k in spatial
        // basis function i. These are moment fields, not fields of QMOM nodes
        // or weights.
        const std::vector<double> initial_coefficients = make_spatially_varying_coefficients(base_moments);

        // ---------------------------------------------------------------------
        // Route A: instantaneous QMOMKit rates, SSP-RK3 owned by the host
        // ---------------------------------------------------------------------
        auto rate_coefficients = initial_coefficients;
        for (std::size_t step = 0; step < time_step_count; ++step)
        {
            const double time = initial_time + static_cast<double>(step) * dt;

            // QMOMKit does not advance time in this route. Each host RK stage
            // needs a new coefficient-space RHS from the current stage state.

            // Stage 1 starts at the current Galerkin coefficient state.
            const std::vector<double> stage_one_local_moments =
                interpolate_coefficients_to_local_moments(element, rate_coefficients, base_moments.size());
            std::vector<double> stage_one_local_rates(quadrature_count * base_moments.size(), 0.0);
            for (std::size_t q = 0; q < quadrature_count; ++q)
            {
                const qmomkit::QmomAdvectionMomentRateResult rate = qmomkit::evaluate_qmom_advection_moment_rates(
                    local_moment_slice(stage_one_local_moments, q, base_moments.size()), node_count, time, velocity, qmom_options);
                if (!rate.ok())
                {
                    throw std::runtime_error("Route A stage 1 rate evaluation failed at spatial quadrature point " + std::to_string(q) + ": " +
                                             rate.message);
                }
                for (std::size_t moment = 0; moment < base_moments.size(); ++moment)
                {
                    stage_one_local_rates[q * base_moments.size() + moment] = rate.rates[moment];
                }
            }
            const std::vector<double> stage_one_rhs = project_local_values_to_coefficients(element, stage_one_local_rates, base_moments.size());
            std::vector<double> stage_one_coefficients(rate_coefficients.size(), 0.0);
            for (std::size_t value = 0; value < rate_coefficients.size(); ++value)
            {
                stage_one_coefficients[value] = rate_coefficients[value] + dt * stage_one_rhs[value];
            }

            // Stage 2 must revisit the spatial Gauss points: its local moment
            // states differ from those in stage 1, so their rates do as well.
            const std::vector<double> stage_two_local_moments =
                interpolate_coefficients_to_local_moments(element, stage_one_coefficients, base_moments.size());
            std::vector<double> stage_two_local_rates(quadrature_count * base_moments.size(), 0.0);
            for (std::size_t q = 0; q < quadrature_count; ++q)
            {
                const qmomkit::QmomAdvectionMomentRateResult rate = qmomkit::evaluate_qmom_advection_moment_rates(
                    local_moment_slice(stage_two_local_moments, q, base_moments.size()), node_count, time + dt, velocity, qmom_options);
                if (!rate.ok())
                {
                    throw std::runtime_error("Route A stage 2 rate evaluation failed at spatial quadrature point " + std::to_string(q) + ": " +
                                             rate.message);
                }
                for (std::size_t moment = 0; moment < base_moments.size(); ++moment)
                {
                    stage_two_local_rates[q * base_moments.size() + moment] = rate.rates[moment];
                }
            }
            const std::vector<double> stage_two_rhs = project_local_values_to_coefficients(element, stage_two_local_rates, base_moments.size());
            std::vector<double> stage_two_coefficients(rate_coefficients.size(), 0.0);
            for (std::size_t value = 0; value < rate_coefficients.size(); ++value)
            {
                stage_two_coefficients[value] = 0.75 * rate_coefficients[value] + 0.25 * (stage_one_coefficients[value] + dt * stage_two_rhs[value]);
            }

            // Stage 3 makes the same fresh local calls before the host finishes
            // its SSP-RK3 update.
            const std::vector<double> stage_three_local_moments =
                interpolate_coefficients_to_local_moments(element, stage_two_coefficients, base_moments.size());
            std::vector<double> stage_three_local_rates(quadrature_count * base_moments.size(), 0.0);
            for (std::size_t q = 0; q < quadrature_count; ++q)
            {
                const qmomkit::QmomAdvectionMomentRateResult rate = qmomkit::evaluate_qmom_advection_moment_rates(
                    local_moment_slice(stage_three_local_moments, q, base_moments.size()), node_count, time + 0.5 * dt, velocity, qmom_options);
                if (!rate.ok())
                {
                    throw std::runtime_error("Route A stage 3 rate evaluation failed at spatial quadrature point " + std::to_string(q) + ": " +
                                             rate.message);
                }
                for (std::size_t moment = 0; moment < base_moments.size(); ++moment)
                {
                    stage_three_local_rates[q * base_moments.size() + moment] = rate.rates[moment];
                }
            }
            const std::vector<double> stage_three_rhs = project_local_values_to_coefficients(element, stage_three_local_rates, base_moments.size());
            for (std::size_t value = 0; value < rate_coefficients.size(); ++value)
            {
                rate_coefficients[value] =
                    (rate_coefficients[value] / 3.0) + (2.0 / 3.0) * (stage_two_coefficients[value] + dt * stage_three_rhs[value]);
            }
        }

        // ---------------------------------------------------------------------
        // Route B: finite local phase-space evolution, SSP-RK3 owned by QMOMKit
        // ---------------------------------------------------------------------
        auto finite_step_coefficients = initial_coefficients;
        for (std::size_t step = 0; step < time_step_count; ++step)
        {
            const double time = initial_time + static_cast<double>(step) * dt;

            // Interpolation still belongs to the Galerkin side. Each spatial
            // quadrature point now defines a separate local phase-space IVP;
            // QMOMKit sees none of the basis or mass-matrix machinery.
            const std::vector<double> local_moments =
                interpolate_coefficients_to_local_moments(element, finite_step_coefficients, base_moments.size());
            std::vector<double> evolved_local_moments(quadrature_count * base_moments.size(), 0.0);

            for (std::size_t q = 0; q < quadrature_count; ++q)
            {
                try
                {
                    // Build the QMOM state for this spatial quadrature point and
                    // advance that local phase-space problem over dt.

                    // The state includes the corresponding QMOM quadrature for
                    // the local calculation. Afterward, only its evolved raw
                    // moments return to the Galerkin basis; QMOM nodes and
                    // weights never become spatial fields.
                    const qmomkit::QmomState state =
                        qmomkit::QmomCase{qmomkit::MomentSet(local_moment_slice(local_moments, q, base_moments.size())), node_count, qmom_options}
                            .initial_state();
                    const qmomkit::QmomPhaseSpaceEvolutionResult evolved =
                        qmomkit::advance_qmom_phase_space(state, node_count, time, dt, velocity, zero_source, evolution_options);
                    if (!evolved.ok())
                    {
                        throw std::runtime_error(evolved.message);
                    }
                    for (std::size_t moment = 0; moment < base_moments.size(); ++moment)
                    {
                        evolved_local_moments[q * base_moments.size() + moment] = evolved.state.moments()[moment];
                    }
                }
                catch (const std::exception& error)
                {
                    throw std::runtime_error("Route B finite-step evolution failed at spatial quadrature point " + std::to_string(q) + ": " +
                                             error.what());
                }
            }

            // Bring the evolved local moments back to the basis. A real host
            // can place this finite phase-space substep wherever it belongs in
            // its outer splitting sequence.
            finite_step_coefficients = project_local_values_to_coefficients(element, evolved_local_moments, base_moments.size());
        }

        std::cout << std::setprecision(6);
        std::cout << "QMOM spatial-coupling demonstrator\n";
        std::cout << "Common SSP-RK3 interval: t = " << initial_time << " to " << final_time << " in " << time_step_count << " steps of " << dt
                  << "\n";
        print_coefficients("Initial raw-moment coefficients:", initial_coefficients, base_moments.size());
        print_coefficients("Route A: host SSP-RK3 using QMOMKit moment rates:", rate_coefficients, base_moments.size());
        print_coefficients("Route B: QMOMKit SSP-RK3 local phase-space evolution:", finite_step_coefficients, base_moments.size());
        std::cout << "The routes are reported independently; equality is not a validation condition.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << "moment-rate element-operator failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
