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

// This is deliberately only a small reference element. The host application,
// not QMOMKit, owns its spatial basis, quadrature rule, and mass matrix.
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

// The host stores moment fields as coefficients[basis][moment]. QMOMKit instead
// consumes one independent raw-moment set, so the host evaluates every moment
// field at each spatial quadrature point before making a local QMOMKit call.
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

// QMOMKit returns values at spatial quadrature points, while the host evolves
// basis coefficients. This is the host's L2 projection: first form the
// weighted basis inner products, then apply M^{-1}. The local values may be
// either instantaneous rates (Route A) or evolved raw moments (Route B).
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
    // The nodal scale field is 1 + 0.15 x. It stays positive on the element,
    // so it scales all weights of the fixed realizable moment sequence by a
    // positive factor without changing its two-node support.
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
        // This manufactured problem has one one-dimensional QMOM population on
        // one quadratic reference element. There is deliberately no mesh or
        // physical-space operator here: a real CG or DG application owns
        // those operations independently of the local QMOMKit calls below.
        constexpr std::size_t node_count = 2;
        const double initial_time = 0.0;
        const double final_time = 0.05;
        const double dt = 0.01;
        const std::size_t time_step_count = 5;
        if (std::abs(initial_time + static_cast<double>(time_step_count) * dt - final_time) > 1.0e-14)
        {
            throw std::runtime_error("The common time interval is not an integer number of time steps.");
        }

        // Both routes use this exact velocity model, this node count, this
        // reconstruction policy, and the same SSP-RK3 time step. The zero
        // source makes this a pure phase-space advection example.
        const qmomkit::QmomOptions qmom_options;
        const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.1}; };
        const qmomkit::MomentSource zero_source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); };

        qmomkit::QmomPhaseSpaceEvolutionOptions evolution_options;
        evolution_options.qmom_options = qmom_options;
        evolution_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
        evolution_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

        const ReferenceElement element = make_quadratic_reference_element();

        // A two-node quadrature with nodes -0.5 and 1.2 and weights 0.4 and
        // 0.6 defines this raw-moment sequence. The spatially varying positive
        // scaling below keeps the local sequences realizable at all three
        // spatial quadrature points.
        const std::vector<double> base_moments{1.0, 0.52, 0.964, 0.9868};

        // Conceptually, initial_coefficients[i][k] is the coefficient of raw
        // moment order k in spatial basis function i. The host has no spatial
        // field of QMOM nodes or weights; those exist only inside local calls.
        const std::vector<double> initial_coefficients = make_spatially_varying_coefficients(base_moments);

        // ---------------------------------------------------------------------
        // Route A: instantaneous QMOMKit rates, SSP-RK3 owned by the host
        // ---------------------------------------------------------------------
        auto rate_coefficients = initial_coefficients;
        for (std::size_t step = 0; step < time_step_count; ++step)
        {
            const double time = initial_time + static_cast<double>(step) * dt;

            // QMOMKit returns only an instantaneous local phase-space rate in
            // this route. The host must form a fresh coefficient-space RHS at
            // every SSP-RK3 stage; a real spatial solver would use its own
            // native time integrator at this same coupling boundary.

            // Stage 1: coefficients -> local moments -> local QMOMKit rates -> projected RHS.
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

            // Stage 2 repeats interpolation, local rate evaluation, and
            // projection for the first SSP-RK3 stage state. Rates from stage 1
            // are not reused because the local raw-moment states changed.
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

            // Stage 3 again makes fresh local QMOMKit calls before completing
            // the host-owned SSP-RK3 update of the coefficient fields.
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

            // The host still interpolates its moment fields. Each resulting
            // point is an independent local phase-space IVP; QMOMKit never
            // receives the basis functions, mass matrix, or spatial rule.
            const std::vector<double> local_moments =
                interpolate_coefficients_to_local_moments(element, finite_step_coefficients, base_moments.size());
            std::vector<double> evolved_local_moments(quadrature_count * base_moments.size(), 0.0);

            for (std::size_t q = 0; q < quadrature_count; ++q)
            {
                try
                {
                    // Construct a local QMOM state from the raw moments at this spatial
                    // quadrature point and advance it over the phase-space substep.

                    // The QMOM state contains the corresponding quadrature, but that quadrature
                    // is only a local phase-space representation. After the substep, the host
                    // takes the evolved raw moments and projects those back to its spatial basis;
                    // QMOM nodes and weights are not represented as spatial fields.
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

            // QMOMKit has advanced each local phase-space IVP over dt. The
            // host projects the evolved raw moments to its basis, then decides
            // where this substep belongs relative to physical-space evolution
            // in an outer operator-splitting scheme.
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
