#include "embedded_runge_kutta.hpp"

#include "numeric_utils.hpp"

#include <cmath>
#include <stdexcept>

namespace qmomkit::internal
{
namespace
{

constexpr EmbeddedRkTableau bogacki_shampine_rk23{"Bogacki-Shampine RK3(2)",
                                                  4,
                                                  3,
                                                  2,
                                                  {0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0, 0.0, 0.0, 0.0},
                                                  {{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                    {1.0 / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                    {0.0, 3.0 / 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                    {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
                                                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}},
                                                  {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
                                                  {7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0, 0.0, 0.0, 0.0},
                                                  true};

constexpr EmbeddedRkTableau dormand_prince_rk45{
    "Dormand-Prince RK5(4)7M",
    7,
    5,
    4,
    {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0},
    {{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
      {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
      {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}}},
    {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0},
    {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0},
    true};

void require_same_size(const std::vector<double>& lhs, const std::vector<double>& rhs, const char* message)
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument(message);
    }
}

[[nodiscard]] std::vector<double> rhs_checked(const OdeRhs& rhs, const double time, const std::vector<double>& state)
{
    std::vector<double> derivative = rhs(time, state);
    require_same_size(derivative, state, "RHS returned a vector with the wrong size.");
    if (!all_finite(derivative))
    {
        throw std::invalid_argument("RHS returned a non-finite value.");
    }
    return derivative;
}

} // namespace

const EmbeddedRkTableau& embedded_rk_tableau(const EmbeddedRkMethod method)
{
    switch (method)
    {
    case EmbeddedRkMethod::rk23:
        return bogacki_shampine_rk23;
    case EmbeddedRkMethod::rk45:
        return dormand_prince_rk45;
    }
    throw std::invalid_argument("Unsupported embedded Runge-Kutta method.");
}

EmbeddedRkStep take_embedded_rk_step(const EmbeddedRkTableau& tableau, const OdeRhs& rhs, const double t, const std::vector<double>& state,
                                     const double step, const std::optional<std::vector<double>>& first_derivative, std::size_t& function_evaluations)
{
    std::array<std::vector<double>, maximum_embedded_rk_stages> derivatives;
    if (first_derivative.has_value())
    {
        require_same_size(*first_derivative, state, "Cached RHS vector has the wrong size.");
        derivatives[0] = *first_derivative;
    }
    else
    {
        derivatives[0] = rhs_checked(rhs, t, state);
        ++function_evaluations;
    }

    // Reused across stages: every stage resets it to the incoming state before
    // accumulating, so only the first iteration needs to allocate.
    std::vector<double> stage_state;
    for (std::size_t stage = 1; stage < tableau.stage_count; ++stage)
    {
        stage_state.assign(state.begin(), state.end());
        for (std::size_t previous = 0; previous < stage; ++previous)
        {
            const double coefficient = step * tableau.a[stage][previous];
            for (std::size_t component = 0; component < stage_state.size(); ++component)
            {
                stage_state[component] += coefficient * derivatives[previous][component];
            }
        }
        derivatives[stage] = rhs_checked(rhs, t + tableau.c[stage] * step, stage_state);
        ++function_evaluations;
    }

    std::vector<double> higher_state = state;
    std::vector<double> embedded_state = state;
    for (std::size_t stage = 0; stage < tableau.stage_count; ++stage)
    {
        for (std::size_t component = 0; component < state.size(); ++component)
        {
            higher_state[component] += step * tableau.b[stage] * derivatives[stage][component];
            embedded_state[component] += step * tableau.b_embedded[stage] * derivatives[stage][component];
        }
    }

    std::vector<double> error(state.size(), 0.0);
    for (std::size_t component = 0; component < state.size(); ++component)
    {
        error[component] = higher_state[component] - embedded_state[component];
    }

    return EmbeddedRkStep{std::move(higher_state), std::move(embedded_state), std::move(error), std::move(derivatives[0]),
                          std::move(derivatives[tableau.stage_count - 1])};
}

double embedded_rk_controller_exponent(const EmbeddedRkTableau& tableau)
{
    return -1.0 / static_cast<double>(tableau.embedded_order + 1);
}

} // namespace qmomkit::internal
