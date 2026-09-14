#include "embedded_runge_kutta.hpp"

#include "ode_rhs_validation.hpp"

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

void take_embedded_rk_step(const EmbeddedRkTableau& tableau, const OdeRhs& rhs, const double t, const std::vector<double>& state, const double step,
                           std::size_t& function_evaluations, EmbeddedRkWorkspace& workspace)
{
    if (workspace.has_first_derivative)
    {
        if (workspace.derivatives[0].size() != state.size())
        {
            throw std::invalid_argument("Cached RHS vector has the wrong size.");
        }
    }
    else
    {
        rhs_checked(rhs, t, state, workspace.derivatives[0]);
        ++function_evaluations;
        workspace.has_first_derivative = true;
    }

    // Every buffer below is assigned its full length before it is read, so a
    // workspace that has already served an attempt reuses the storage from that
    // attempt and only the first attempt of an integration allocates.
    for (std::size_t stage = 1; stage < tableau.stage_count; ++stage)
    {
        workspace.stage_state.assign(state.begin(), state.end());
        for (std::size_t previous = 0; previous < stage; ++previous)
        {
            const double coefficient = step * tableau.a[stage][previous];
            for (std::size_t component = 0; component < workspace.stage_state.size(); ++component)
            {
                workspace.stage_state[component] += coefficient * workspace.derivatives[previous][component];
            }
        }
        rhs_checked(rhs, t + tableau.c[stage] * step, workspace.stage_state, workspace.derivatives[stage]);
        ++function_evaluations;
    }

    workspace.higher_state.assign(state.begin(), state.end());
    workspace.embedded_state.assign(state.begin(), state.end());
    for (std::size_t stage = 0; stage < tableau.stage_count; ++stage)
    {
        for (std::size_t component = 0; component < state.size(); ++component)
        {
            workspace.higher_state[component] += step * tableau.b[stage] * workspace.derivatives[stage][component];
            workspace.embedded_state[component] += step * tableau.b_embedded[stage] * workspace.derivatives[stage][component];
        }
    }

    workspace.error.resize(state.size());
    for (std::size_t component = 0; component < state.size(); ++component)
    {
        workspace.error[component] = workspace.higher_state[component] - workspace.embedded_state[component];
    }
}

double embedded_rk_controller_exponent(const EmbeddedRkTableau& tableau)
{
    return -1.0 / static_cast<double>(tableau.embedded_order + 1);
}

} // namespace qmomkit::internal
