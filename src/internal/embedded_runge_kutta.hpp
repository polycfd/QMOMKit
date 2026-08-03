#pragma once

#include "qmomkit/runge_kutta.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace qmomkit::internal
{

constexpr std::size_t maximum_embedded_rk_stages = 7;

struct EmbeddedRkTableau
{
    const char* name = "";
    std::size_t stage_count = 0;
    std::size_t higher_order = 0;
    std::size_t embedded_order = 0;
    std::array<double, maximum_embedded_rk_stages> c{};
    std::array<std::array<double, maximum_embedded_rk_stages>, maximum_embedded_rk_stages> a{};
    std::array<double, maximum_embedded_rk_stages> b{};
    std::array<double, maximum_embedded_rk_stages> b_embedded{};
    bool fsal = false;
};

struct EmbeddedRkStep
{
    std::vector<double> higher_state;
    std::vector<double> embedded_state;
    std::vector<double> error;
    std::vector<double> first_derivative;
    std::vector<double> final_derivative;
};

[[nodiscard]] const EmbeddedRkTableau& embedded_rk_tableau(EmbeddedRkMethod method);

[[nodiscard]] EmbeddedRkStep take_embedded_rk_step(const EmbeddedRkTableau& tableau, const OdeRhs& rhs, double t, const std::vector<double>& state,
                                                   double step, const std::optional<std::vector<double>>& first_derivative,
                                                   std::size_t& function_evaluations);

[[nodiscard]] double embedded_rk_controller_exponent(const EmbeddedRkTableau& tableau);

} // namespace qmomkit::internal
