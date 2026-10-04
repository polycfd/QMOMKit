#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

template <class Function> void assert_throws_invalid_argument(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
}

template <class Function> void assert_throws_invalid_argument_message(Function&& function, const std::string& expected_message)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument& error)
    {
        threw = true;
        assert(error.what() == expected_message);
    }
    assert(threw);
}

template <class Function> void assert_throws_domain_error_message(Function&& function, const std::string& expected_message)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::domain_error& error)
    {
        threw = true;
        assert(error.what() == expected_message);
    }
    assert(threw);
}

template <class Function> void assert_throws_overflow_error_message(Function&& function, const std::string& expected_message)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::overflow_error& error)
    {
        threw = true;
        assert(error.what() == expected_message);
    }
    assert(threw);
}

std::vector<double> manufactured_qmom_moments()
{
    return {1.0, 0.75, 3.625, 6.5625, 27.53125, 66.515625};
}

std::vector<double> manufactured_terminal_moments(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                                  const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                double power = 1.0;
                for (std::size_t order = 0; order < moment_index[dimension]; ++order)
                {
                    power *= abscissas[dimension][node];
                }
                contribution *= power;
            }
            moments[flat_moment] += contribution;
        }
    }
    return moments;
}

void compare_cqmom_moment_data_reconstruction_to_moment_tensor(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<double>& weights,
                                                               const std::vector<std::vector<double>>& abscissas, const double tolerance)
{
    const std::vector<std::size_t> extents = hierarchy.moment_extents();
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());
    assert(cqmom.quadrature.has_rectangular_storage());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentTensor moment_tensor_reconstructed = qmomkit::reconstruct_raw_moments(cqmom.quadrature, extents);
    const qmomkit::CqmomMomentData required_moment_data = qmomkit::extract_cqmom_moments(moment_tensor_reconstructed, layout);
    const qmomkit::CqmomMomentData moment_data_reconstructed = qmomkit::reconstruct_cqmom_moments(cqmom.quadrature, layout);

    assert(moment_data_reconstructed.layout.hierarchy.node_counts() == layout.hierarchy.node_counts());
    assert(moment_data_reconstructed.layout.hierarchy.conditioning_order() == layout.hierarchy.conditioning_order());
    assert(moment_data_reconstructed.layout.required_indices == layout.required_indices);
    assert_vector_near(moment_data_reconstructed.values, required_moment_data.values, tolerance);
}

void test_qmom_reconstructs_manufactured_moments()
{
    const std::vector<double> moments = manufactured_qmom_moments();
    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 3);
    assert(qmom.ok());

    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments({qmom.abscissas, qmom.weights}, moments.size());
    assert_vector_near(reconstructed, moments, 1.0e-11);
}

void test_qmom_evaluates_generalized_moment()
{
    const qmomkit::Quadrature1D quadrature{{1.0, 8.0}, {2.0, 3.0}};

    const double moment = qmomkit::evaluate_generalized_moment(quadrature, 1.0 / 3.0);

    assert_near(moment, 8.0, 1.0e-12);
}

void test_generalized_moment_zero_exponent_returns_total_weight()
{
    const qmomkit::Quadrature1D quadrature{{1.0, 4.0, 9.0}, {2.0, 3.0, 5.0}};

    const double moment = qmomkit::evaluate_generalized_moment(quadrature, 0.0);

    assert_near(moment, 10.0, 1.0e-15);
}

void test_generalized_moment_integer_exponents_match_raw_reconstruction()
{
    const qmomkit::Quadrature1D quadrature{{1.0, 8.0}, {2.0, 3.0}};
    const std::vector<double> raw_moments = qmomkit::reconstruct_raw_moments(quadrature, 3);

    assert_near(qmomkit::evaluate_generalized_moment(quadrature, 1.0), raw_moments[1], 1.0e-15);
    assert_near(qmomkit::evaluate_generalized_moment(quadrature, 2.0), raw_moments[2], 1.0e-15);
}

void test_generalized_moment_accepts_zero_abscissa_with_positive_fractional_exponent()
{
    const qmomkit::Quadrature1D quadrature{{0.0, 4.0}, {2.0, 3.0}};

    const double moment = qmomkit::evaluate_generalized_moment(quadrature, 0.5);

    assert_near(moment, 6.0, 1.0e-15);
}

void test_cqmom_reconstructs_fixed_non_identity_phase_order_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());
    assert(cqmom.quadrature.has_rectangular_storage());

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(cqmom.quadrature, extents);
    assert_vector_near(reconstructed.values, moments, 1.0e-10);
}

void test_cqmom_evaluates_generalized_moment()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{2.0, 3.0, 5.0, 7.0};
    const std::vector<std::vector<double>> abscissas{{1.0, 1.0, 4.0, 4.0}, {9.0, 16.0, 9.0, 16.0}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    const std::vector<double> exponents{0.5, 0.25};
    double expected = 0.0;
    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        expected += weights[node] * std::pow(abscissas[0][node], exponents[0]) * std::pow(abscissas[1][node], exponents[1]);
    }

    const double moment = qmomkit::evaluate_generalized_moment(cqmom.quadrature, exponents);
    assert_near(moment, expected, 1.0e-9);
}

void test_cqmom_generalized_moment_rejects_exponent_count_mismatch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{2.0, 3.0, 5.0, 7.0};
    const std::vector<std::vector<double>> abscissas{{1.0, 1.0, 4.0, 4.0}, {9.0, 16.0, 9.0, 16.0}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    assert_throws_invalid_argument_message([&] { (void)qmomkit::evaluate_generalized_moment(cqmom.quadrature, {0.5}); },
                                           "Generalized moment exponent count must match the CQMOM phase-space dimension count.");
}

void test_cqmom_generalized_moment_rejects_negative_abscissa_for_fractional_exponent()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    assert_throws_domain_error_message([&] { (void)qmomkit::evaluate_generalized_moment(cqmom.quadrature, {0.5, 0.0}); },
                                       "Generalized moment evaluation requires non-negative abscissas for non-integer exponents.");

    assert_near(qmomkit::evaluate_generalized_moment(cqmom.quadrature, {1.0, 0.0}), moments[4], 1.0e-10);
}

void test_cqmom_reconstructs_selected_requested_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    const qmomkit::RawMomentRequest request{{{1, 0}, {0, 2}, {1, 1}}};
    const qmomkit::SelectedRawMomentResult reconstructed = qmomkit::reconstruct_raw_moments(cqmom.quadrature, request);

    assert(reconstructed.ok());
    assert(reconstructed.indices == request.indices);
    assert_vector_near(reconstructed.values, {moments[4], moments[2], moments[5]}, 1.0e-10);
}

void test_cqmom_selected_reconstruction_preserves_duplicate_requests()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    const qmomkit::RawMomentRequest request{{{1, 0}, {1, 0}, {0, 0}}};
    const qmomkit::SelectedRawMomentResult reconstructed = qmomkit::reconstruct_raw_moments(cqmom.quadrature, request);

    assert(reconstructed.ok());
    assert(reconstructed.indices == request.indices);
    assert_vector_near(reconstructed.values, {moments[4], moments[4], moments[0]}, 1.0e-10);
}

void test_cqmom_moment_tensor_reconstruction_rejects_rank_mismatch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    assert_throws_invalid_argument_message([&] { (void)qmomkit::reconstruct_raw_moments(cqmom.quadrature, std::vector<std::size_t>{4}); },
                                           "Moment reconstruction extents must match the CQMOM phase-space dimension count.");
}

void test_cqmom_selected_reconstruction_rejects_rank_mismatch()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    assert_throws_invalid_argument([&] { (void)qmomkit::reconstruct_raw_moments(cqmom.quadrature, qmomkit::RawMomentRequest{{{1}}}); });
}

void test_cqmom_moment_data_reconstructs_2d_required_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};

    compare_cqmom_moment_data_reconstruction_to_moment_tensor(hierarchy, weights, abscissas, 1.0e-10);
}

void test_cqmom_moment_data_reconstructs_3d_required_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<double> weights{0.06, 0.09, 0.14, 0.21, 0.10, 0.15, 0.10, 0.15};
    const std::vector<std::vector<double>> abscissas{
        {-1.0, -1.0, -1.0, -1.0, 1.5, 1.5, 1.5, 1.5}, {-0.75, -0.75, 1.25, 1.25, -0.75, -0.75, 1.25, 1.25}, {0.2, 1.4, 0.2, 1.4, 0.2, 1.4, 0.2, 1.4}};

    compare_cqmom_moment_data_reconstruction_to_moment_tensor(hierarchy, weights, abscissas, 1.0e-9);
}

void test_cqmom_moment_data_preserves_non_identity_conditioning_order()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});
    const std::vector<double> weights{0.06, 0.09, 0.14, 0.21, 0.10, 0.15, 0.10, 0.15};
    const std::vector<std::vector<double>> abscissas{
        {-1.0, -1.0, -1.0, -1.0, 1.5, 1.5, 1.5, 1.5}, {-0.75, -0.75, 1.25, 1.25, -0.75, -0.75, 1.25, 1.25}, {0.2, 1.4, 0.2, 1.4, 0.2, 1.4, 0.2, 1.4}};

    compare_cqmom_moment_data_reconstruction_to_moment_tensor(hierarchy, weights, abscissas, 1.0e-9);
}

void test_cqmom_moment_data_reconstructs_unequal_node_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 1, 2}, {2, 0, 1});
    const std::vector<double> weights{0.12, 0.28, 0.18, 0.42};
    const std::vector<std::vector<double>> abscissas{{-1.2, -1.2, 1.8, 1.8}, {0.75, 0.75, 0.75, 0.75}, {-0.4, 1.6, -0.4, 1.6}};

    compare_cqmom_moment_data_reconstruction_to_moment_tensor(hierarchy, weights, abscissas, 1.0e-9);
}

void test_cqmom_reconstructs_adaptive_terminal_moments()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom =
        qmomkit::compute_cqmom(qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy)),
                                                  {},
                                                  {true, 1.0e-10, 1.0e-10, 1.0e-12},
                                                  true});
    assert(cqmom.ok());
    assert(!cqmom.quadrature.has_rectangular_storage());

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(cqmom.quadrature, extents);
    assert_vector_near(reconstructed.values, moments, 1.0e-10);
}

void test_cqmom_moment_data_reconstructs_adaptive_quadrature()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom =
        qmomkit::compute_cqmom(qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy)),
                                                  {},
                                                  {true, 1.0e-10, 1.0e-10, 1.0e-12},
                                                  true});
    assert(cqmom.ok());
    assert(!cqmom.quadrature.has_rectangular_storage());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(cqmom.quadrature, layout);
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments({extents, moments}, layout);
    assert_vector_near(reconstructed.values, expected.values, 1.0e-10);
}

void test_cqmom_moment_data_rejects_mismatched_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    const qmomkit::CqmomMomentLayout wrong_layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 1}));
    assert_throws_invalid_argument_message([&] { (void)qmomkit::reconstruct_cqmom_moments(cqmom.quadrature, wrong_layout); },
                                           "CQMOM moment layout hierarchy must match the quadrature hierarchy.");
}

void test_cqmom_moment_data_rejects_noncanonical_fixed_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const std::vector<double> moments = manufactured_terminal_moments(extents, weights, abscissas);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, moments}, qmomkit::create_cqmom_moment_layout(hierarchy))});
    assert(cqmom.ok());

    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    layout.required_indices.pop_back();
    assert_throws_invalid_argument_message([&] { (void)qmomkit::reconstruct_cqmom_moments(cqmom.quadrature, layout); },
                                           "CQMOM moment layout must be the fixed required layout for its hierarchy.");
}

void test_reconstruction_rejects_invalid_qmom_storage()
{
    bool threw = false;
    try
    {
        (void)qmomkit::reconstruct_raw_moments(qmomkit::Quadrature1D{{1.0}, {0.5, 0.5}}, 2);
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
}

void test_reconstruction_reports_nonfinite_output()
{
    bool threw = false;
    try
    {
        (void)qmomkit::reconstruct_raw_moments(qmomkit::Quadrature1D{{1.0e200}, {1.0}}, 4);
    }
    catch (const std::overflow_error&)
    {
        threw = true;
    }
    assert(threw);
}

void test_generalized_moment_rejects_negative_abscissa_for_non_integer_exponent()
{
    const qmomkit::Quadrature1D quadrature{{-1.0}, {1.0}};

    assert_throws_domain_error_message([&] { (void)qmomkit::evaluate_generalized_moment(quadrature, 0.5); },
                                       "Generalized moment evaluation requires non-negative abscissas for non-integer exponents.");
}

void test_generalized_moment_reports_zero_abscissa_negative_exponent_overflow()
{
    const qmomkit::Quadrature1D quadrature{{0.0}, {1.0}};

    assert_throws_overflow_error_message([&] { (void)qmomkit::evaluate_generalized_moment(quadrature, -1.0); },
                                         "Generalized moment evaluation produced a non-finite value.");
}

void test_generalized_moment_rejects_nonfinite_exponent()
{
    const qmomkit::Quadrature1D quadrature{{1.0}, {1.0}};

    assert_throws_invalid_argument_message([&] { (void)qmomkit::evaluate_generalized_moment(quadrature, std::numeric_limits<double>::infinity()); },
                                           "Generalized moment exponent must be finite.");
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"qmom_reconstructs_manufactured_moments", test_qmom_reconstructs_manufactured_moments},
            {"qmom_evaluates_generalized_moment", test_qmom_evaluates_generalized_moment},
            {"generalized_moment_zero_exponent_returns_total_weight", test_generalized_moment_zero_exponent_returns_total_weight},
            {"generalized_moment_integer_exponents_match_raw_reconstruction", test_generalized_moment_integer_exponents_match_raw_reconstruction},
            {"generalized_moment_accepts_zero_abscissa_with_positive_fractional_exponent",
             test_generalized_moment_accepts_zero_abscissa_with_positive_fractional_exponent},
            {"cqmom_reconstructs_fixed_non_identity_phase_order_moments", test_cqmom_reconstructs_fixed_non_identity_phase_order_moments},
            {"cqmom_evaluates_generalized_moment", test_cqmom_evaluates_generalized_moment},
            {"cqmom_generalized_moment_rejects_exponent_count_mismatch", test_cqmom_generalized_moment_rejects_exponent_count_mismatch},
            {"cqmom_generalized_moment_rejects_negative_abscissa_for_fractional_exponent",
             test_cqmom_generalized_moment_rejects_negative_abscissa_for_fractional_exponent},
            {"cqmom_reconstructs_selected_requested_moments", test_cqmom_reconstructs_selected_requested_moments},
            {"cqmom_selected_reconstruction_preserves_duplicate_requests", test_cqmom_selected_reconstruction_preserves_duplicate_requests},
            {"cqmom_moment_tensor_reconstruction_rejects_rank_mismatch", test_cqmom_moment_tensor_reconstruction_rejects_rank_mismatch},
            {"cqmom_selected_reconstruction_rejects_rank_mismatch", test_cqmom_selected_reconstruction_rejects_rank_mismatch},
            {"cqmom_moment_data_reconstructs_2d_required_moments", test_cqmom_moment_data_reconstructs_2d_required_moments},
            {"cqmom_moment_data_reconstructs_3d_required_moments", test_cqmom_moment_data_reconstructs_3d_required_moments},
            {"cqmom_moment_data_preserves_non_identity_conditioning_order", test_cqmom_moment_data_preserves_non_identity_conditioning_order},
            {"cqmom_moment_data_reconstructs_unequal_node_counts", test_cqmom_moment_data_reconstructs_unequal_node_counts},
            {"cqmom_reconstructs_adaptive_terminal_moments", test_cqmom_reconstructs_adaptive_terminal_moments},
            {"cqmom_moment_data_reconstructs_adaptive_quadrature", test_cqmom_moment_data_reconstructs_adaptive_quadrature},
            {"cqmom_moment_data_rejects_mismatched_layout", test_cqmom_moment_data_rejects_mismatched_layout},
            {"cqmom_moment_data_rejects_noncanonical_fixed_layout", test_cqmom_moment_data_rejects_noncanonical_fixed_layout},
            {"reconstruction_rejects_invalid_qmom_storage", test_reconstruction_rejects_invalid_qmom_storage},
            {"reconstruction_reports_nonfinite_output", test_reconstruction_reports_nonfinite_output},
            {"generalized_moment_rejects_negative_abscissa_for_non_integer_exponent",
             test_generalized_moment_rejects_negative_abscissa_for_non_integer_exponent},
            {"generalized_moment_reports_zero_abscissa_negative_exponent_overflow",
             test_generalized_moment_reports_zero_abscissa_negative_exponent_overflow},
            {"generalized_moment_rejects_nonfinite_exponent", test_generalized_moment_rejects_nonfinite_exponent},
        });
}
