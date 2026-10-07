#include "cqmom_case_matrix.hpp"

#include <cmath>
#include <optional>
#include <utility>

namespace qmomkit::case_matrix
{
namespace
{

[[nodiscard]] double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t count = 0; count < exponent; ++count)
    {
        value *= base;
    }
    return value;
}

/*
 * Build the tensor-product node list from per-dimension abscissas and weights.
 * A tensor product of one-dimensional quadratures is itself a quadrature, so
 * the moments derived from it are realizable by construction.
 */
[[nodiscard]] std::vector<QuadratureNode> product_nodes(const std::vector<std::vector<double>>& abscissas,
                                                        const std::vector<std::vector<double>>& weights)
{
    std::vector<QuadratureNode> nodes{QuadratureNode{1.0, {}}};
    for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
    {
        std::vector<QuadratureNode> expanded;
        expanded.reserve(nodes.size() * abscissas[dimension].size());
        for (const QuadratureNode& node : nodes)
        {
            for (std::size_t index = 0; index < abscissas[dimension].size(); ++index)
            {
                QuadratureNode child = node;
                child.weight *= weights[dimension][index];
                child.coordinate.push_back(abscissas[dimension][index]);
                expanded.push_back(std::move(child));
            }
        }
        nodes = std::move(expanded);
    }
    return nodes;
}

// Two well-separated nodes per dimension, offset per dimension so that no two
// coordinates share a value and the conditioning order is observable.
[[nodiscard]] StaticCase two_node_case(std::string name, const std::size_t dimensions, std::vector<std::size_t> conditioning_order)
{
    std::vector<std::vector<double>> abscissas;
    std::vector<std::vector<double>> weights;
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        const double offset = 0.25 * static_cast<double>(dimension);
        abscissas.push_back({-1.0 - offset, 1.0 + offset});
        weights.push_back({0.4 + 0.05 * static_cast<double>(dimension), 0.6 - 0.05 * static_cast<double>(dimension)});
    }

    StaticCase case_;
    case_.name = std::move(name);
    case_.node_counts.assign(dimensions, 2);
    case_.conditioning_order = std::move(conditioning_order);
    case_.nodes = product_nodes(abscissas, weights);
    return case_;
}

[[nodiscard]] std::vector<std::size_t> identity_order(const std::size_t dimensions)
{
    std::vector<std::size_t> order(dimensions, 0);
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        order[dimension] = dimension;
    }
    return order;
}

/*
 * A case whose second coordinate carries only one distinct abscissa. With
 * adaptive branch counts the accepted topology uses fewer nodes there than the
 * configured maximum, which is what makes the adaptive path observable.
 */
[[nodiscard]] StaticCase collapsed_second_dimension_case(std::string name, const std::size_t dimensions)
{
    std::vector<std::vector<double>> abscissas;
    std::vector<std::vector<double>> weights;
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        if (dimension == 1)
        {
            abscissas.push_back({0.5});
            weights.push_back({1.0});
            continue;
        }
        const double offset = 0.25 * static_cast<double>(dimension);
        abscissas.push_back({-1.0 - offset, 1.0 + offset});
        weights.push_back({0.45, 0.55});
    }

    StaticCase case_;
    case_.name = std::move(name);
    case_.node_counts.assign(dimensions, 2);
    case_.conditioning_order = identity_order(dimensions);
    case_.nodes = product_nodes(abscissas, weights);
    case_.adaptive_branch_counts = true;
    return case_;
}

/*
 * The solver sweep uses its own node positions rather than reusing the plain
 * three-dimensional case, so that no solver case duplicates a case from the
 * dimension sweep. The first coordinate is deliberately clustered: a narrow
 * spread there makes the conditional Vandermonde systems less trivial than the
 * widely separated positions used elsewhere.
 */
[[nodiscard]] StaticCase solver_case(const std::string& solver_name, const VandermondeSolver solver)
{
    StaticCase case_;
    case_.name = "static_3d_solver_" + solver_name;
    case_.node_counts = {2, 2, 2};
    case_.conditioning_order = identity_order(3);
    case_.nodes = product_nodes({{-0.15, 0.15}, {-1.0, 1.0}, {-1.25, 1.25}}, {{0.45, 0.55}, {0.4, 0.6}, {0.5, 0.5}});
    case_.vandermonde_solver = solver;
    return case_;
}

} // namespace

std::vector<StaticCase> static_cases()
{
    std::vector<StaticCase> cases;

    // Phase-space dimension coverage.
    cases.push_back(two_node_case("static_2d_fixed_identity", 2, {0, 1}));
    cases.push_back(two_node_case("static_3d_fixed_identity", 3, {0, 1, 2}));
    cases.push_back(two_node_case("static_5d_fixed_identity", 5, identity_order(5)));

    // Conditioning order is distinct from phase-space dimension order, so a
    // reversed order must be exercised separately.
    cases.push_back(two_node_case("static_2d_fixed_reverse_order", 2, {1, 0}));
    cases.push_back(two_node_case("static_3d_fixed_reverse_order", 3, {2, 1, 0}));

    // Unequal node counts across dimensions.
    {
        StaticCase case_;
        case_.name = "static_2d_fixed_unequal_counts";
        case_.node_counts = {3, 2};
        case_.conditioning_order = {0, 1};
        case_.nodes = product_nodes({{-1.5, 0.25, 1.75}, {-0.75, 1.25}}, {{0.3, 0.45, 0.25}, {0.55, 0.45}});
        cases.push_back(std::move(case_));
    }

    // Adaptive branch counts, where the accepted topology is narrower than the
    // configured maximum.
    cases.push_back(collapsed_second_dimension_case("static_3d_adaptive_collapsed_branch", 3));
    cases.push_back(collapsed_second_dimension_case("static_5d_adaptive_collapsed_branch", 5));

    // Bounded coordinate support on one coordinate.
    {
        StaticCase case_;
        case_.name = "static_3d_bounded_support";
        case_.node_counts = {2, 2, 2};
        case_.conditioning_order = {0, 1, 2};
        case_.nodes = product_nodes({{-1.0, 1.0}, {0.2, 0.8}, {-1.25, 1.25}}, {{0.4, 0.6}, {0.5, 0.5}, {0.45, 0.55}});
        case_.coordinate_support = {std::nullopt, AbscissaBounds{0.0, 1.0}, std::nullopt};
        cases.push_back(std::move(case_));
    }

    // Every Vandermonde solver policy, each selecting the linear solver for the
    // conditional systems.
    cases.push_back(solver_case("automatic", VandermondeSolver::automatic));
    cases.push_back(solver_case("dense", VandermondeSolver::dense));
    cases.push_back(solver_case("pivoted_qr", VandermondeSolver::pivoted_qr));
    cases.push_back(solver_case("svd", VandermondeSolver::svd));
    cases.push_back(solver_case("moore_penrose", VandermondeSolver::moore_penrose));
    cases.push_back(solver_case("bjorck_pereyra", VandermondeSolver::bjorck_pereyra));

    return cases;
}

std::vector<EvolutionCase> evolution_cases()
{
    std::vector<EvolutionCase> cases;

    // All three splitting schemes over the same reconstruction, so that a
    // difference between them is attributable to the splitting alone.
    for (const auto& [scheme, suffix] :
         {std::pair{PhaseSpaceEvolutionScheme::strang, "strang"}, std::pair{PhaseSpaceEvolutionScheme::advsrc, "advsrc"},
          std::pair{PhaseSpaceEvolutionScheme::srcadv, "srcadv"}})
    {
        EvolutionCase case_;
        case_.reconstruction = two_node_case(std::string("evolution_3d_") + suffix, 3, {0, 1, 2});
        case_.scheme = scheme;
        case_.time_step = 0.05;
        case_.step_count = 3;
        cases.push_back(std::move(case_));
    }

    {
        EvolutionCase case_;
        case_.reconstruction = collapsed_second_dimension_case("evolution_3d_adaptive_strang", 3);
        case_.scheme = PhaseSpaceEvolutionScheme::strang;
        case_.time_step = 0.05;
        case_.step_count = 3;
        cases.push_back(std::move(case_));
    }

    {
        EvolutionCase case_;
        case_.reconstruction = two_node_case("evolution_5d_strang", 5, identity_order(5));
        case_.scheme = PhaseSpaceEvolutionScheme::strang;
        case_.time_step = 0.05;
        case_.step_count = 3;
        cases.push_back(std::move(case_));
    }

    return cases;
}

PhaseSpaceVelocity matrix_velocity()
{
    // Linear in the first coordinate and constant otherwise: enough to move
    // every node, cheap enough to stay exactly reproducible.
    return make_phase_space_velocity(
        [](double, const std::vector<double>& coordinate)
        {
            std::vector<double> rates(coordinate.size(), 0.0);
            for (std::size_t dimension = 0; dimension < coordinate.size(); ++dimension)
            {
                rates[dimension] = 0.03 - 0.01 * static_cast<double>(dimension) + 0.02 * coordinate[0];
            }
            return rates;
        });
}

CqmomMomentSource matrix_source()
{
    // A pure decay of every moment except the zeroth: it exercises the source
    // stage without changing total number density.
    return [](double, const CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t moment = 0; moment < context.size(); ++moment)
        {
            std::size_t total_order = 0;
            const MultiIndexView index = context.index(moment);
            for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
            {
                total_order += index[dimension];
            }
            if (total_order == 0)
            {
                continue;
            }
            rates[moment] = -0.1 * static_cast<double>(total_order) * context.value(moment);
        }
        return CqmomMomentRate{std::move(rates)};
    };
}

namespace
{
PhaseSpaceHierarchy hierarchy_of(const StaticCase& case_)
{
    return PhaseSpaceHierarchy(case_.node_counts, case_.conditioning_order);
}

RawMomentTensor moment_tensor_of(const StaticCase& case_)
{
    const std::vector<std::size_t> extents = hierarchy_of(case_).moment_extents();
    std::vector<double> values(product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        const MultiIndex order = unflatten_multi_index(flat, extents);
        double moment = 0.0;
        for (const QuadratureNode& node : case_.nodes)
        {
            double contribution = node.weight;
            for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
            {
                contribution *= integer_power(node.coordinate[dimension], order[dimension]);
            }
            moment += contribution;
        }
        values[flat] = moment;
    }
    return RawMomentTensor{extents, std::move(values)};
}

CqmomMomentData moments_of(const StaticCase& case_)
{
    const CqmomMomentLayout layout = create_cqmom_moment_layout(hierarchy_of(case_));
    return extract_cqmom_moments(moment_tensor_of(case_), layout);
}
} // namespace

CqmomCase cqmom_case_of(const StaticCase& case_)
{
    // CQMOM requires the QMOM adaptive option to agree with the branch-count flag.
    QmomOptions qmom_options;
    qmom_options.adaptive = case_.adaptive_branch_counts;
    return CqmomCase{moments_of(case_),       {}, qmom_options, case_.adaptive_branch_counts, CqmomDiagnosticLevel::basic, case_.vandermonde_solver,
                     case_.coordinate_support};
}

} // namespace qmomkit::case_matrix
