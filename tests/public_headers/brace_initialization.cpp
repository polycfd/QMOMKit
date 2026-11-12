#include <qmomkit/qmomkit.hpp>
#include <vector>

// Callers brace-initialize the input aggregates and QmomResult, omitting
// trailing members. This translation unit is compiled with
// -Wmissing-field-initializers, so a member without a default initializer
// fails the -Werror CI build.
namespace public_header_brace_initialization
{

[[maybe_unused]] void initialize_with_trailing_members_omitted()
{
    [[maybe_unused]] const qmomkit::QmomCase qmom_case{qmomkit::MomentSet{1.0, 0.0}, 1};
    // The family follows the four older options, which keep their positions.
    [[maybe_unused]] const qmomkit::QmomOptions fixed{false};
    [[maybe_unused]] const qmomkit::QmomOptions hyperbolic{true, 1.0e-10, 1.0e-10, 1.0e-30, qmomkit::QmomFamily::hyperbolic};
    // The settings of the generalized family follow the family.
    [[maybe_unused]] const qmomkit::QmomOptions generalized{true, 1.0e-10, 1.0e-10, 1.0e-30, qmomkit::QmomFamily::generalized, {9}};
    // The base node count follows the four older members of the result.
    [[maybe_unused]] const qmomkit::QmomResult result{{0.5}, {1.0}, qmomkit::QmomStatus::success, ""};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({1, 1}));
    [[maybe_unused]] const qmomkit::CqmomCase cqmom_case{qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size(), 1.0)}};
    [[maybe_unused]] const qmomkit::RealizabilityPolicy policy{qmomkit::RealizabilityCheckMode::FinalOnly};
    const qmomkit::SspRk3Options ssp_rk3_options{true};
    [[maybe_unused]] const qmomkit::MomentSourceUpdateOptions source_options{ssp_rk3_options};
    const qmomkit::NodeAdvectionOptions node_options{qmomkit::NodeAdvectionScheme::rk45};
    [[maybe_unused]] const qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options{node_options};
    [[maybe_unused]] const qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options{node_options};
    [[maybe_unused]] const qmomkit::OutputOptions output_options{"output", true};
    [[maybe_unused]] const qmomkit::RawMomentTensor tensor{{2}};
    [[maybe_unused]] const qmomkit::config::IndependentProductInitialization initialization{qmomkit::IndependentProductDistribution{}};
}

} // namespace public_header_brace_initialization
