#include <qmomkit/qmomkit.hpp>
#include <vector>

// Callers brace-initialize the input aggregates and omit trailing members. This
// translation unit is compiled with -Wmissing-field-initializers, so a member
// without a default initializer fails the -Werror CI build.
namespace public_header_brace_initialization
{

[[maybe_unused]] void initialize_with_trailing_members_omitted()
{
    [[maybe_unused]] const qmomkit::QmomCase qmom_case{qmomkit::MomentSet{1.0, 0.0}, 1};
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
