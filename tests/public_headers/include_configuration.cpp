// Only the header under test is included. The pins below use std::filesystem,
// std::string and std::string_view through it, so this file also checks that
// the header provides the standard types its own declarations need.
#include <qmomkit/configuration.hpp>

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature. A
// changed parameter or return type fails to compile here. The pins have
// external linkage, so each one is emitted and its function must be defined
// when this test is linked. In an unnamed namespace the compiler may drop them,
// as Clang does even without optimization, and the definitions would go
// unchecked.

using QmomMomentInitializer = qmomkit::config::ConfigurationResult<qmomkit::MomentSet> (*)(const qmomkit::config::InitialConditionSpecification&,
                                                                                           std::string);

using CqmomDenseMomentInitializer =
    qmomkit::config::ConfigurationResult<qmomkit::RawMomentTensor> (*)(const qmomkit::config::InitialConditionSpecification&, std::string);

using CqmomMomentDataInitializer = qmomkit::config::ConfigurationResult<qmomkit::CqmomMomentData> (*)(
    const qmomkit::config::InitialConditionSpecification&, const qmomkit::PhaseSpaceHierarchy&, std::string);

using CqmomLeanCaseFactory = qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase> (*)(const qmomkit::config::InitialConditionSpecification&,
                                                                                              const qmomkit::PhaseSpaceHierarchy&,
                                                                                              const qmomkit::QmomOptions&, bool, std::string,
                                                                                              qmomkit::CqmomCoordinateSupport);

using CqmomLeanCaseParser = qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase> (*)(std::string_view, std::string);

using CqmomPreparedCaseParser = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(std::string_view, std::string);

using CqmomLeanCaseLoader = qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase> (*)(const std::filesystem::path&);

using CqmomPreparedCaseLoader = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const std::filesystem::path&);

using CqmomInitialMomentReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const qmomkit::CqmomPreparedCase&,
                                                                                                           qmomkit::RawMomentTensor);

using CqmomLeanInitialMomentReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const qmomkit::CqmomPreparedCase&,
                                                                                                               qmomkit::CqmomMomentData);

QmomMomentInitializer qmom_moment_initializer = qmomkit::config::initialize_moments;

CqmomDenseMomentInitializer cqmom_dense_moment_initializer = qmomkit::config::initialize_cqmom_dense_moments;

CqmomMomentDataInitializer cqmom_moment_data_initializer = qmomkit::config::initialize_cqmom_moment_data;

CqmomLeanCaseFactory cqmom_lean_case_factory = qmomkit::config::create_cqmom_lean_case;

CqmomLeanCaseParser cqmom_lean_case_parser = qmomkit::config::parse_cqmom_lean;

CqmomPreparedCaseParser cqmom_prepared_case_parser = qmomkit::config::parse_cqmom_prepared;

CqmomLeanCaseLoader cqmom_lean_case_loader = qmomkit::config::load_cqmom_lean;

CqmomPreparedCaseLoader cqmom_prepared_case_loader = qmomkit::config::load_cqmom_prepared;

CqmomInitialMomentReplacement cqmom_initial_moment_replacement = qmomkit::config::with_initial_moments;

CqmomLeanInitialMomentReplacement cqmom_lean_initial_moment_replacement = qmomkit::config::with_initial_moments;

} // namespace public_header_pins
