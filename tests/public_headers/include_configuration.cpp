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

using CqmomMomentTensorInitializer =
    qmomkit::config::ConfigurationResult<qmomkit::RawMomentTensor> (*)(const qmomkit::config::InitialConditionSpecification&, std::string);

using CqmomMomentInitializer = qmomkit::config::ConfigurationResult<qmomkit::CqmomMomentData> (*)(
    const qmomkit::config::InitialConditionSpecification&, const qmomkit::PhaseSpaceHierarchy&, std::string);

using CqmomCaseFactory = qmomkit::config::ConfigurationResult<qmomkit::CqmomCase> (*)(const qmomkit::config::InitialConditionSpecification&,
                                                                                      const qmomkit::PhaseSpaceHierarchy&,
                                                                                      const qmomkit::QmomOptions&, std::string,
                                                                                      qmomkit::CqmomAbscissaBounds);

using CqmomCaseParser = qmomkit::config::ConfigurationResult<qmomkit::CqmomCase> (*)(std::string_view, std::string);

using CqmomPreparedCaseParser = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(std::string_view, std::string);

using CqmomCaseLoader = qmomkit::config::ConfigurationResult<qmomkit::CqmomCase> (*)(const std::filesystem::path&);

using CqmomPreparedCaseLoader = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const std::filesystem::path&);

using CqmomCaseMomentTensorReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomCase> (*)(const qmomkit::CqmomCase&,
                                                                                                      qmomkit::RawMomentTensor);

using CqmomCaseMomentDataReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomCase> (*)(const qmomkit::CqmomCase&,
                                                                                                    qmomkit::CqmomMomentData);

using CqmomPreparedCaseMomentTensorReplacement =
    qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const qmomkit::CqmomPreparedCase&, qmomkit::RawMomentTensor);

using CqmomPreparedCaseMomentDataReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const qmomkit::CqmomPreparedCase&,
                                                                                                                    qmomkit::CqmomMomentData);

QmomMomentInitializer qmom_moment_initializer = qmomkit::config::initialize_qmom_moments;

CqmomMomentTensorInitializer cqmom_moment_tensor_initializer = qmomkit::config::initialize_cqmom_moment_tensor;

CqmomMomentInitializer cqmom_moment_initializer = qmomkit::config::initialize_cqmom_moments;

CqmomCaseFactory cqmom_case_factory = qmomkit::config::create_cqmom_case;

CqmomCaseParser cqmom_case_parser = qmomkit::config::parse_cqmom;

CqmomPreparedCaseParser cqmom_prepared_case_parser = qmomkit::config::parse_cqmom_prepared;

CqmomCaseLoader cqmom_case_loader = qmomkit::config::load_cqmom;

CqmomPreparedCaseLoader cqmom_prepared_case_loader = qmomkit::config::load_cqmom_prepared;

CqmomCaseMomentTensorReplacement cqmom_case_moment_tensor_replacement = qmomkit::config::with_initial_moments;

CqmomCaseMomentDataReplacement cqmom_case_moment_data_replacement = qmomkit::config::with_initial_moments;

CqmomPreparedCaseMomentTensorReplacement cqmom_prepared_case_moment_tensor_replacement = qmomkit::config::with_initial_moments;

CqmomPreparedCaseMomentDataReplacement cqmom_prepared_case_moment_data_replacement = qmomkit::config::with_initial_moments;

} // namespace public_header_pins
