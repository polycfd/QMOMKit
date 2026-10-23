#include <qmomkit/moment_reconstruction.hpp>

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature, as
// in include_cqmom.cpp: a changed parameter or return type fails to compile, and
// because the pins have external linkage, a missing definition fails to link.

using QmomMomentReconstruction = std::vector<double> (*)(const qmomkit::Quadrature1D&, std::size_t);

using CqmomTensorReconstruction = qmomkit::RawMomentTensor (*)(const qmomkit::CqmomQuadrature&, const std::vector<std::size_t>&);

using CqmomSelectedReconstruction = qmomkit::SelectedRawMomentResult (*)(const qmomkit::CqmomQuadrature&, const qmomkit::RawMomentRequest&);

using CqmomLayoutReconstruction = qmomkit::CqmomMomentData (*)(const qmomkit::CqmomQuadrature&, const qmomkit::CqmomMomentLayout&);

using QuadratureGeneralizedMoment = double (*)(const qmomkit::Quadrature1D&, double);

using QmomGeneralizedMoment = double (*)(const qmomkit::QmomResult&, double);

using CqmomGeneralizedMoment = double (*)(const qmomkit::CqmomQuadrature&, const std::vector<double>&);

QmomMomentReconstruction qmom_moment_reconstruction = qmomkit::reconstruct_raw_moments;

CqmomTensorReconstruction cqmom_tensor_reconstruction = qmomkit::reconstruct_raw_moments;

CqmomSelectedReconstruction cqmom_selected_reconstruction = qmomkit::reconstruct_raw_moments;

CqmomLayoutReconstruction cqmom_layout_reconstruction = qmomkit::reconstruct_cqmom_moments;

QuadratureGeneralizedMoment quadrature_generalized_moment = qmomkit::evaluate_generalized_moment;

QmomGeneralizedMoment qmom_generalized_moment = qmomkit::evaluate_generalized_moment;

CqmomGeneralizedMoment cqmom_generalized_moment = qmomkit::evaluate_generalized_moment;

} // namespace public_header_pins
