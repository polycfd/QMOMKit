#include "qmomkit/output.hpp"

#include "internal/cqmom_validation.hpp"
#include "internal/numeric_utils.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace qmomkit
{
namespace
{

[[nodiscard]] std::string qmom_label(const std::size_t order)
{
    return "M[" + std::to_string(order) + "]";
}

[[nodiscard]] std::string cqmom_label(const MultiIndex& index)
{
    std::ostringstream label;
    label << "M[";
    for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
    {
        if (dimension != 0)
        {
            label << ',';
        }
        label << index[dimension];
    }
    label << ']';
    return label.str();
}

[[nodiscard]] std::vector<std::size_t> sequential_indices(const std::size_t count)
{
    if (count == 0)
    {
        throw std::invalid_argument("Moment output requires at least one moment.");
    }

    std::vector<std::size_t> indices;
    indices.reserve(count);
    for (std::size_t index = 0; index < count; ++index)
    {
        indices.push_back(index);
    }
    return indices;
}

[[nodiscard]] std::vector<std::string> qmom_labels(const std::vector<std::size_t>& orders)
{
    if (orders.empty())
    {
        throw std::invalid_argument("Selected QMOM moment output requires at least one moment order.");
    }

    std::vector<std::string> labels;
    labels.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        labels.push_back(qmom_label(order));
    }
    return labels;
}

[[nodiscard]] std::vector<std::string> cqmom_labels(const std::vector<MultiIndex>& indices)
{
    if (indices.empty())
    {
        throw std::invalid_argument("Selected CQMOM moment output requires at least one moment index.");
    }

    std::vector<std::string> labels;
    labels.reserve(indices.size());
    for (const MultiIndex& index : indices)
    {
        labels.push_back(cqmom_label(index));
    }
    return labels;
}

[[nodiscard]] std::vector<std::size_t> flatten_selected_indices(const std::vector<std::size_t>& extents, const std::vector<MultiIndex>& indices)
{
    std::vector<std::size_t> flat_indices;
    flat_indices.reserve(indices.size());
    for (const MultiIndex& index : indices)
    {
        flat_indices.push_back(flatten_multi_index(index, extents));
    }
    return flat_indices;
}

[[nodiscard]] std::vector<std::string> all_cqmom_labels(const std::vector<std::size_t>& extents)
{
    return cqmom_labels(enumerate_multi_indices(extents));
}

[[nodiscard]] std::string weight_label(const std::size_t node)
{
    return "w" + std::to_string(node);
}

[[nodiscard]] std::string qmom_abscissa_label(const std::size_t node)
{
    return "x" + std::to_string(node);
}

[[nodiscard]] std::string cqmom_abscissa_label(const std::size_t dimension, const std::size_t node)
{
    return "x" + std::to_string(dimension) + "_" + std::to_string(node);
}

[[nodiscard]] std::vector<std::string> qmom_quadrature_labels(const std::size_t node_count)
{
    if (node_count == 0)
    {
        throw std::invalid_argument("QMOM quadrature output requires at least one node.");
    }

    std::vector<std::string> labels;
    labels.reserve(2 * node_count);
    for (std::size_t node = 0; node < node_count; ++node)
    {
        labels.push_back(weight_label(node));
        labels.push_back(qmom_abscissa_label(node));
    }
    return labels;
}

[[nodiscard]] std::vector<std::string> cqmom_quadrature_labels(const std::size_t dimension_count, const std::size_t node_count)
{
    if (dimension_count == 0)
    {
        throw std::invalid_argument("CQMOM quadrature output requires at least one phase-space dimension.");
    }
    if (node_count == 0)
    {
        throw std::invalid_argument("CQMOM quadrature output requires at least one terminal node.");
    }

    std::vector<std::string> labels;
    labels.reserve(node_count * (dimension_count + 1));
    for (std::size_t node = 0; node < node_count; ++node)
    {
        labels.push_back(weight_label(node));
        for (std::size_t dimension = 0; dimension < dimension_count; ++dimension)
        {
            labels.push_back(cqmom_abscissa_label(dimension, node));
        }
    }
    return labels;
}

} // namespace

MomentWriter MomentWriter::for_qmom(const std::filesystem::path& path, const std::size_t moment_count)
{
    const std::vector<std::size_t> indices = sequential_indices(moment_count);
    return MomentWriter(path, MomentLayout::qmom, {moment_count}, indices, qmom_labels(indices));
}

MomentWriter MomentWriter::for_qmom(const std::filesystem::path& path, std::vector<std::size_t> selected_orders)
{
    std::vector<std::size_t> extents;
    return MomentWriter(path, MomentLayout::qmom, std::move(extents), selected_orders, qmom_labels(selected_orders));
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents)
{
    const std::vector<std::size_t> indices = sequential_indices(product_extents(extents));
    std::vector<std::string> labels = all_cqmom_labels(extents);
    return MomentWriter(path, MomentLayout::cqmom_dense, std::move(extents), indices, std::move(labels));
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents, std::vector<MultiIndex> selected_moments)
{
    const std::vector<std::size_t> flat_indices = flatten_selected_indices(extents, selected_moments);
    std::vector<std::string> labels = cqmom_labels(selected_moments);
    return MomentWriter(path, MomentLayout::cqmom_dense, std::move(extents), flat_indices, std::move(labels));
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, const CqmomMomentLayout& layout)
{
    return MomentWriter(path, MomentLayout::cqmom_moment_data, layout, cqmom_labels(layout.required_indices));
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, const RawMomentRequest& request)
{
    return MomentWriter(path, MomentLayout::cqmom_requested, request, cqmom_labels(request.indices));
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, std::vector<std::size_t> extents,
                           std::vector<std::size_t> flat_indices, std::vector<std::string> labels)
    : output_(path), layout_(layout), extents_(std::move(extents)), flat_indices_(std::move(flat_indices)), labels_(std::move(labels))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open moment output file: " + path.string());
    }
    if (flat_indices_.empty() || flat_indices_.size() != labels_.size())
    {
        throw std::invalid_argument("Moment output columns are inconsistent.");
    }
    write_header();
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, CqmomMomentLayout cqmom_layout,
                           std::vector<std::string> labels)
    : output_(path), layout_(layout), cqmom_layout_(std::move(cqmom_layout)), labels_(std::move(labels))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open moment output file: " + path.string());
    }
    if (layout_ != MomentLayout::cqmom_moment_data || !cqmom_layout_.has_value() || cqmom_layout_->required_indices.size() != labels_.size())
    {
        throw std::invalid_argument("CQMOM moment data output columns are inconsistent.");
    }
    write_header();
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, RawMomentRequest request, std::vector<std::string> labels)
    : output_(path), layout_(layout), raw_moment_request_(std::move(request)), labels_(std::move(labels))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open moment output file: " + path.string());
    }
    if (layout_ != MomentLayout::cqmom_requested || raw_moment_request_.indices.size() != labels_.size())
    {
        throw std::invalid_argument("Requested CQMOM moment output columns are inconsistent.");
    }
    write_header();
}

void MomentWriter::write_header()
{
    output_ << "# time";
    for (const std::string& label : labels_)
    {
        output_ << ' ' << label;
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing moment output header.");
    }
}

void MomentWriter::write_row(const double time, const std::vector<double>& values)
{
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("Moment output time must be finite.");
    }
    if (!internal::all_finite(values))
    {
        throw std::invalid_argument("Moment output values must be finite.");
    }

    output_ << std::setprecision(17) << time;
    for (const std::size_t flat_index : flat_indices_)
    {
        if (flat_index >= values.size())
        {
            throw std::out_of_range("Selected moment output index is outside the supplied moments.");
        }
        output_ << ' ' << std::setprecision(17) << values[flat_index];
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing moment output row.");
    }
}

void MomentWriter::write_aligned_row(const double time, const std::vector<double>& values)
{
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("Moment output time must be finite.");
    }
    if (!internal::all_finite(values))
    {
        throw std::invalid_argument("Moment output values must be finite.");
    }
    if (values.size() != labels_.size())
    {
        throw std::invalid_argument("Moment output row does not match the writer column count.");
    }

    output_ << std::setprecision(17) << time;
    for (const double value : values)
    {
        output_ << ' ' << std::setprecision(17) << value;
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing moment output row.");
    }
}

void MomentWriter::write(const double time, const MomentSet& moments)
{
    if (layout_ != MomentLayout::qmom)
    {
        throw std::invalid_argument("QMOM moment output cannot be written with a CQMOM moment writer.");
    }
    if (!extents_.empty() && moments.size() != extents_.front())
    {
        throw std::invalid_argument("QMOM moment output row does not match the writer moment count.");
    }
    write_row(time, moments.values());
}

void MomentWriter::write(const double time, const RawMomentTensor& moments)
{
    if (layout_ != MomentLayout::cqmom_dense)
    {
        throw std::invalid_argument("Dense CQMOM moment output cannot be written with this moment writer.");
    }
    if (moments.extents != extents_)
    {
        throw std::invalid_argument("CQMOM moment output tensor extents do not match the writer extents.");
    }
    if (moments.values.size() != product_extents(moments.extents))
    {
        throw std::invalid_argument("CQMOM moment output tensor size does not match its extents.");
    }
    write_row(time, moments.values);
}

void MomentWriter::write(const double time, const CqmomMomentData& moments)
{
    if (layout_ != MomentLayout::cqmom_moment_data)
    {
        throw std::invalid_argument("CQMOM moment data output cannot be written with this moment writer.");
    }
    if (!cqmom_layout_.has_value() || !internal::same_cqmom_layout(moments.layout, *cqmom_layout_))
    {
        throw std::invalid_argument("CQMOM moment data output layout does not match the writer layout.");
    }
    write_aligned_row(time, moments.values);
}

void MomentWriter::write(const double time, const CqmomQuadrature& quadrature)
{
    if (layout_ != MomentLayout::cqmom_requested)
    {
        throw std::invalid_argument("Requested CQMOM moment output cannot be written with this moment writer.");
    }
    const SelectedRawMomentResult reconstructed = reconstruct_raw_moments(quadrature, raw_moment_request_);
    if (!reconstructed.ok())
    {
        throw std::invalid_argument("Requested CQMOM moment reconstruction failed before output: " + reconstructed.message);
    }
    write_aligned_row(time, reconstructed.values);
}

QuadratureWriter QuadratureWriter::for_qmom(const std::filesystem::path& path, const std::size_t node_count)
{
    return QuadratureWriter(path, QuadratureLayout::qmom, 1, node_count, qmom_quadrature_labels(node_count));
}

QuadratureWriter QuadratureWriter::for_cqmom(const std::filesystem::path& path, const std::size_t dimension_count, const std::size_t node_count)
{
    return QuadratureWriter(path, QuadratureLayout::cqmom, dimension_count, node_count, cqmom_quadrature_labels(dimension_count, node_count));
}

QuadratureWriter::QuadratureWriter(const std::filesystem::path& path, const QuadratureLayout layout, const std::size_t dimension_count,
                                   const std::size_t node_count, std::vector<std::string> labels)
    : output_(path), layout_(layout), dimension_count_(dimension_count), node_count_(node_count), labels_(std::move(labels))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open quadrature output file: " + path.string());
    }
    if (labels_.empty())
    {
        throw std::invalid_argument("Quadrature output requires at least one column.");
    }
    write_header();
}

void QuadratureWriter::write_header()
{
    output_ << "# time";
    for (const std::string& label : labels_)
    {
        output_ << ' ' << label;
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing quadrature output header.");
    }
}

void QuadratureWriter::write_row(const double time, const std::vector<double>& values)
{
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("Quadrature output time must be finite.");
    }
    if (!internal::all_finite(values))
    {
        throw std::invalid_argument("Quadrature output values must be finite.");
    }
    if (values.size() != labels_.size())
    {
        throw std::invalid_argument("Quadrature output row does not match the writer column count.");
    }

    output_ << std::setprecision(17) << time;
    for (const double value : values)
    {
        output_ << ' ' << std::setprecision(17) << value;
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing quadrature output row.");
    }
}

void QuadratureWriter::write_qmom_row(const double time, const std::vector<double>& weights, const std::vector<double>& abscissas)
{
    if (layout_ != QuadratureLayout::qmom)
    {
        throw std::invalid_argument("QMOM quadrature output cannot be written with a CQMOM quadrature writer.");
    }
    if (weights.size() != node_count_ || abscissas.size() != node_count_)
    {
        throw std::invalid_argument("QMOM quadrature output row does not match the writer node count.");
    }

    std::vector<double> values;
    values.reserve(2 * node_count_);
    for (std::size_t node = 0; node < node_count_; ++node)
    {
        values.push_back(weights[node]);
        values.push_back(abscissas[node]);
    }
    write_row(time, values);
}

void QuadratureWriter::write(const double time, const Quadrature1D& quadrature)
{
    write_qmom_row(time, quadrature.weights, quadrature.abscissas);
}

void QuadratureWriter::write(const double time, const QmomResult& quadrature)
{
    write_qmom_row(time, quadrature.weights, quadrature.abscissas);
}

void QuadratureWriter::write(const double time, const CqmomQuadrature& quadrature)
{
    if (layout_ != QuadratureLayout::cqmom)
    {
        throw std::invalid_argument("CQMOM quadrature output cannot be written with a QMOM quadrature writer.");
    }
    if (quadrature.number_of_dimensions() != dimension_count_)
    {
        throw std::invalid_argument("CQMOM quadrature output dimension count does not match the writer dimension count.");
    }
    if (quadrature.total_nodes() != node_count_ || quadrature.weights().size() != node_count_)
    {
        throw std::invalid_argument("CQMOM quadrature output node count does not match the writer node count.");
    }
    if (quadrature.abscissas().size() != dimension_count_)
    {
        throw std::invalid_argument("CQMOM quadrature output abscissa dimension count is inconsistent.");
    }

    std::vector<double> values;
    values.reserve(node_count_ * (dimension_count_ + 1));
    for (std::size_t node = 0; node < node_count_; ++node)
    {
        values.push_back(quadrature.weights()[node]);
        for (std::size_t dimension = 0; dimension < dimension_count_; ++dimension)
        {
            if (quadrature.abscissas()[dimension].size() != node_count_)
            {
                throw std::invalid_argument("CQMOM quadrature output abscissa node count is inconsistent.");
            }
            values.push_back(quadrature.abscissas()[dimension][node]);
        }
    }
    write_row(time, values);
}

} // namespace qmomkit
