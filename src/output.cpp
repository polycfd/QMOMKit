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

[[nodiscard]] std::string centered_column_label(const std::string& label, const std::size_t width)
{
    if (label.size() >= width)
    {
        return label;
    }

    const std::size_t padding = width - label.size();
    const std::size_t left_padding = padding / 2;
    return std::string(left_padding, ' ') + label + std::string(padding - left_padding, ' ');
}

[[nodiscard]] std::size_t column_width_for(const std::size_t precision)
{
    return std::max(precision + 8, std::size_t(6));
}

void write_column_header(std::ofstream& output, const std::size_t column_width, const std::vector<std::string>& labels, const char* context)
{
    output << centered_column_label("# time", column_width);
    for (const std::string& label : labels)
    {
        output << ' ' << centered_column_label(label, column_width);
    }
    output << '\n';
    if (!output)
    {
        throw std::runtime_error(std::string("Failed while writing ") + context + " header.");
    }
}

[[nodiscard]] std::string cqmom_label(const MultiIndexView index)
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

[[nodiscard]] std::string generalized_moment_label(const std::vector<double>& orders)
{
    std::ostringstream label;
    label << "G[" << std::setprecision(17);
    for (std::size_t dimension = 0; dimension < orders.size(); ++dimension)
    {
        if (dimension != 0)
        {
            label << ',';
        }
        label << orders[dimension];
    }
    label << ']';
    return label.str();
}

void append_generalized_labels(std::vector<std::string>& labels, const std::vector<std::vector<double>>& orders)
{
    labels.reserve(labels.size() + orders.size());
    for (const std::vector<double>& order : orders)
    {
        labels.push_back(generalized_moment_label(order));
    }
}

void validate_generalized_moment_orders(const std::vector<std::vector<double>>& orders, const std::size_t dimensions)
{
    for (const std::vector<double>& order : orders)
    {
        if (order.size() != dimensions ||
            !std::all_of(order.begin(), order.end(), [](const double value) { return std::isfinite(value) && value >= 0.0; }))
        {
            throw std::invalid_argument("Generalized output moment orders must be finite, non-negative, and match the phase-space dimension.");
        }
    }
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

[[nodiscard]] std::vector<std::string> cqmom_labels(const MultiIndexTable& indices)
{
    if (indices.empty())
    {
        throw std::invalid_argument("Selected CQMOM moment output requires at least one moment index.");
    }

    std::vector<std::string> labels;
    labels.reserve(indices.size());
    for (const MultiIndexView index : indices)
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
    return cqmom_labels(MultiIndexTable(enumerate_multi_indices(extents)));
}

[[nodiscard]] std::string cqmom_abscissa_label(const std::size_t dimension)
{
    return "x" + std::to_string(dimension);
}

[[nodiscard]] std::vector<std::string> qmom_quadrature_labels(const std::size_t maximum_node_count)
{
    if (maximum_node_count == 0)
    {
        throw std::invalid_argument("QMOM quadrature output requires at least one node.");
    }

    return {"node_index", "weight", "x0"};
}

[[nodiscard]] std::vector<std::string> cqmom_quadrature_labels(const std::size_t dimension_count)
{
    if (dimension_count == 0)
    {
        throw std::invalid_argument("CQMOM quadrature output requires at least one phase-space dimension.");
    }

    std::vector<std::string> labels;
    labels.reserve(2 * dimension_count + 2);
    labels.push_back("node_index");
    for (std::size_t level = 0; level < dimension_count; ++level)
    {
        labels.push_back("path_" + std::to_string(level));
    }
    labels.push_back("weight");
    for (std::size_t dimension = 0; dimension < dimension_count; ++dimension)
    {
        labels.push_back(cqmom_abscissa_label(dimension));
    }
    return labels;
}

} // namespace

MomentWriter MomentWriter::for_qmom(const std::filesystem::path& path, const std::size_t moment_count, const std::size_t precision)
{
    const std::vector<std::size_t> indices = sequential_indices(moment_count);
    return MomentWriter(path, MomentLayout::qmom, {moment_count}, indices, qmom_labels(indices), {}, precision);
}

MomentWriter MomentWriter::for_qmom(const std::filesystem::path& path, std::vector<std::size_t> selected_orders, const std::size_t precision)
{
    std::vector<std::size_t> extents;
    return MomentWriter(path, MomentLayout::qmom, std::move(extents), selected_orders, qmom_labels(selected_orders), {}, precision);
}

MomentWriter MomentWriter::for_qmom(const std::filesystem::path& path, const std::size_t moment_count, const OutputOptions& output_options,
                                    const std::size_t precision)
{
    const std::vector<std::size_t> indices = sequential_indices(moment_count);
    validate_generalized_moment_orders(output_options.generalized_moment_orders, 1);
    std::vector<std::string> labels = qmom_labels(indices);
    append_generalized_labels(labels, output_options.generalized_moment_orders);
    return MomentWriter(path, MomentLayout::qmom, {moment_count}, indices, std::move(labels), output_options.generalized_moment_orders, precision);
}

MomentWriter MomentWriter::for_qmom(const std::filesystem::path& path, const QmomState& state, const OutputOptions& output_options,
                                    const std::size_t precision)
{
    const std::vector<std::size_t> indices = sequential_indices(state.moments().size());
    validate_generalized_moment_orders(output_options.generalized_moment_orders, 1);
    std::vector<std::string> labels = qmom_labels(indices);
    append_generalized_labels(labels, output_options.generalized_moment_orders);
    return MomentWriter(path, MomentLayout::qmom_state, {state.moments().size()}, indices, std::move(labels),
                        output_options.generalized_moment_orders, precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents, const std::size_t precision)
{
    const std::vector<std::size_t> indices = sequential_indices(product_extents(extents));
    std::vector<std::string> labels = all_cqmom_labels(extents);
    return MomentWriter(path, MomentLayout::cqmom_moment_tensor, std::move(extents), indices, std::move(labels), {}, precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents, std::vector<MultiIndex> selected_moments,
                                     const std::size_t precision)
{
    return for_cqmom(path, std::move(extents), std::move(selected_moments), {}, precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents, std::vector<MultiIndex> selected_moments,
                                     std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision)
{
    validate_generalized_moment_orders(generalized_moment_orders, extents.size());
    const std::vector<std::size_t> flat_indices = flatten_selected_indices(extents, selected_moments);
    std::vector<std::string> labels = cqmom_labels(MultiIndexTable(selected_moments));
    append_generalized_labels(labels, generalized_moment_orders);
    return MomentWriter(path, MomentLayout::cqmom_moment_tensor, std::move(extents), flat_indices, std::move(labels),
                        std::move(generalized_moment_orders), precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, const CqmomMomentLayout& layout, const std::size_t precision)
{
    return for_cqmom(path, layout, {}, precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, const CqmomMomentLayout& layout,
                                     std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision)
{
    validate_generalized_moment_orders(generalized_moment_orders, layout.hierarchy.number_of_dimensions());
    std::vector<std::string> labels = cqmom_labels(layout.required_indices);
    append_generalized_labels(labels, generalized_moment_orders);
    return MomentWriter(path, MomentLayout::cqmom_moment_data, layout, std::move(labels), std::move(generalized_moment_orders), precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, const CqmomState& state, const OutputOptions& output_options,
                                     const std::size_t precision)
{
    std::vector<MultiIndex> indices = state.moment_indices();
    std::vector<std::string> labels = cqmom_labels(MultiIndexTable(indices));
    validate_generalized_moment_orders(output_options.generalized_moment_orders, state.hierarchy().number_of_dimensions());
    append_generalized_labels(labels, output_options.generalized_moment_orders);
    return MomentWriter(path, MomentLayout::cqmom_state, std::move(indices), std::move(labels), output_options.generalized_moment_orders, precision);
}

MomentWriter MomentWriter::for_cqmom(const std::filesystem::path& path, const RawMomentRequest& request, const std::size_t precision)
{
    return MomentWriter(path, MomentLayout::cqmom_requested, request, cqmom_labels(MultiIndexTable(request.indices)), precision);
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, std::vector<std::size_t> extents,
                           std::vector<std::size_t> flat_indices, std::vector<std::string> labels,
                           std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision)
    : output_(path), layout_(layout), extents_(std::move(extents)), flat_indices_(std::move(flat_indices)),
      generalized_moment_orders_(std::move(generalized_moment_orders)), labels_(std::move(labels)), precision_(precision),
      column_width_(column_width_for(precision))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open moment output file: " + path.string());
    }
    if (flat_indices_.empty() || flat_indices_.size() + generalized_moment_orders_.size() != labels_.size())
    {
        throw std::invalid_argument("Moment output columns are inconsistent.");
    }
    write_header();
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, CqmomMomentLayout cqmom_layout,
                           std::vector<std::string> labels, std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision)
    : output_(path), layout_(layout), cqmom_layout_(std::move(cqmom_layout)), generalized_moment_orders_(std::move(generalized_moment_orders)),
      labels_(std::move(labels)), precision_(precision), column_width_(column_width_for(precision))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open moment output file: " + path.string());
    }
    if (layout_ != MomentLayout::cqmom_moment_data || !cqmom_layout_.has_value() ||
        cqmom_layout_->required_indices.size() + generalized_moment_orders_.size() != labels_.size())
    {
        throw std::invalid_argument("CQMOM moment data output columns are inconsistent.");
    }
    write_header();
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, std::vector<MultiIndex> state_indices,
                           std::vector<std::string> labels, std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision)
    : output_(path), layout_(layout), state_indices_(std::move(state_indices)), generalized_moment_orders_(std::move(generalized_moment_orders)),
      labels_(std::move(labels)), precision_(precision), column_width_(column_width_for(precision))
{
    if (!output_)
    {
        throw std::runtime_error("Could not open moment output file: " + path.string());
    }
    if (layout_ != MomentLayout::cqmom_state || state_indices_.empty() || state_indices_.size() + generalized_moment_orders_.size() != labels_.size())
    {
        throw std::invalid_argument("CQMOM state moment output columns are inconsistent.");
    }
    write_header();
}

MomentWriter::MomentWriter(const std::filesystem::path& path, const MomentLayout layout, RawMomentRequest request, std::vector<std::string> labels,
                           const std::size_t precision)
    : output_(path), layout_(layout), raw_moment_request_(std::move(request)), labels_(std::move(labels)), precision_(precision),
      column_width_(column_width_for(precision))
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

std::string MomentWriter::layout_mismatch_message(const char* const overload, const MomentLayout required, const MomentLayout actual)
{
    const auto factory = [](const MomentLayout layout) -> std::string
    {
        switch (layout)
        {
        case MomentLayout::qmom:
            return "MomentWriter::for_qmom(path, moment_count[, output_options]) or MomentWriter::for_qmom(path, selected_orders)";
        case MomentLayout::qmom_state:
            return "MomentWriter::for_qmom(path, qmom_state, output_options)";
        case MomentLayout::cqmom_moment_tensor:
            return "MomentWriter::for_cqmom(path, extents, ...)";
        case MomentLayout::cqmom_moment_data:
            return "MomentWriter::for_cqmom(path, cqmom_moment_layout, ...)";
        case MomentLayout::cqmom_state:
            return "MomentWriter::for_cqmom(path, cqmom_state, output_options)";
        case MomentLayout::cqmom_requested:
            return "MomentWriter::for_cqmom(path, raw_moment_request)";
        }
        return "an unknown factory";
    };

    return "MomentWriter::" + std::string(overload) + " requires a writer created by " + factory(required) + ", but this writer was created by " +
           factory(actual) + ".";
}

void MomentWriter::write_header()
{
    write_column_header(output_, column_width_, labels_, "moment output");
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

    output_ << std::scientific << std::setprecision(precision_) << std::setw(column_width_) << time;
    for (const std::size_t flat_index : flat_indices_)
    {
        if (flat_index >= values.size())
        {
            throw std::out_of_range("Selected moment output index is outside the supplied moments.");
        }
        output_ << ' ' << std::setprecision(precision_) << std::setw(column_width_) << values[flat_index];
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing moment output row.");
    }
}

void MomentWriter::write_aligned_row(const double time, const std::vector<double>& values, const std::vector<double>& generalized_values)
{
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("Moment output time must be finite.");
    }
    if (!internal::all_finite(values) || !internal::all_finite(generalized_values))
    {
        throw std::invalid_argument("Moment output values must be finite.");
    }
    if (values.size() + generalized_values.size() != labels_.size())
    {
        throw std::invalid_argument("Moment output row does not match the writer column count.");
    }

    output_ << std::scientific << std::setprecision(precision_) << std::setw(column_width_) << time;
    for (const double value : values)
    {
        output_ << ' ' << std::setprecision(precision_) << std::setw(column_width_) << value;
    }
    for (const double value : generalized_values)
    {
        output_ << ' ' << std::setprecision(precision_) << std::setw(column_width_) << value;
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
        throw std::invalid_argument(layout_mismatch_message("write(time, MomentSet)", MomentLayout::qmom, layout_));
    }
    if (!extents_.empty() && moments.size() != extents_.front())
    {
        throw std::invalid_argument("QMOM moment output row does not match the writer moment count.");
    }
    if (!generalized_moment_orders_.empty())
    {
        throw std::invalid_argument(
            "MomentWriter::write(time, MomentSet) cannot supply the writer's generalized-moment columns; use write(time, MomentSet, "
            "generalized_values) or write(time, MomentSet, Quadrature1D).");
    }
    write_row(time, moments.values());
}

void MomentWriter::write(const double time, const MomentSet& moments, const std::vector<double>& generalized_values)
{
    if (layout_ != MomentLayout::qmom)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, MomentSet, generalized_values)", MomentLayout::qmom, layout_));
    }
    if (!extents_.empty() && moments.size() != extents_.front())
    {
        throw std::invalid_argument("QMOM moment output row does not match the writer moment count.");
    }
    std::vector<double> values;
    values.reserve(flat_indices_.size());
    for (const std::size_t index : flat_indices_)
    {
        if (index >= moments.size())
        {
            throw std::out_of_range("Selected moment output index is outside the supplied moments.");
        }
        values.push_back(moments[index]);
    }
    write_aligned_row(time, values, generalized_values);
}

void MomentWriter::write(const double time, const MomentSet& moments, const Quadrature1D& quadrature)
{
    if (layout_ != MomentLayout::qmom)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, MomentSet, Quadrature1D)", MomentLayout::qmom, layout_));
    }
    if (!extents_.empty() && moments.size() != extents_.front())
    {
        throw std::invalid_argument("QMOM moment output row does not match the writer moment count.");
    }
    std::vector<double> generalized_values;
    generalized_values.reserve(generalized_moment_orders_.size());
    for (const std::vector<double>& orders : generalized_moment_orders_)
    {
        generalized_values.push_back(evaluate_generalized_moment(quadrature, orders.front()));
    }
    write(time, moments, generalized_values);
}

void MomentWriter::write(const double time, const QmomState& state, const QmomResult& quadrature)
{
    if (layout_ != MomentLayout::qmom_state)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, QmomState, QmomResult)", MomentLayout::qmom_state, layout_));
    }
    if (extents_.size() != 1 || state.moments().size() != extents_.front())
    {
        throw std::invalid_argument("QMOM state moment count does not match the writer columns.");
    }
    std::vector<double> generalized_values;
    generalized_values.reserve(generalized_moment_orders_.size());
    for (const std::vector<double>& orders : generalized_moment_orders_)
    {
        generalized_values.push_back(evaluate_generalized_moment(quadrature, orders.front()));
    }
    write_aligned_row(time, state.moments(), generalized_values);
}

void MomentWriter::write(const double time, const RawMomentTensor& moment_tensor)
{
    if (layout_ != MomentLayout::cqmom_moment_tensor)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, RawMomentTensor)", MomentLayout::cqmom_moment_tensor, layout_));
    }
    if (moment_tensor.extents != extents_)
    {
        throw std::invalid_argument("CQMOM moment output tensor extents do not match the writer extents.");
    }
    if (moment_tensor.values.size() != product_extents(moment_tensor.extents))
    {
        throw std::invalid_argument("CQMOM moment output tensor size does not match its extents.");
    }
    if (!generalized_moment_orders_.empty())
    {
        throw std::invalid_argument(
            "MomentWriter::write(time, RawMomentTensor) cannot supply the writer's generalized-moment columns; use write(time, "
            "RawMomentTensor, generalized_values).");
    }
    write_row(time, moment_tensor.values);
}

void MomentWriter::write(const double time, const RawMomentTensor& moment_tensor, const std::vector<double>& generalized_values)
{
    if (layout_ != MomentLayout::cqmom_moment_tensor)
    {
        throw std::invalid_argument(
            layout_mismatch_message("write(time, RawMomentTensor, generalized_values)", MomentLayout::cqmom_moment_tensor, layout_));
    }
    if (moment_tensor.extents != extents_)
    {
        throw std::invalid_argument("CQMOM moment output tensor extents do not match the writer extents.");
    }
    if (moment_tensor.values.size() != product_extents(moment_tensor.extents))
    {
        throw std::invalid_argument("CQMOM moment output tensor size does not match its extents.");
    }
    std::vector<double> values;
    values.reserve(flat_indices_.size());
    for (const std::size_t flat_index : flat_indices_)
    {
        if (flat_index >= moment_tensor.values.size())
        {
            throw std::out_of_range("Selected moment output index is outside the supplied moments.");
        }
        values.push_back(moment_tensor.values[flat_index]);
    }
    write_aligned_row(time, values, generalized_values);
}

void MomentWriter::write(const double time, const CqmomMomentData& moments)
{
    if (layout_ != MomentLayout::cqmom_moment_data)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, CqmomMomentData)", MomentLayout::cqmom_moment_data, layout_));
    }
    if (!cqmom_layout_.has_value() || !internal::same_cqmom_layout(moments.layout, *cqmom_layout_))
    {
        throw std::invalid_argument("CQMOM moment data output layout does not match the writer layout.");
    }
    if (!generalized_moment_orders_.empty())
    {
        throw std::invalid_argument(
            "MomentWriter::write(time, CqmomMomentData) cannot supply the writer's generalized-moment columns; use write(time, "
            "CqmomMomentData, generalized_values).");
    }
    write_aligned_row(time, moments.values);
}

void MomentWriter::write(const double time, const CqmomMomentData& moments, const std::vector<double>& generalized_values)
{
    if (layout_ != MomentLayout::cqmom_moment_data)
    {
        throw std::invalid_argument(
            layout_mismatch_message("write(time, CqmomMomentData, generalized_values)", MomentLayout::cqmom_moment_data, layout_));
    }
    if (!cqmom_layout_.has_value() || !internal::same_cqmom_layout(moments.layout, *cqmom_layout_))
    {
        throw std::invalid_argument("CQMOM moment data output layout does not match the writer layout.");
    }
    write_aligned_row(time, moments.values, generalized_values);
}

void MomentWriter::write(const double time, const CqmomState& state, const CqmomQuadrature& quadrature)
{
    if (layout_ != MomentLayout::cqmom_state)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, CqmomState, CqmomQuadrature)", MomentLayout::cqmom_state, layout_));
    }
    if (state.moment_indices() != state_indices_)
    {
        throw std::invalid_argument("CQMOM state moment indices do not match the writer columns.");
    }
    std::vector<double> generalized_values;
    generalized_values.reserve(generalized_moment_orders_.size());
    for (const std::vector<double>& orders : generalized_moment_orders_)
    {
        generalized_values.push_back(evaluate_generalized_moment(quadrature, orders));
    }
    write_aligned_row(time, state.moment_values(), generalized_values);
}

void MomentWriter::write(const double time, const CqmomQuadrature& quadrature)
{
    if (layout_ != MomentLayout::cqmom_requested)
    {
        throw std::invalid_argument(layout_mismatch_message("write(time, CqmomQuadrature)", MomentLayout::cqmom_requested, layout_));
    }
    const SelectedRawMomentResult reconstructed = reconstruct_raw_moments(quadrature, raw_moment_request_);
    if (!reconstructed.ok())
    {
        throw std::invalid_argument("Requested CQMOM moment reconstruction failed before output: " + reconstructed.message);
    }
    write_aligned_row(time, reconstructed.values);
}

QuadratureWriter QuadratureWriter::for_qmom(const std::filesystem::path& path, const std::size_t maximum_node_count, const std::size_t precision)
{
    return QuadratureWriter(path, QuadratureLayout::qmom, 1, maximum_node_count, qmom_quadrature_labels(maximum_node_count), precision);
}

QuadratureWriter QuadratureWriter::for_cqmom(const std::filesystem::path& path, const std::size_t dimension_count, const std::size_t precision)
{
    return QuadratureWriter(path, QuadratureLayout::cqmom, dimension_count, 0, cqmom_quadrature_labels(dimension_count), precision);
}

QuadratureWriter::QuadratureWriter(const std::filesystem::path& path, const QuadratureLayout layout, const std::size_t dimension_count,
                                   const std::size_t maximum_node_count, std::vector<std::string> labels, const std::size_t precision)
    : output_(path), layout_(layout), dimension_count_(dimension_count), maximum_node_count_(maximum_node_count), labels_(std::move(labels)),
      precision_(precision), column_width_(column_width_for(precision))
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
    write_column_header(output_, column_width_, labels_, "quadrature output");
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

    output_ << std::scientific << std::setprecision(precision_) << std::setw(column_width_) << time;
    for (const double value : values)
    {
        output_ << ' ' << std::setprecision(precision_) << std::setw(column_width_) << value;
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
    if (weights.empty())
    {
        throw std::invalid_argument("QMOM quadrature output requires at least one active node.");
    }
    if (weights.size() != abscissas.size())
    {
        throw std::invalid_argument("QMOM quadrature output weight and abscissa counts must match.");
    }
    if (weights.size() > maximum_node_count_)
    {
        throw std::invalid_argument("QMOM quadrature output active node count exceeds the writer maximum.");
    }

    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        std::vector<double> values;
        values.reserve(3);
        values.push_back(static_cast<double>(node));
        values.push_back(weights[node]);
        values.push_back(abscissas[node]);
        write_row(time, values);
    }
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
    const std::size_t node_count = quadrature.weights().size();
    if (node_count == 0)
    {
        throw std::invalid_argument("CQMOM quadrature output requires at least one active terminal node.");
    }
    if (quadrature.abscissas().size() != dimension_count_)
    {
        throw std::invalid_argument("CQMOM quadrature output abscissa dimension count is inconsistent.");
    }
    if (quadrature.total_nodes() != node_count)
    {
        throw std::invalid_argument("CQMOM quadrature output terminal-node storage is inconsistent.");
    }
    if (quadrature.terminal_paths().size() != node_count)
    {
        throw std::invalid_argument("CQMOM quadrature output requires one terminal path per active terminal node.");
    }
    if (!internal::all_finite(quadrature.weights()))
    {
        throw std::invalid_argument("CQMOM quadrature output weights must be finite.");
    }
    for (const std::vector<double>& dimension_abscissas : quadrature.abscissas())
    {
        if (dimension_abscissas.size() != node_count)
        {
            throw std::invalid_argument("CQMOM quadrature output abscissa node count is inconsistent.");
        }
        if (!internal::all_finite(dimension_abscissas))
        {
            throw std::invalid_argument("CQMOM quadrature output abscissas must be finite.");
        }
    }

    for (std::size_t node = 0; node < node_count; ++node)
    {
        const ConditioningPathView path = quadrature.terminal_paths()[node];
        if (path.size() != dimension_count_)
        {
            throw std::invalid_argument("CQMOM quadrature output terminal path depth does not match the writer dimension count.");
        }
        // Discards the computed index; called only for its out_of_range check that every
        // conditioning-path node index is within the hierarchy's bounds at each level.
        (void)recursive_path_to_node_index(quadrature.hierarchy(), path);

        std::vector<double> values;
        values.reserve(2 * dimension_count_ + 2);
        values.push_back(static_cast<double>(node));
        for (std::size_t level = 0; level < dimension_count_; ++level)
        {
            values.push_back(static_cast<double>(path.node_index_at_level(level)));
        }
        values.push_back(quadrature.weights()[node]);
        for (std::size_t dimension = 0; dimension < dimension_count_; ++dimension)
        {
            values.push_back(quadrature.abscissas()[dimension][node]);
        }
        write_row(time, values);
    }
}

} // namespace qmomkit
