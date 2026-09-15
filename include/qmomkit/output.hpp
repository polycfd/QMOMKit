#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

/**
 * Optional application-controlled output configuration.
 *
 * The numerical algorithms do not consume this structure directly. By default,
 * QMOMKit performs no file output and no terminal output.
 */
struct OutputOptions
{
    std::filesystem::path directory = ".";
    bool write_moments = false;
    bool write_quadrature = false;
    std::optional<double> interval;
    std::vector<std::vector<double>> generalized_moment_orders;
};

/**
 * Write raw moment histories as space-separated text.
 *
 * The writer is application-controlled: construction opens exactly one file,
 * writes a header, and each `write(...)` call appends one time row. The writer
 * does not schedule output and is not used by numerical algorithms.
 *
 * The factory fixes the column schema, so it also fixes which `write(...)`
 * overloads the writer accepts. The pairing is one-to-one in the data argument:
 * `MomentSet` rows belong to the `for_qmom(path, moment_count)` and
 * `for_qmom(path, selected_orders)` writers, `QmomState` rows to
 * `for_qmom(path, state, output_options)`, `RawMomentTensor` rows to the
 * `for_cqmom(path, extents, ...)` writers, `CqmomMomentData` rows to the
 * `for_cqmom(path, layout, ...)` writers, `CqmomState` rows to
 * `for_cqmom(path, state, output_options)`, and a bare `CqmomQuadrature` row to
 * `for_cqmom(path, request)`. Passing a row to a writer built by a different
 * factory throws `std::invalid_argument` naming the factory the caller needs.
 * `docs/developer/output_design.md` tabulates the same pairing.
 */
class MomentWriter
{
public:
    MomentWriter(const MomentWriter&) = delete;
    MomentWriter& operator=(const MomentWriter&) = delete;
    MomentWriter(MomentWriter&&) noexcept = default;
    MomentWriter& operator=(MomentWriter&&) noexcept = default;
    ~MomentWriter() = default;

    /** Write all one-dimensional moments M[0], M[1], ... M[moment_count - 1]. */
    [[nodiscard]] static MomentWriter for_qmom(const std::filesystem::path& path, std::size_t moment_count, const std::size_t precision = 17);

    /** Write selected one-dimensional moments in the supplied order. */
    [[nodiscard]] static MomentWriter for_qmom(const std::filesystem::path& path, std::vector<std::size_t> selected_orders,
                                               const std::size_t precision = 17);

    /** Advanced raw-moment output with configured generalized columns. */
    [[nodiscard]] static MomentWriter for_qmom(const std::filesystem::path& path, std::size_t moment_count, const OutputOptions& output_options,
                                               const std::size_t precision = 17);

    /** Write the moments represented by a QMOM state and configured generalized columns. */
    [[nodiscard]] static MomentWriter for_qmom(const std::filesystem::path& path, const QmomState& state, const OutputOptions& output_options,
                                               const std::size_t precision = 17);

    /** Write all tensor moments in row-major multi-index order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                const std::size_t precision = 17);

    /** Write selected tensor moments in the supplied order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                std::vector<MultiIndex> selected_moments, const std::size_t precision = 17);

    /** Write selected tensor moments followed by supplied generalized-moment columns. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                std::vector<MultiIndex> selected_moments, std::vector<std::vector<double>> generalized_moment_orders,
                                                const std::size_t precision = 17);

    /** Write current CQMOM moment data in the layout's required-index order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, const CqmomMomentLayout& layout, const std::size_t precision = 17);

    /** Write layout-ordered CQMOM moment data followed by supplied generalized-moment columns. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, const CqmomMomentLayout& layout,
                                                std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision = 17);

    /** Write the moments represented by a CQMOM state and configured generalized columns. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, const CqmomState& state, const OutputOptions& output_options,
                                                const std::size_t precision = 17);

    /** Write requested CQMOM moments reconstructed from quadrature in request order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, const RawMomentRequest& request, const std::size_t precision = 17);

    /** Append one row for a one-dimensional raw moment sequence. */
    void write(double time, const MomentSet& moments);

    /** Append selected one-dimensional moments followed by supplied generalized-moment values. */
    void write(double time, const MomentSet& moments, const std::vector<double>& generalized_values);

    /** Append raw moments followed by configured generalized moments from a quadrature rule. */
    void write(double time, const MomentSet& moments, const Quadrature1D& quadrature);

    /** Append state moments followed by configured generalized moments from the current quadrature. */
    void write(double time, const QmomState& state, const QmomResult& quadrature);

    /** Append one row for a multidimensional raw moment tensor. */
    void write(double time, const RawMomentTensor& moments);

    /** Append tensor moments followed by generalized-moment values in configured request order. */
    void write(double time, const RawMomentTensor& moments, const std::vector<double>& generalized_values);

    /** Append one row for CQMOM required moment data. */
    void write(double time, const CqmomMomentData& moments);

    /** Append layout-ordered moments followed by generalized-moment values in configured request order. */
    void write(double time, const CqmomMomentData& moments, const std::vector<double>& generalized_values);

    /** Append state moments followed by requested generalized moments from the current quadrature. */
    void write(double time, const CqmomState& state, const CqmomQuadrature& quadrature);

    /** Append one row for requested moments reconstructed from CQMOM quadrature. */
    void write(double time, const CqmomQuadrature& quadrature);

private:
    enum class MomentLayout
    {
        qmom,
        qmom_state,
        cqmom_dense,
        cqmom_moment_data,
        cqmom_state,
        cqmom_requested
    };

    MomentWriter(const std::filesystem::path& path, MomentLayout layout, std::vector<std::size_t> extents, std::vector<std::size_t> flat_indices,
                 std::vector<std::string> labels, std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision);
    MomentWriter(const std::filesystem::path& path, MomentLayout layout, CqmomMomentLayout cqmom_layout, std::vector<std::string> labels,
                 std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision);
    MomentWriter(const std::filesystem::path& path, MomentLayout layout, std::vector<MultiIndex> state_indices, std::vector<std::string> labels,
                 std::vector<std::vector<double>> generalized_moment_orders, const std::size_t precision);
    MomentWriter(const std::filesystem::path& path, MomentLayout layout, RawMomentRequest request, std::vector<std::string> labels,
                 const std::size_t precision);

    /** Report that `overload` needs a `required` writer but was called on an `actual` one, naming both factories. */
    [[nodiscard]] static std::string layout_mismatch_message(const char* overload, MomentLayout required, MomentLayout actual);

    void write_header();
    void write_row(double time, const std::vector<double>& values);
    void write_aligned_row(double time, const std::vector<double>& values, const std::vector<double>& generalized_values = {});

    std::ofstream output_;
    MomentLayout layout_ = MomentLayout::qmom;
    std::vector<std::size_t> extents_;
    std::vector<std::size_t> flat_indices_;
    std::optional<CqmomMomentLayout> cqmom_layout_;
    std::vector<MultiIndex> state_indices_;
    std::vector<std::vector<double>> generalized_moment_orders_;
    RawMomentRequest raw_moment_request_;
    std::vector<std::string> labels_;
    // Both are established in every constructor's member initializer list.
    std::size_t precision_;
    std::size_t column_width_;
};

/**
 * Write quadrature histories as space-separated text.
 *
 * The writer is application-controlled: construction opens exactly one file,
 * writes a header, and each `write(...)` call appends the full quadrature. QMOM
 * and CQMOM quadrature output use one row per active node.
 */
class QuadratureWriter
{
public:
    QuadratureWriter(const QuadratureWriter&) = delete;
    QuadratureWriter& operator=(const QuadratureWriter&) = delete;
    QuadratureWriter(QuadratureWriter&&) noexcept = default;
    QuadratureWriter& operator=(QuadratureWriter&&) noexcept = default;
    ~QuadratureWriter() = default;

    /** Write one-dimensional nodes with at most `maximum_node_count` active nodes. */
    [[nodiscard]] static QuadratureWriter for_qmom(const std::filesystem::path& path, std::size_t maximum_node_count,
                                                   const std::size_t precision = 17);

    /** Write terminal multidimensional nodes as one row per active terminal node. */
    [[nodiscard]] static QuadratureWriter for_cqmom(const std::filesystem::path& path, std::size_t dimension_count, const std::size_t precision = 17);

    /** Append one row for a one-dimensional quadrature rule. */
    void write(double time, const Quadrature1D& quadrature);

    /** Append one row for a public QMOM result. */
    void write(double time, const QmomResult& quadrature);

    /** Append one row per active terminal node for a public CQMOM quadrature. */
    void write(double time, const CqmomQuadrature& quadrature);

private:
    enum class QuadratureLayout
    {
        qmom,
        cqmom
    };

    QuadratureWriter(const std::filesystem::path& path, QuadratureLayout layout, std::size_t dimension_count, std::size_t maximum_node_count,
                     std::vector<std::string> labels, const std::size_t precision);

    void write_header();
    void write_row(double time, const std::vector<double>& values);
    void write_qmom_row(double time, const std::vector<double>& weights, const std::vector<double>& abscissas);

    std::ofstream output_;
    QuadratureLayout layout_ = QuadratureLayout::qmom;
    std::size_t dimension_count_ = 1;
    std::size_t maximum_node_count_ = 0;
    std::vector<std::string> labels_;
    // Both are established in every constructor's member initializer list.
    std::size_t precision_;
    std::size_t column_width_;
};

} // namespace qmomkit
