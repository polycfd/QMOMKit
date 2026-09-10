#pragma once

#include "qmomkit/conditional_moments.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"

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
};

/**
 * Write raw moment histories as space-separated text.
 *
 * The writer is application-controlled: construction opens exactly one file,
 * writes a header, and each `write(...)` call appends one time row. The writer
 * does not schedule output and is not used by numerical algorithms.
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

    /** Write all tensor moments in row-major multi-index order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                const std::size_t precision = 17);

    /** Write selected tensor moments in the supplied order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                std::vector<MultiIndex> selected_moments, const std::size_t precision = 17);

    /** Write current CQMOM moment data in the layout's required-index order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, const CqmomMomentLayout& layout, const std::size_t precision = 17);

    /** Write requested CQMOM moments reconstructed from quadrature in request order. */
    [[nodiscard]] static MomentWriter for_cqmom(const std::filesystem::path& path, const RawMomentRequest& request, const std::size_t precision = 17);

    /** Append one row for a one-dimensional raw moment sequence. */
    void write(double time, const MomentSet& moments);

    /** Append one row for a multidimensional raw moment tensor. */
    void write(double time, const RawMomentTensor& moments);

    /** Append one row for CQMOM required moment data. */
    void write(double time, const CqmomMomentData& moments);

    /** Append one row for requested moments reconstructed from CQMOM quadrature. */
    void write(double time, const CqmomQuadrature& quadrature);

private:
    enum class MomentLayout
    {
        qmom,
        cqmom_dense,
        cqmom_moment_data,
        cqmom_requested
    };

    MomentWriter(const std::filesystem::path& path, MomentLayout layout, std::vector<std::size_t> extents, std::vector<std::size_t> flat_indices,
                 std::vector<std::string> labels, const std::size_t precision);
    MomentWriter(const std::filesystem::path& path, MomentLayout layout, CqmomMomentLayout cqmom_layout, std::vector<std::string> labels,
                 const std::size_t precision);
    MomentWriter(const std::filesystem::path& path, MomentLayout layout, RawMomentRequest request, std::vector<std::string> labels,
                 const std::size_t precision);

    void write_header();
    void write_row(double time, const std::vector<double>& values);
    void write_aligned_row(double time, const std::vector<double>& values);

    std::ofstream output_;
    MomentLayout layout_ = MomentLayout::qmom;
    std::vector<std::size_t> extents_;
    std::vector<std::size_t> flat_indices_;
    std::optional<CqmomMomentLayout> cqmom_layout_;
    RawMomentRequest raw_moment_request_;
    std::vector<std::string> labels_;
    std::size_t precision;
    std::size_t column_width;
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
    std::size_t precision;
    std::size_t column_width;
};

} // namespace qmomkit
