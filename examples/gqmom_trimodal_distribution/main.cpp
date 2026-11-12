#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>

int main()
{
    // The case: the moments of the trimodal distribution, the node count N and
    // the generalized family with its N_G nodes, all from configuration.qmom.
    const std::filesystem::path configuration_file{"configuration.qmom"};
    const auto configured_case = qmomkit::config::load_qmom(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    const auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::QmomCase& gqmom_case = configured_case.value();
    const qmomkit::OutputOptions& output = configured_output.value();

    // GQMOM: N_G nodes that reproduce m0..m(2N).
    const qmomkit::QmomResult gqmom = qmomkit::compute_qmom(gqmom_case);
    if (!gqmom.ok())
    {
        std::cerr << "GQMOM failed: " << gqmom.message << '\n';
        return 1;
    }

    // Classical QMOM of the same moments, for comparison: N nodes that
    // reproduce m0..m(2N-1).
    qmomkit::QmomCase qmom_case = gqmom_case;
    qmom_case.qmom_options.family = qmomkit::QmomFamily::classical;
    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(qmom_case);
    if (!qmom.ok())
    {
        std::cerr << "QMOM failed: " << qmom.message << '\n';
        return 1;
    }

    if (output.write_quadrature)
    {
        std::filesystem::create_directories(output.directory);
        qmomkit::QuadratureWriter::for_qmom(output.directory / "gqmom_quadrature.txt", gqmom.weights.size()).write(0.0, gqmom);
        qmomkit::QuadratureWriter::for_qmom(output.directory / "qmom_quadrature.txt", qmom.weights.size()).write(0.0, qmom);
    }

    std::cout << "m0 = " << gqmom_case.moments[0] << '\n';
    std::cout << "GQMOM: " << gqmom.weights.size() << " nodes on a base of " << gqmom.base_node_count << ", from xi = " << gqmom.abscissas.front()
              << " to " << gqmom.abscissas.back() << '\n';
    std::cout << "QMOM: " << qmom.weights.size() << " nodes\n";
    for (std::size_t node = 0; node < qmom.weights.size(); ++node)
    {
        std::cout << "  xi = " << qmom.abscissas[node] << ", weight = " << qmom.weights[node] << '\n';
    }
    if (output.write_quadrature)
    {
        std::cout << "Quadratures written to " << output.directory.string() << "/gqmom_quadrature.txt and qmom_quadrature.txt\n";
    }
    return 0;
}
