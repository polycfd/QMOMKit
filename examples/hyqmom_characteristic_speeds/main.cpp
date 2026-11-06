#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <string>
#include <vector>

namespace
{

// One moment set of Fox & Laurent (2022), Fig. 1: m0 = 1, mean 0 and variance
// 1, so that the raw moments are the standardized moments S_k, with
// S_3, ..., S_{2n-1} fixed and the Hankel determinant H_{2n} from 0 to its
// largest value in the figure.
struct FigureMomentSet
{
    std::size_t node_count;
    std::vector<double> standardized_moments;
    double largest_hankel_determinant;
};

void print_speeds(const std::size_t node_count, const double hankel_determinant, const std::size_t active_node_count,
                  const std::vector<double>& speeds)
{
    std::cout << "  H_" << 2 * node_count << " = " << hankel_determinant << ", " << active_node_count << " nodes:";
    for (const double speed : speeds)
    {
        std::cout << ' ' << speed;
    }
    std::cout << '\n';
}

} // namespace

int main()
{
    const std::vector<FigureMomentSet> figure{
        {2, {-1.0}, 20.0},
        {3, {-1.0, 5.0, -8.0}, 700.0},
        {4, {-1.0, 5.0, -8.0, 67.3, -100.0}, 1.0e6},
        {5, {-1.0, 5.0, -8.0, 67.3, -100.0, 3000.0, 0.0}, 1.0e11},
    };
    const std::size_t point_count = 101;
    const std::filesystem::path output_directory = "output";
    std::filesystem::create_directories(output_directory);
    std::cout << std::setprecision(6);

    for (const FigureMomentSet& set : figure)
    {
        const std::size_t node_count = set.node_count;
        std::vector<double> moments{1.0, 0.0, 1.0};
        moments.insert(moments.end(), set.standardized_moments.begin(), set.standardized_moments.end());

        // On the boundary of the realizable set, H_{2n} = 0, m_{2n} is that of
        // the classical quadrature of n nodes of m_0, ..., m_{2n-1}.
        qmomkit::QmomOptions fixed_options;
        fixed_options.adaptive = false;
        const qmomkit::QmomResult classical_result = qmomkit::compute_qmom(qmomkit::MomentSet(moments), node_count, fixed_options);
        if (!classical_result.ok())
        {
            std::cerr << "QMOM failed: " << classical_result.message << '\n';
            return 1;
        }
        qmomkit::Quadrature1D classical_quadrature;
        classical_quadrature.abscissas = classical_result.abscissas;
        classical_quadrature.weights = classical_result.weights;
        const double boundary_moment = qmomkit::reconstruct_raw_moments(classical_quadrature, 2 * node_count + 1)[2 * node_count];

        // With m0 = 1, H_{2n-2} is the product over k = 1, ..., n - 1 of
        // b_1 ... b_k, and H_{2n} / H_{2n-2} = b_1 ... b_n = m_{2n} - boundary_moment.
        const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(qmomkit::MomentSet(moments), node_count);
        if (!recurrence.ok())
        {
            std::cerr << "Wheeler recurrence failed: " << recurrence.message << '\n';
            return 1;
        }
        double previous_hankel_determinant = 1.0;
        double recurrence_product = 1.0;
        for (std::size_t k = 1; k < node_count; ++k)
        {
            recurrence_product *= recurrence.b[k];
            previous_hankel_determinant *= recurrence_product;
        }

        // The case sets the family, node count and options of the states made
        // from it; state_candidate(...) does not use the moments of the case.
        moments.push_back(boundary_moment);
        qmomkit::QmomOptions hyperbolic_options;
        hyperbolic_options.family = qmomkit::QmomFamily::hyperbolic;
        const qmomkit::QmomCase qmom_case{qmomkit::MomentSet(moments), node_count, hyperbolic_options};

        const std::filesystem::path file_path = output_directory / ("speeds_n" + std::to_string(node_count) + ".txt");
        std::ofstream file(file_path);
        file.imbue(std::locale::classic());
        file << std::setprecision(17);
        file << "# H_" << 2 * node_count << " active_nodes";
        for (std::size_t speed = 0; speed < 2 * node_count + 1; ++speed)
        {
            file << " speed_" << speed;
        }
        file << '\n';

        std::cout << "n = " << node_count << ": m_" << 2 * node_count << " = " << boundary_moment << " on the boundary, H_" << 2 * node_count - 2
                  << " = " << previous_hankel_determinant << '\n';
        for (std::size_t point = 0; point < point_count; ++point)
        {
            const double hankel_determinant = set.largest_hankel_determinant * static_cast<double>(point) / static_cast<double>(point_count - 1);
            moments.back() = boundary_moment + hankel_determinant / previous_hankel_determinant;

            // The state of these moments, as a host creates it for moments it
            // computed itself, and the characteristic speeds of its moments.
            const qmomkit::QmomStateCandidate candidate = qmom_case.state_candidate(qmomkit::MomentSet(moments));
            if (!candidate.ok())
            {
                std::cerr << "QMOM state failed: " << candidate.message << '\n';
                return 1;
            }
            const qmomkit::HyqmomCharacteristicSpeeds result = qmomkit::compute_hyqmom_characteristic_speeds(*candidate.state);
            if (!result.ok())
            {
                std::cerr << "Characteristic speeds failed: " << result.message << '\n';
                return 1;
            }

            file << hankel_determinant << ' ' << candidate.state->active_node_count();
            for (const double speed : result.speeds)
            {
                file << ' ' << speed;
            }
            file << '\n';
            if (point == 0 || point == point_count - 1)
            {
                print_speeds(node_count, hankel_determinant, candidate.state->active_node_count(), result.speeds);
            }
        }
    }

    std::cout << "Output written to " << output_directory.string() << "/\n";
    return 0;
}
