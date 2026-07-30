#include "utilities.hpp"

#include <filesystem>

namespace qmom_1d_explicit_strang_splitting
{

void ensure_output_directory(const std::filesystem::path& output_directory)
{
    std::filesystem::create_directories(output_directory);
}

} // namespace qmom_1d_explicit_strang_splitting
