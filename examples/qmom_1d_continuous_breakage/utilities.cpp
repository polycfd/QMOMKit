#include "utilities.hpp"

#include <iostream>

namespace qmom_1d_continuous_breakage
{

void print_configuration_errors(const std::vector<qmomkit::config::ConfigurationError>& errors)
{
    for (const qmomkit::config::ConfigurationError& error : errors)
    {
        std::cerr << error.source_name << ':' << error.line << ": " << error.message << '\n';
    }
}

} // namespace qmom_1d_continuous_breakage
