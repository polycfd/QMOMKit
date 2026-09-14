#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] inline bool all_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

// A callback that fills a library-supplied buffer rather than returning one may
// write only part of it, and the library would then use whatever it happened to
// leave there. Size and finiteness cannot detect that: the buffer is the right
// size and zero is a perfectly ordinary value.
//
// Debug builds therefore fill the buffer with a recognisable quiet-NaN payload
// before handing it over, so that any entry still holding it is a missing write.
// Release builds compile the marking away and keep the plain zero
// initialisation, so the check costs nothing where it is not wanted. Comparison
// is by bit pattern because NaN never compares equal to itself.
#ifndef NDEBUG
inline constexpr std::uint64_t output_buffer_poison_bits = 0x7ff8000051434b50ULL;

[[nodiscard]] inline double output_buffer_poison() noexcept
{
    double value = 0.0;
    std::memcpy(&value, &output_buffer_poison_bits, sizeof(value));
    return value;
}

[[nodiscard]] inline bool is_output_buffer_poison(const double value) noexcept
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits == output_buffer_poison_bits;
}
#endif

/** Mark every entry as unwritten so a missing write is distinguishable from a
 *  deliberate zero. A no-op in release builds. */
inline void poison_output_buffer([[maybe_unused]] std::vector<double>& values) noexcept
{
#ifndef NDEBUG
    std::fill(values.begin(), values.end(), output_buffer_poison());
#endif
}

} // namespace qmomkit::internal
