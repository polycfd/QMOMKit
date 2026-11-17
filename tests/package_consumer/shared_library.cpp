#include <cstddef>
#include <qmomkit/qmomkit.hpp>

// A shared library of a host that calls into the installed static library. On
// Linux it links only if QMOMKit was compiled as position-independent code.
std::size_t package_consumer_shared_library_node_count()
{
    const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet({1.0, 0.8, 2.8, 4.4}), 2);
    return result.ok() ? result.weights.size() : 0;
}
