#include <pybind11/pybind11.h>
#include <vil/vil_smart_ptr.h>

// Define vil_smart_pointer for vil_image_resource_sptr
PYBIND11_DECLARE_HOLDER_TYPE(vil_class, vil_smart_ptr<vil_class>, true);

// pybind11 assumes .get() member is available
namespace pybind11 { namespace detail {
    template <typename vil_class>
    struct holder_helper<vil_smart_ptr<vil_class>> { // <-- specialization
        static const vil_class *get(const vil_smart_ptr<vil_class> &p) { return p.ptr(); }
    };
}}
