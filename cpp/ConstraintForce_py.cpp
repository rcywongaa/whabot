#include "pybind11/eigen.h"
#include "pybind11/eval.h"
#include "pybind11/pybind11.h"

#include "drake/bindings/pydrake/common/cpp_template_pybind.h"
//#include "drake/bindings/pydrake/common/default_scalars_pybind.h"
//#include "drake/bindings/pydrake/common/deprecation_pybind.h"
//#include "drake/bindings/pydrake/common/drake_optional_pybind.h"
//#include "drake/bindings/pydrake/common/type_pack.h"
//#include "drake/bindings/pydrake/common/type_safe_index_pybind.h"
//#include "drake/bindings/pydrake/documentation_pybind.h"
#include "drake/bindings/pydrake/pydrake_pybind.h"

#include "ConstraintForce.hpp"

PYBIND11_MODULE(tree, m) {
  using namespace drake::multibody;

  m.doc() = "Bindings for custom ConstraintForce";

  py::module::import("pydrake.common.eigen_geometry");
  py::module::import("pydrake.multibody.math");

  {
    using Class = LinearSpringDamper<T>;
    constexpr auto& cls_doc = doc.LinearSpringDamper;
    auto cls = DefineTemplateClassWithDefault<Class, ForceElement<T>>(
        m, "LinearSpringDamper", param, cls_doc.doc);
    cls  // BR
        .def(py::init<const Body<T>&, const Vector3<double>&, const Body<T>&,
                 const Vector3<double>&, double, double, double>(),
            py::arg("bodyA"), py::arg("p_AP"), py::arg("bodyB"),
            py::arg("p_BQ"), py::arg("free_length"), py::arg("stiffness"),
            py::arg("damping"), cls_doc.ctor.doc);
  }
}

