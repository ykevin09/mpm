#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "gimp_explicit.h"

namespace mpm {
// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
  mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
  mpm_explicit_3d("MPMExplicit3D");

// 2D Explicit GIMP
static Register<mpm::MPM, mpm::GIMPExplicit<2>, const std::shared_ptr<mpm::IO>&>
  gimp_explicit_2d("GIMPExplicit2D");

// 3D Explicit GIMP
static Register<mpm::MPM, mpm::GIMPExplicit<3>, const std::shared_ptr<mpm::IO>&>
  gimp_explicit_3d("GIMPExplicit3D");

}  // namespace mpm
