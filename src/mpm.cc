#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "ddmp_explicit.h"

namespace mpm {
// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_3d("MPMExplicit3D");

// 2D Explicit DDMP
static Register<mpm::MPM, mpm::DDMPExplicit<2>, const std::shared_ptr<mpm::IO>&>
    ddmp_explicit_2d("DDMPExplicit2D");

// 2D Explicit DDMP
static Register<mpm::MPM, mpm::DDMPExplicit<3>, const std::shared_ptr<mpm::IO>&>
    ddmp_explicit_3d("DDMPExplicit3D");
}  // namespace mpm
