#ifndef MPM_DDMP_EXPLICIT_H_
#define MPM_DDMP_EXPLICIT_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_base.h"

namespace mpm {

//! DDMPExplicit class
//! \brief A class that implements the fully explicit one phase ddmp
//! \details A single-phase explicit DDMP
//! \tparam Tdim Dimension
template <unsigned Tdim>
class DDMPExplicit : public MPMBase<Tdim> {
public:
  //! Default constructor
  DDMPExplicit(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  //! \param[in] phase Phase to smooth pressure
  void compute_stress_strain(unsigned phase);

protected:
  //! Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Output stress
  using mpm::MPMBase<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMBase<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMBase<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMBase<Tdim>::post_process_;
  //! Logger
  using mpm::MPMBase<Tdim>::console_;
  //! MPM Scheme
  using mpm::MPMBase<Tdim>::mpm_scheme_;
  //! Stress update method
  using mpm::MPMBase<Tdim>::stress_update_;
  //! Interface scheme
  using mpm::MPMBase<Tdim>::contact_;

#ifdef USE_GRAPH_PARTITIONING
  //! Graph
  using mpm::MPMBase<Tdim>::graph_;
#endif

  //! Velocity update
  using mpm::MPMBase<Tdim>::velocity_update_;
  //! Gravity
  using mpm::MPMBase<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMBase<Tdim>::mesh_;
  //! Materials
  using mpm::MPMBase<Tdim>::materials_;
  //! Damping type
  using mpm::MPMBase<Tdim>::damping_type_;
  //! Damping factor
  using mpm::MPMBase<Tdim>::damping_factor_;
  //! Locate particles
  using mpm::MPMBase<Tdim>::locate_particles_;

private:
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  //! Interface
  bool interface_{false};
}; //! DDMPExplicit class
} // namespace mpm

#include "ddmp_explicit.tcc"

#endif // MPM_DDMP_EXPLICIT_H_