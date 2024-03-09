//! Constructor
template <unsigned Tdim>
mpm::DDMPExplicit<Tdim>::DDMPExplicit(const std::shared_ptr<IO>& io)
  : mpm::MPMBase<Tdim>(io)
{
  //! Logger
  console_ = spdlog::get("DDMPExplicit");
  //! Stress update
  if (this->stress_update_ == "usl")
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSL<Tdim>>(mesh_, dt_);
  else
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSF<Tdim>>(mesh_, dt_);

  //! Interface scheme
  if (this->interface_)
    contact_ = std::make_shared<mpm::ContactFriction<Tdim>>(mesh_);
  else
    contact_ = std::make_shared<mpm::Contact<Tdim>>(mesh_);
}

//! DDMP Explicit compute stress strain
template <unsigned Tdim>
void mpm::DDMPExplicit<Tdim>::compute_stress_strain(unsigned phase) {

}

//! DDMP Explicit solver
template <unsigned Tdim>
bool mpm::DDMPExplicit<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialize MPI rank and size
  int mpi_rank = 0, mpi_size = 0;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Phase
  const unsigned phase = 0;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Enable repartitioning if resume is done with particles generated outside
  // the MPM code.
  bool repartition = false;
  if (analysis_.find("resume") != analysis_.end() &&
      analysis_["resume"].find("repartition") != analysis_["resume"].end())
    repartition = analysis_["resume"]["repartition"].template get<bool>();

  // Pressure smoothing
  pressure_smoothing_ = io_->analysis_bool("pressure_smoothing");

  // Initialise material
  this->initialise_materials();

  // Initialise mesh
  this->initialise_mesh();
  mesh_->find_node_neighbours();

  return status;
}


