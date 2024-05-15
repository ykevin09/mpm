//! Constructor
template <unsigned Tdim>
mpm::GIMPExplicit<Tdim>::GIMPExplicit(const std::shared_ptr<IO>& io)
  : mpm::MPMBase<Tdim>(io)
{
  //! Logger
  console_ = spdlog::get("GIMPExplicit");
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

//! GIMP Explicit compute stress strain
template <unsigned Tdim>
void mpm::GIMPExplicit<Tdim>::compute_stress_strain(unsigned int phase) {

}

//! GIMP Explicit solver
template <unsigned Tdim>
bool mpm::GIMPExplicit<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialize MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

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

  // Interface
  interface_ = io_->analysis_bool("interface");

  // Initialize material
  this->initialise_materials();

  // Initialize mesh
  this->initialise_mesh();

  // Initialize particles
  if (!resume) this->initialise_particles();

  // Create nodal properties
  if (interface_) mesh_->create_nodal_properties();

  // Compute mass
  if (!resume)
    mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  bool initial_step = !resume;
  // Check point resume
  if (resume) {
    this->checkpoint_resume();
    if (repartition) {
      this->mpi_domain_decompose(initial_step);
    } else {
      mesh_->resume_domain_cell_ranks();
#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
  } else {
    // Domain decompose
    this->mpi_domain_decompose(initial_step);
  }

  //! Particle entity sets and velocity constraints
  if (resume) {
    this->particle_entity_sets(false);
    this->particle_velocity_constraints();
  }

  // Initialise loading conditions
  this->initialise_loads();

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {

    if (mpi_rank == 0 && step_ % output_steps_ == 0) {
      console_->info("Step: {} of {}.", step_, nsteps_);
    }

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    // Run load balancer at a specified frequency
    if (step_ % nload_balance_steps_ == 0 && step_ != 0)
      this->mpi_domain_decompose(false);
#endif
#endif

    // Inject particles: particle nodes
    mesh_->inject_particles(step_*dt_);

    // Initialize nodes, cells, and shape functions
    mesh_->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

    auto cells = mesh_->cells();
    for (auto citr=cells.cbegin(); citr != cells.cend(); ++citr) {
      if ((*citr)->nparticles() > 0) {
        for (auto node: (*citr)->nodes()) {
          node->assign_status(true);
        }

        for (auto neighbour: (*citr)->neighbours()) {
          auto cell = mesh_->cell(neighbour);
          for (auto node: cell->nodes()) {
            node->assign_status(true);
          }
        }
      }
    }

    mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));


    break;
  }

  return status;
}