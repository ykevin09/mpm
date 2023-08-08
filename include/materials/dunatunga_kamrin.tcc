//! Constructor with material properties
template <unsigned Tdim>
mpm::DunatungaKamrin<Tdim>::DunatungaKamrin(unsigned int id,
                                            const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // General parameters
    // Density
    density_ = material_properties.at("density").template get<double>();
    // Young's modulus
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    // Critical density
    density_c_ = material_properties.at("density_c").template get<double>();
    // Static friction coefficient at zero I
    mu_s_ = material_properties.at("mu_s").template get<double>();
    // Limiting friction coefficient at high I
    mu_2_ = material_properties.at("mu_2").template get<double>();
    // Material constant
    I_0_ = material_properties.at("I_0").template get<double>();
    // Density of solid grains
    density_s_ = material_properties.at("density_s").template get<double>();
    // Mean particle size
    diameter_ = material_properties.at("diameter").template get<double>();

    // Calculate shear modulus
    shear_modulus_ = youngs_modulus_/(2.*(1.+poisson_ratio_));

    // Properties
    properties_ = material_properties;
    // Set elastic tensor
    this->compute_elastic_tensor();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::DunatungaKamrin<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {
      // Failure state
      {"failure_state", 0},
      // Equivalent plastic shear strain rate
      {"eq_dot_gamma_p", 0.},
      // Mean pressure
      {"pressure", 0.}
  };
  return state_vars;
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::DunatungaKamrin<Tdim>::state_variables() const {
  const std::vector<std::string> state_vars = {
    "failure_state", "eq_dot_gamma_p", "pressure"
  };
  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::DunatungaKamrin<Tdim>::compute_elastic_tensor() {
  // Bulk modulus
  const double G = shear_modulus_;
  const double K = youngs_modulus_/(3.*(1.-2.*poisson_ratio_));
  const double a1 = K + 4./3.*G;
  const double a2 = K - 2./3.*G;
  // compute elastic stiffness matrix
  de_ << a1, a2, a2, 0., 0., 0.,
         a2, a1, a2, 0., 0., 0.,
         a2, a2, a1, 0., 0., 0.,
         0., 0., 0.,  G, 0., 0.,
         0., 0., 0., 0.,  G, 0.,
         0., 0., 0., 0., 0., 0.;

  return true;
}

//! Compute stress invariants
template <unsigned Tdim>
void mpm::DunatungaKamrin<Tdim>::compute_stress_invariants(const Vector6d& stress,
                                                           double* p,
                                                           double* tau_bar) {
  // Compute mean stress p
  *p = mpm::materials::p(stress);
  // Compute tau_bar
  *tau_bar = mpm::materials::q(stress) / sqrt(3.);
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::DunatungaKamrin<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  // Particle density (already updated)
  const double rho = ptr->mass_density();
  // Trial stress - elastic
  Vector6d stress_tr = stress + (this->de_ * dstrain);
  // Get stress invariants
  const double p_tr = -mpm::materials::p(stress_tr);
  const double tau_bar_tr = mpm::materials::q(stress_tr)/sqrt(3.);
  // Update stress
  Vector6d updated_stress;
  updated_stress << 0., 0., 0., 0., 0., 0.;

  double eq_dot_gammap = 0.;
  (*state_vars).at("eq_dot_gamma_p") = eq_dot_gammap;
  (*state_vars).at("pressure") = p_tr;

  if (rho < density_c_ || p_tr <= 0.) {
    (*state_vars).at("failure_state") = 0;
    (*state_vars).at("pressure") = 0.;

    return updated_stress;
  }
  else {
    const double S_0 = mu_s_*p_tr;
    if (tau_bar_tr <= S_0) {
      (*state_vars).at("failure_state") = 1;

      return stress_tr;
    }

    const double S_2 = mu_2_*p_tr;
    const double xi = I_0_/sqrt(diameter_*diameter_*density_s_);
    const double alpha = xi*shear_modulus_*sqrt(p_tr)*ptr->dt();
    const double B = S_2 + tau_bar_tr + alpha;
    const double H = S_2*tau_bar_tr + S_0*alpha;
    const double tmp = B*B - 4.*H;
    if (tmp <= 0.)
      console_->error("DunatungaKamrin: B*B-4.*H={} <= 0.", tmp);
    const double updated_tau_bar = 2.*H/(B+sqrt(tmp));
    eq_dot_gammap = 1./(shear_modulus_*ptr->dt())*(tau_bar_tr-updated_tau_bar);
    if (eq_dot_gammap < 0.)
      console_->error("DunatungaKamrin: eq_dot_gammap={} < 0.", eq_dot_gammap);

    Vector6d I; I << 1., 1., 1., 0., 0., 0.;
    Vector6d s_tr = stress_tr + p_tr*I;
    updated_stress = updated_tau_bar/tau_bar_tr*s_tr - p_tr*I;

    (*state_vars).at("failure_state") = 2;
    (*state_vars).at("eq_dot_gamma_p") = eq_dot_gammap;

    return updated_stress;
  }
}