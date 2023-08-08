#ifndef MPM_MATERIAL_DUNATUNGA_KAMRIN_H_
#define MPM_MATERIAL_DUNATUNGA_KAMRIN_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! DunatungaKamrin class
//! \brief Drucker-Prager + mu(I) rheology (Dunatunga and Kamrin (2015))
//! \details DunatungaKamrin class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class DunatungaKamrin : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  DunatungaKamrin(unsigned id, const Json& material_properties);

  //! Destructor
  ~DunatungaKamrin() override = default;

  //! Delete copy constructor
  DunatungaKamrin(const DunatungaKamrin&) = delete;

  //! Delete assignment operator
  DunatungaKamrin& operator=(const DunatungaKamrin&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! State variables
  std::vector<std::string> state_variables() const override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain increment
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Compute elastic tensor
  bool compute_elastic_tensor();

  //! Compute stress invariants (p, tau_bar)
  //! \param[in] stress Stress
  //! \param[in|out] p Mean stress
  //! \param[in|out] tau_bar Square root of J2=1/2*s*s
  void compute_stress_invariants(const Vector6d& stress, double *p, double *tau_bar);

  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Critical density
  double density_c_{std::numeric_limits<double>::max()};
  //! Static friction coefficient at zero I
  double mu_s_{std::numeric_limits<double>::max()};
  //! Limiting friction coefficient at high I
  double mu_2_{std::numeric_limits<double>::max()};
  //! Material constant
  double I_0_{std::numeric_limits<double>::max()};
  //! Density of solid grains
  double density_s_{std::numeric_limits<double>::max()};
  //! Mean particle size
  double diameter_{std::numeric_limits<double>::max()};
  //! Shear modulus
  double shear_modulus_{std::numeric_limits<double>::max()};
}; // DunatungaKamrin class

} // namespace mpm

#include "dunatunga_kamrin.tcc"

#endif  // MPM_MATERIAL_DUNATUNGA_KAMRIN_H_
