#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit_twophase.h"
#include "write_mesh_particles_unitcell.h"

// Check MPM Explicit USF
TEST_CASE("MPM 2D Explicit TwoPhase USF implementation is checked in unitcells",
          "[MPM][2D][USF][Explicit][2Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit-twophase-usf";
  const std::string analysis = "MPMExplicitTwoPhase2D";
  const std::string stress_update = "usf";
  REQUIRE(mpm_test::write_json_unitcell_twophase(2, analysis, stress_update,
                                                 fname) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_2d_unitcell() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_2d_unitcell() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-twophase-usf-2d-unitcell.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE(mpm->initialise_materials() == true);

    // Initialise mesh and particles
    REQUIRE(mpm->initialise_mesh() == true);
    REQUIRE(mpm->initialise_particles() == true);

    // Initialise external loading
    REQUIRE(mpm->initialise_loads() == true);

    // Renitialise materials
    REQUIRE(mpm->initialise_materials() == false);
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

// Check MPM Explicit
TEST_CASE("MPM 3D Explicit TwoPhase USF implementation is checked in unitcells",
          "[MPM][3D][Explicit][USF][2Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit-twophase-usf";
  const std::string analysis = "MPMExplicitTwoPhase3D";
  const std::string stress_update = "usf";
  REQUIRE(mpm_test::write_json_unitcell_twophase(3, analysis, stress_update,
                                                 fname) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_3d_unitcell() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_3d_unitcell() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-twophase-usf-3d-unitcell.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE(mpm->initialise_materials() == true);

    // Initialise mesh and particles
    REQUIRE(mpm->initialise_mesh() == true);
    REQUIRE(mpm->initialise_particles() == true);

    // Renitialise materials
    REQUIRE(mpm->initialise_materials() == false);
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}