#include "material.h"
#include "bingham.h"
#include "linear_elastic.h"
#include "modified_cam_clay.h"
#include "mohr_coulomb.h"
#include "newtonian.h"
#include "norsand.h"
#include "dunatunga_kamrin.h"

// Bingham 2D
static Register<mpm::Material<2>, mpm::Bingham<2>, unsigned, const Json&>
    bingham_2d("Bingham2D");

// Bingham 3D
static Register<mpm::Material<3>, mpm::Bingham<3>, unsigned, const Json&>
    bingham_3d("Bingham3D");

// LinearElastic 2D
static Register<mpm::Material<2>, mpm::LinearElastic<2>, unsigned, const Json&>
    linear_elastic_2d("LinearElastic2D");

// LinearElastic 3D
static Register<mpm::Material<3>, mpm::LinearElastic<3>, unsigned, const Json&>
    linear_elastic_3d("LinearElastic3D");

// ModifiedCamClay 2D
static Register<mpm::Material<2>, mpm::ModifiedCamClay<2>, unsigned,
                const Json&>
    modified_cam_clay_2d("ModifiedCamClay2D");

// ModifiedCamClay 3D
static Register<mpm::Material<3>, mpm::ModifiedCamClay<3>, unsigned,
                const Json&>
    modified_cam_clay_3d("ModifiedCamClay3D");

// MohrCoulomb 2D
static Register<mpm::Material<2>, mpm::MohrCoulomb<2>, unsigned, const Json&>
    mohr_coulomb_2d("MohrCoulomb2D");

// MohrCoulomb 3D
static Register<mpm::Material<3>, mpm::MohrCoulomb<3>, unsigned, const Json&>
    mohr_coulomb_3d("MohrCoulomb3D");

// Newtonian 2D
static Register<mpm::Material<2>, mpm::Newtonian<2>, unsigned, const Json&>
    newtonian_2d("Newtonian2D");

// Newtonian 3D
static Register<mpm::Material<3>, mpm::Newtonian<3>, unsigned, const Json&>
    newtonian_3d("Newtonian3D");

// Norsand 2D
static Register<mpm::Material<2>, mpm::NorSand<2>, unsigned, const Json&>
    nor_sand_2d("NorSand2D");

// Norsand 3D
static Register<mpm::Material<3>, mpm::NorSand<3>, unsigned, const Json&>
    nor_sand_3d("NorSand3D");

// DunatungaKamrin 2D
static Register<mpm::Material<2>, mpm::DunatungaKamrin<2>, unsigned, const Json&>
    dunatunga_kamrin_2d("DunatungaKamrin2D");

// DunatungaKamrin 3D
static Register<mpm::Material<3>, mpm::DunatungaKamrin<3>, unsigned, const Json&>
    dunatunga_kamrin_3d("DunatungaKamrin3D");