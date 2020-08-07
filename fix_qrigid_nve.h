/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qrigid/nve,FixQRigidNVE)

#else

#ifndef LMP_FIX_QRIGID_NVE_H
#define LMP_FIX_QRIGID_NVE_H

#include "fix_qrigid_nh.h"

namespace LAMMPS_NS {

class FixQRigidNVE : public FixQRigidNH {
 public:
  FixQRigidNVE(class LAMMPS *, int, char **);
  ~FixQRigidNVE() {}
};

}

#endif
#endif
