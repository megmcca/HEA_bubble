/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(bcc_test,DiagBccTest)

#else

#ifndef SPK_DIAG_BCC_TEST_H
#define SPK_DIAG_BCC_TEST_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagBccTest : public Diag {

 public:
  DiagBccTest(class SPPARKS *, int, char **);
  ~DiagBccTest() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  // class AppBccTest *appdiff;
  class AppBccSelfdiffusion *appdiff;
  int diag_naccept_danni;
  int diag_naccept_vanni;
  int diag_naccept_dvanni;
  int diag_naccept_rot;
  int diag_naccept_nntr;
  int diag_naccept_nnt;
  int diag_naccept_nnntr;
  int diag_naccept_VDex;
  int diag_naccept_Vnn;
  int diag_NumD;
  int diag_NumV;

  int diag_naccept_danni_all;
  int diag_naccept_vanni_all;
  int diag_naccept_dvanni_all;
  int diag_naccept_rot_all;
  int diag_naccept_nntr_all;
  int diag_naccept_nnt_all;
  int diag_naccept_nnntr_all;
  int diag_naccept_VDex_all;
  int diag_naccept_Vnn_all;
  int diag_NumD_all;
  int diag_NumV_all;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag_style bcc_test requires app_style bcc_test

Self-explanatory.

*/
