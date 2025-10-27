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

#include "mpi.h"
#include "string.h"
#include "diag_bcc_test.h"
#include "app_bcc_selfdiffusion.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagBccTest::DiagBccTest(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (strcmp(app->style,"bcc/selfdiffusion") != 0)
    error->all(FLERR,"Diag_style bcc/test requires app_style bcc_test");
}

/* ---------------------------------------------------------------------- */

void DiagBccTest::init()
{
  appdiff = (AppBccSelfdiffusion *) app;
  diag_naccept_danni = diag_naccept_danni_all = 12345;
  diag_naccept_vanni = diag_naccept_vanni_all  = 42;
  diag_naccept_dvanni = diag_naccept_dvanni_all  = 22;
  diag_NumD = diag_NumD_all = 0;
  diag_NumV = diag_NumV_all = 0;
}

/* ---------------------------------------------------------------------- */

void DiagBccTest::compute()
{
  // collect first group of tracked events 

  diag_naccept_danni  = appdiff->naccept_danni;
  diag_naccept_vanni  = appdiff->naccept_vanni;
  diag_naccept_dvanni = appdiff->naccept_dvanni;

  MPI_Allreduce(&diag_naccept_danni,&diag_naccept_danni_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&diag_naccept_vanni,&diag_naccept_vanni_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&diag_naccept_dvanni,&diag_naccept_dvanni_all,1,MPI_INT,MPI_SUM,world);

  // collect next group of tracked events 

  diag_naccept_rot   = appdiff->naccept_rot;
  diag_naccept_nntr  = appdiff->naccept_nntr;
  diag_naccept_nnt   = appdiff->naccept_nnt;
  diag_naccept_nnntr = appdiff->naccept_nnntr;
  
  MPI_Allreduce(&diag_naccept_rot,&diag_naccept_rot_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&diag_naccept_nntr,&diag_naccept_nntr_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&diag_naccept_nnt,&diag_naccept_nnt_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&diag_naccept_nnntr,&diag_naccept_nnntr_all,1,MPI_INT,MPI_SUM,world);

  // collect final tracked events 

  diag_naccept_Vnn = appdiff->naccept_Vnn;

  MPI_Allreduce(&diag_naccept_Vnn,&diag_naccept_Vnn_all,1,MPI_INT,MPI_SUM,world);

  // collect number of D and V

  diag_NumD = appdiff->NumD;
  diag_NumV = appdiff->NumV;

  MPI_Allreduce(&diag_NumD,&diag_NumD_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&diag_NumV,&diag_NumV_all,1,MPI_INT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void DiagBccTest::stats(char *strtmp)
{

  // print collected variables across all processors 

  sprintf(strtmp," %7i %7i %7i %7i %8i %7i %7i %7i ",
          diag_NumD_all,diag_NumV_all,
          diag_naccept_danni_all,diag_naccept_vanni_all,diag_naccept_dvanni_all,
	  diag_naccept_rot_all,diag_naccept_nntr_all,diag_naccept_Vnn_all);
}

/* ---------------------------------------------------------------------- */

void DiagBccTest::stats_header(char *strtmp)
{
  sprintf(strtmp," %7s %7s %7s %7s %8s %7s %7s %7s",
	  "NumD","NumV","#d-anni","#v-anni","#dv-anni","#rot","#nntr","#Vnn");
}

