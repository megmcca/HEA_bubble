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

#ifdef APP_CLASS
AppStyle(bcc/selfdiffusion,AppBccSelfdiffusion)

#else

#ifndef SPK_APP_BCC_SELFDIFFUSION_H
#define SPK_APP_BCC_SELFDIFFUSION_H

#include "app_lattice.h"
#include <map>
#include <set>
#include <string>
#include <utility>

namespace SPPARKS_NS {

class AppBccSelfdiffusion : public AppLattice {

 public:
  AppBccSelfdiffusion(class SPPARKS *, int, char **);
  ~AppBccSelfdiffusion();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

  int naccept_danni;
  int naccept_vanni;
  int naccept_dvanni;
  int naccept_rot;
  int naccept_nntr;
  int naccept_nnt;
  int naccept_nnntr;
  int naccept_VFex;

 private:
  int engstyle;
  int allocated;
  int *esites;
  int *echeck;

  int dimension;
  int *sitetype;
  int *lattice;
  double total_R;

  struct Event {           // one event for an owned site
    double propensity;     // propensity of this event
    int destination;       // local ID of destination site
    int style;		   // annihilation, rotation,...
			   // ...translation+rotation,exchange
    int next;              // index of next event for this site
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void parse_bccselfdiffusion(int narg, char **arg);

  // These variables hold the rates for various events.
  // They are initialized with default values.
  double V_anniAtSinks = 0.0001; // Vacancy Annihilation at Sinks
  double D_anniAtSinks = 0.001;  // Dumbbell Annihilation at Sinks
  double VD_mutualAnni = 0.5;	 // Vacancy-Dumbbell Annihhilation rate
  double Rrot = 0.01;            // Onsite Rotation Rate
  double RTrans=(1.0-RTRansRot); // for nn=1 to 8, Translation Rate. 
  double RTransRot = 0;          // RTransRot (is ~dumbbell diff coeff).
  double RnnnTrans;		 // Translation + rot to nnn.
  double V_DiffCoeff = 0.0001;	 // Vacancy Diffusion Coeff 
				 // as fraction of Dumbbell Diffusion Coeff.

  double P_rot;
  double P_nntr;
  double P_nnt = 0.0;
  double P_nnntr = 0.0;
  double P_vdiff;

  // These variables are bond energies between DD, VV and DV.
  // They are initialized with default values.
  double VV = 1.0;	//Vacancy-Vacancy neighbor
  double DD = 2.0;	//Dumbbell-Dumbbell neighbor
  double DV = -1.0;	//Dumbbell-Vacancy neighbor

  void site_event_linear(int, class RandomPark *);

  void clear_events(int);

  void add_event(int, int, double, int);

  void allocate_data();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot use %s command until sites exist

This command requires sites exist before using it in an input script.

E: Can only use ecoord command with app_style diffusion nonlinear

Self-explanatory.

E: Cannot define Schwoebel barrier without Schwoebel model

Self-explanatory.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: Cannot perform deposition in parallel

UNDOCUMENTED

E: Cannot perform deposition with multiple sectors

UNDOCUMENTED

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Did not reach event propensity threshhold

UNDOCUMENTED

E: BAD DONE

UNDOCUMENTED

*/
