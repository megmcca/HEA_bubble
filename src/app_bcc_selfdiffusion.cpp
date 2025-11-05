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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "solve.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "app_bcc_selfdiffusion.h"

using namespace SPPARKS_NS;
using std::map;
using std::set;

enum{LINEAR};
enum{ANNIATSINK,ANNI,ONSITEROT,NNTRANSROT,EXCHANGE};

#define DELTAEVENT 100000

// This app is based on the bcc_selfdiffusion app 
// These are the significant changes 
// 1. Six components: empty interstial sites, vacancies, Metal (M) atoms, 
//    M-M dumbbells, He atoms on regular sites (Hr), He atoms on 
//    interstial sites.
// 2. Vacancies can exchange with He and M and they can be annihilated
//    by Hi hopping into them.
// 3. Hi can hop into empty interstitial sites or into vacancies (#2)
// 4. Dumbbells can exchange with M. 
// 5. Vacancies and M-M dumbbells can be annihilated

/* ---------------------------------------------------------------------- */

AppBccSelfdiffusion::AppBccSelfdiffusion(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  // need to double check these values

  ninteger = 1;
  ndouble = 2;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 0;
  numrandom = 1;
  clustersizecol = 0;

  // no args for this app

  if (narg > 1) error->all(FLERR,"Illegal app_style command");

  engstyle = LINEAR;

  create_arrays();
  esites = NULL;
  echeck = NULL;
  maxevent = 0;
  events = NULL;
  firstevent = NULL;

  allocated = 0;

  naccept_danni = naccept_vanni = naccept_dvanni = 0;
  naccept_rot = naccept_nntr = naccept_nnt = naccept_nnntr = 0;
  naccept_Vnn = 0;
  NumD = NumV = NumHr = NumHi = 0;

}

/* ---------------------------------------------------------------------- */

AppBccSelfdiffusion::~AppBccSelfdiffusion()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->destroy(firstevent);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::input_app(char *command, int narg, char **arg)
{
  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(FLERR,str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

  if (strcmp(command,"bcc/selfdiffusion") == 0)
    parse_bccselfdiffusion(narg,arg);
  else 
   error->all(FLERR,"Unrecognized command in parse");
}

/* ---------------------------------------------------------------------- */

void AppBccSelfdiffusion::parse_bccselfdiffusion(int narg, char **arg)
{
   // 2 args: bcc/selfdiffusion rate <double value>
   // 2 args: bcc/selfdiffusion bondenergy <double value> 
   // 3 args: bcc/selfdiffusion actenergy <double value> Do <double value>

   if (narg < 2)
     error->all(FLERR,"Illegal bcc/selfdiffusion command");

   double THz = 1e12;;
   if (strcmp(arg[0], "rate") == 0) {
       if (narg != 3) {
            error->all(FLERR, "Illegal bcc/selfdiffusion rate command: "
                               "Expected 'rate <name> <value>'");
        } 
        std::string rate_name = arg[1]; // The name of the rate variable 
					// (e.g., "V_anniAtSinks")
        double rate_value = std::atof(arg[2]); // Using std::atof as requested

        // Assign the value to the corresponding member variable
        if (rate_name == "V_anniAtSinks") {
            V_anniAtSinks = rate_value;
        } else if (rate_name == "D_anniAtSinks") {
            D_anniAtSinks = rate_value;
        } else if (rate_name == "VD_mutualAnni") {
            VD_mutualAnni = rate_value;
        } else if (rate_name == "Freq_DRot") {
            Rrot = rate_value * THz;
        } else if (rate_name == "Freq_DDiff") {
            RTransRot = rate_value * THz;
        } else if (rate_name == "Freq_VDiff") {
            RVDiff = rate_value * THz;
        } else if (rate_name == "Freq_HrDiff") {
            RHrDiff = rate_value * THz;
        } else if (rate_name == "Freq_HiDiff") {
            RHiDiff = rate_value * THz;
        } else {
            error->all(FLERR, ("Unknown bcc/selfdiffusion rate name: " + rate_name).c_str());
        }
    } 
    else if (strcmp(arg[0], "actenergy") == 0) {
       if (narg != 4) {
            error->all(FLERR, "Illegal bcc/selfdiffusion actenergy command: "
                               "Expected 'actenergy <name> <value> <value>'");
        } 
        std::string actenergy = arg[1]; // The name of the energy variable (e.g., "QRot")
        double actenergy_value = std::atof(arg[2]); // Using std::atof as requested
        double G_value = std::atof(arg[3]);

        // Assign the value to the corresponding member variable
        if (actenergy == "QDRot") {
            Q_DRot = actenergy_value;
            G_DRot = G_value * THz;
        } else if (actenergy == "QDDiff") {
            Q_DDiff = actenergy_value;
            G_DDiff = G_value * THz;
        } else if (actenergy == "QVDiff") {
            Q_VDiff = actenergy_value;
            G_VDiff = G_value * THz;
        } else if (actenergy == "QHrDiff") {
            Q_HrDiff = actenergy_value;
            G_HrDiff = G_value * THz;
        } else if (actenergy == "QHiDiff") {
            Q_HiDiff = actenergy_value;
            G_HiDiff = G_value * THz;
        } else {
            error->all(FLERR, ("Unknown bcc/selfdiffusion actenergy name: "+std::to_string(actenergy_value)).c_str());
        }
    }
    else if (strcmp(arg[0], "bondenergy") == 0) {
       if (narg != 3) {
            error->all(FLERR, "Illegal bcc/selfdiffusion bondenergy command: "
                               "Expected 'bondenergy <name> <value>'");
        } 
        std::string energy = arg[1]; // The name of the energy variable (e.g., "V-V")
        double bondenergy_value = std::atof(arg[2]); // Using std::atof as requested

        // Assign the value to the corresponding member variable
        if (energy == "VV") {
            VV = bondenergy_value;
        } else if (energy == "DD") {
            DD = bondenergy_value;
        } else if (energy == "DV") {
            DV = bondenergy_value;
        } else if (energy == "VHi") {
            VHi = bondenergy_value;
        } else if (energy == "VHr") {
            VHr = bondenergy_value;
        } else if (energy == "HrHr") {
            HrHr = bondenergy_value;
        } else if (energy == "HiHi") {
            HiHi = bondenergy_value;
        } else {
            error->all(FLERR, ("Unknown bcc/selfdiffusion energy name: "+std::to_string(bondenergy_value)).c_str());
        }
    }
    else {
        error->all(FLERR, "Illegal bcc/selfdiffusion command: "
                           "Expected 'bondenergy', 'actenergy' or 'rate'");
    }
} 


/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::grow_app()
{
  lattice = iarray[0];  //empty=1, vac=2, dumbbell=3-8, M=9-18, Hr=19, Hi=20: i1
  bsize = darray[0];	//size of cluster that site is in
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::init_app()
{
   if (!allocated) allocate_data();
   allocated = 1;

   dimension = domain->dimension;
   dt_sweep = 1.0;

   if(G_DRot == 0 && G_DDiff == 0 && G_VDiff == 0 && Rrot == 0 && RTransRot == 0 && RVDiff == 0)
        error->all(FLERR, "Enter activation energies and attempt frequencies or "
	                  "frequencies, both cannot be 0");

   if((G_DRot != 0.0 || G_DDiff != 0.0 || G_VDiff != 0.0) && (Rrot != 0 || RTransRot != 0 || RVDiff != 0))
        error->all(FLERR, "Either enter activation energies and attempt frequencies "
			  "or frequencies, but not both");

   long double P_hi = 0.0;
   long double sweep = 0.0;

   // Count the number of dumbbells and vacancies to adjust the jump frequency accordingly
   int site;
   for (site = 0; site < nlocal; site++){
       if (lattice[site] == 2) NumV++;
       else if (lattice[site] > 2 && lattice[site] < 9) NumD++;
       else if (lattice[site] == 19) NumHr++;
       else if (lattice[site] == 20) NumHi++;
   }
   // For verification
   //printf("NumD = %i, NumV = %i\n", NumD, NumV);

   if(Rrot || RTransRot || RVDiff){
     P_rot = Rrot;
     P_nntr = RTransRot;
     P_vdiff = RVDiff;
     P_Hrdiff = RHrDiff;
     P_Hidiff = RHiDiff;
   }
   else {
     if(temperature <= 0.0)
        error->all(FLERR, "When using activation energies, "
		          "temperature must be > 0");

     P_rot =       G_DRot  * exp(-Q_DRot/temperature);
     P_nntr =      G_DDiff * exp(-Q_DDiff/temperature);
     P_vdiff =     G_VDiff * exp(-Q_VDiff/temperature);
     if (logfile) { 
            double THz = 1e12; //Since P's are not normalized, convert to THz
            fprintf(logfile, "Event frequencies in THz are:\n");
            fprintf(logfile, "    EF_rot = %f, EF_nntr = %f and EF_vdiff = %f\n", P_rot/THz,P_nntr/THz,P_vdiff/THz);
     }
   }
    // If there are no dumbbells and/or no vacancies, then set corresponding
    // jump frequency, P_* =0.0. 
     if (!NumD) {
        P_rot = 0.0;
        P_nntr = 0.0;
     }
     if (!NumV) {
        P_vdiff = 0.0;
     }
     if (!NumHr) {
        P_Hrdiff = 0.0;
     }
     if (!NumHi) {
        P_Hidiff = 0.0;
     }

     // Rescale rates to make more efficient by setting the highest to have
     // rate  = 1.0
     P_hi = P_nntr;
     if (P_rot > P_hi)
	P_hi = P_rot;
     if (P_vdiff > P_hi)
	P_hi = P_vdiff;
     if (P_Hrdiff > P_hi);
	P_hi = P_Hrdiff;
     if (P_Hidiff > P_hi);
	P_hi = P_Hidiff;
  
     sweep = 1.0/P_hi;

     // Rescale jump frequencies to make more efficient by setting the highest to have
     // jump frequency = 1.0
     P_rot = P_rot/P_hi;
     P_nntr = P_nntr/P_hi;
     P_vdiff = P_vdiff/P_hi;
     P_Hrdiff = P_Hrdiff/P_hi;
     P_Hidiff = P_Hidiff/P_hi;

   if (logfile) { 
            fprintf(logfile, "To obtain the event frequency in THz, multiply the individual P_*'s by the ");
            fprintf(logfile, "Rate of attempts = %.6Le THz\n", P_hi/1e12);
            fprintf(logfile, "P_rot = %f, P_nntr = %f P_vdiff = %f and P_hi = %.6Le\n", P_rot,P_nntr,P_vdiff);
            fprintf(logfile, "1 sweep = %.6Le sec\n\n", sweep);
   }
}

/* ----------------------------------------------------------------------
   setup before each run

------------------------------------------------------------------------- */

void AppBccSelfdiffusion::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppBccSelfdiffusion::site_energy(int i)
{
  // energy of site = sum of bond weights
  // default values are JDD = 2, JVV = 1, JVD = -1,   


   double energy = 0.0;
   int i_V, i_D, i_M, i_Hr, i_Hi;
   int j_V, j_D, j_M, j_Hr, j_Hi;

   i_V = i_D = i_M = i_Hr = i_Hi = 0;
   j_V = j_D = j_M = j_Hr = j_Hi = 0;
  
   int ip = lattice[i];
   if (ip == 2) i_V = 1;
   else if (ip > 2 && ip < 9) i_D = 1;
   else if (ip > 8 && ip < 19) i_M = 1;
   else if (ip == 19) i_Hr = 1;
   else if (ip == 20) i_Hi = 1;

   for (int j = 0; j < maxneigh; j++){
      int nj = neighbor[i][j];
      int jp = lattice[nj];
      if (jp == 2) j_V = 1;
      else if (jp > 2 && jp < 9) j_D = 1;
      else if (jp > 8 && jp < 19) j_M = 1;
      else if (jp == 19) j_Hr = 1;
      else if (jp == 20) j_Hi = 1;

      if (i_D && j_D) //if neighboring dumbbells
          energy += DD;
      else if (i_V && j_V) // if neighboring vacancies
          energy += VV;
      else if ((i_D && j_V) || (i_V && j_D)) //if dumbbell and vac
          energy += DV;
      }
   return 0.5*energy;
}

/* ----------------------------------------------------------------------
   rKMC Method
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::site_event_rejection(int i, RandomPark *random)
{

  // This app assigns spin = 0 for empty interstitials, 1 for empty 
  // interstitial, 2 for vacancy, 3 to 8 for dumbbells in 6 <110> 
  // orientations, 9-18 for M, 19 for Hr and 20 for Hi. 

  // It attempts diffusion of vacancies (vac) by exchanging with metal (M) 
  // and with Hr, He atoms on regular sites. It also annihilates both 
  // vacancy and dumbbell when they meet with rate of VD_mutualAnni.

  // Interstitial He (Hi) diffuse by hoping into neighboring empty sites

  // He on regular sites (Hr) can diffuse by exchanging with M and Vac.

  // Dumbbells can also diffusion and change configuration by rotating  
  // onsite, or by translating to nn M sites & rotating.

  // Both dumbbell interstitials (D) and vacancies can be annihilated at sinks
  // at specified rates.  Note, no mass conservation is enforced, they
  // both simply become regular M sites.

  double einitial,edelta;
  int i_old, j_old, i_new, j_new;
  int i_V, i_D, i_M, i_Hr, i_Hi;
  int j_V, j_D, j_M, j_Hr, j_Hi;
  int j, nbor;
  double P, P_event;

  i_old = lattice[i];

  // if site is empty, return
  if (i_old == 1) return;

  i_V = i_D = i_M = i_Hr = i_Hi = 0;
  j_V = j_D = j_M = j_Hr = j_Hi = 0;
  
  if (i_old == 2) i_V = 1;
  else if (i_old > 2 && i_old < 9) i_D = 1;
  else if (i_old > 8 && i_old < 19) i_M = 1;
  else if (i_old == 19) i_Hr = 1;
  else if (i_old == 20) i_Hi = 1;

  // if site is a vacancy
  if (i_V) {
    // annihilation of vacancies at sinks
    P = random->uniform();
    if (P <= V_anniAtSinks) { 
	lattice[i] = 9;
	naccept++;
	naccept_vanni++;
        NumV--;
	return;
    }

    int nbor = (int) (maxneigh*random->uniform());
    if (nbor >= maxneigh) nbor = maxneigh-1;
    int j = neighbor[i][nbor];
    j_old = lattice[j];
      
    // if another vacancy, return
    if (i_old == j_old) 
	return;
    
    if (j_old > 2 && j_old < 9) j_D = 1;
    else if (j_old > 8 && j_old < 19) j_M = 1;
    else if (j_old == 19) j_Hr = 1;
    else if (j_old == 20) j_Hi = 1;

    // if neigh of V is D, annihilate or return.
    else if (j_D) {
       P = random->uniform();
       if(P <= VD_mutualAnni){
           lattice[i] = 9;
           lattice[j] = 9;
           naccept++;
	   naccept_dvanni++;
           NumD--;
           NumV--;
           return;
       }
       return; 
    }
    // DEBUG
    // if(lattice[j] != 2) printf("Error in vacancy loop\n");

    // else if neighbor is M, attempt exchange
    else if (j_M){
       // Vacancy diffusion rate is P_vdiff 
       P_event = random->uniform();
       if (P_event <= P_vdiff) {
          i_old = lattice[i];
          einitial = site_energy(i)+site_energy(j);

          lattice[i] = j_old;
	  lattice[j] = i_old;
          edelta = site_energy(i)+site_energy(j) - einitial;

          if(edelta <= 0.0){
	    naccept++;
	    naccept_Vnn++;
	    return;
          } 
       else if (temperature > 0.0){ 
          P = random->uniform();
	  if(P <= exp(-1*edelta*t_inverse)) {
	    naccept++;
	    naccept_Vnn++;
	    return;
          }
       }
       else{
	  lattice[i] = i_old; 
	  lattice[j] = j_old;
	  return;
       } 
     }
   }
 } 
  // if i is D, then attempt anni, rot, rot+trans
  else if (i_D){

    //Attempt annihilation of dumbbells at sinks
    P = random->uniform();
    if (P <= D_anniAtSinks) { 
	lattice[i] = 9;
	naccept++;
	naccept_danni++;
        NumD--;
	return;
    }
    
    // Attempt to rotate dumbbells to another randomly selected 
    // orientation that is 60 degrees
    P_event = random->uniform();
    if (P_event <= P_rot) {
	while(!0){
	  i_new = (int) (6*random->uniform() + 3);
	  if (i_new != i_old) {
	   if (i_old == 3){
	     if (i_new != 4) break;
	   } else if (i_old == 4) {
	     if (i_new != 3) break;
	   } else if (i_old == 5) {
	     if (i_new != 6) break;
	   } else if (i_old == 6) {
	     if (i_new != 5) break;
	   } else if (i_old == 7) {
	     if (i_new != 8) break;
	   } else if (i_old == 8) {
	     if (i_new != 7) break;
	     }
           }
	}
	  //DEBUG
          //printf("Onsite rotation, old_spin[%i] = %i, ",i,lattice[i]);
	  lattice[i] = i_new;
          naccept++;
	  naccept_rot++;
    }

    // Attempt translation or translation + rotation with randomly 
    // selected neighbor
    P_event = random->uniform();
    if (P_event <= P_nntr){
      nbor = (int) (4*random->uniform());
      if (nbor >= 4) nbor = 3;

      // Find the randomly selected neighbor and randomly select
      // one of the two possible orientations it can have.
      if(i_old==3){ // [1-10] can rotate to [01-1]=6 & [10-1]=8
        if(nbor == 0) { j = neighbor[i][2];
        } else if(nbor == 1) { j = neighbor[i][3];
        } else if(nbor == 2) { j = neighbor[i][4];
        } else if(nbor == 3) { j = neighbor[i][5];
          }
        j_new = 6;
        if (random->uniform() <= 0.5) j_new = 8;
      } else if (i_old==4) {
        if(nbor == 0) { j = neighbor[i][0];
        } else if(nbor == 1) { j = neighbor[i][1];
        } else if(nbor == 2) { j = neighbor[i][6];
        } else if(nbor == 3) { j = neighbor[i][7];
          }
        j_new = 5;
        if (random->uniform() <= 0.5) j_new = 7;
      } else if (i_old==5) {
        if(nbor == 0) { j = neighbor[i][0];
        } else if(nbor == 1) { j = neighbor[i][3];
        } else if(nbor == 2) { j = neighbor[i][4];
        } else if(nbor == 3) { j = neighbor[i][7];
          }
        j_new = 4;
        if (random->uniform() <= 0.5) j_new = 7;
      } else if (i_old==6) {
        if(nbor == 0) { j = neighbor[i][1];
        } else if(nbor == 1) { j = neighbor[i][2];
        } else if(nbor == 2) { j = neighbor[i][5];
        } else if(nbor == 3) { j = neighbor[i][6];
          }
        j_new = 3;
        if (random->uniform() <= 0.5) j_new = 8;
      } else if (i_old==7) {
        if(nbor == 0) { j = neighbor[i][0];
        } else if(nbor == 1) { j = neighbor[i][2];
        } else if(nbor == 2) { j = neighbor[i][5];
        } else if(nbor == 3) { j = neighbor[i][7];
          }
        j_new = 4;
        if (random->uniform() <= 0.5) j_new = 5;
      } else if (i_old==8) {
        if(nbor == 0) { j = neighbor[i][1];
        } else if(nbor == 1) { j = neighbor[i][3];
        } else if(nbor == 2) { j = neighbor[i][4];
        } else if(nbor == 3) { j = neighbor[i][6];
          }
        j_new = 3;
        if (random->uniform() <= 0.5) j_new = 6;
      }

      j_old = lattice[j];
      if (j_old == 2) j_V = 1;
      else if (j_old > 2 && j_old < 9) j_D = 1;
      else if (j_old > 8 && j_old < 19) j_M = 1;
      else if (j_old == 19) j_Hr = 1;

      //if neighbor is another dumbbell, return
      if (j_D)
	 return;

      // If neighbor is a vacancy, i&j can become regular M sites by
      // mutual annihilation.  If not annihilated, return.
      if (j_V){
        P = random->uniform();
	if(P <= VD_mutualAnni) {
           //DEBUG
           //printf("annihilate, old_spin_i[%i] = %i, old_spin_j[%i] = %i, ", i,lattice[i],j,lattice[j]);
           lattice[i] = 9;
	   lattice[j] = 9;
           naccept++;
	   naccept_dvanni++;
           NumD--;
           NumV--;
	   return;
	}
	return;
      }
      // Neighbor is M, so attempt exchange with the previously selected
      // new orientation for the dummbell.
      if(j_M) {
         i_old = lattice[i];
         j_old = lattice[j];
         einitial = site_energy(i)+site_energy(j);
         i_new = j_old;

         lattice[i] = i_new;
         lattice[j] = j_new;
         edelta = site_energy(i) + site_energy(j) - einitial;

         if(edelta <= 0.0){
           naccept++;
           naccept_nntr++;
           return;
         }
         else if (temperature > 0.0){ 
            P = random->uniform();
	    if(P <= exp(-1*edelta*t_inverse)) {
	       naccept++;
               naccept_nntr++;
	       return;
            }
         }
         else {
	  lattice[i] = i_old;
	  lattice[j] = j_old;
	  return;
        }
     }
   }
 }
  // if i is Hr, attempt exchange with M or mutual anni with D
  else if (i_Hr){
    int nbor = (int) (maxneigh*random->uniform());
    if (nbor >= maxneigh) nbor = maxneigh-1;
    int j = neighbor[i][nbor];
    j_old = lattice[j];
    if (j_old > 8 && j_old < 19) 
       j_M = 1;

    // if neigh is M, attempt exchange
    if (j_M){
        P = random->uniform();
	if(P <= P_Hrdiff) {
          einitial = site_energy(i)+site_energy(j);
          i_new = j_old;
          j_new = i_old;

          lattice[i] = i_new;
          lattice[j] = j_new;
          edelta = site_energy(i) + site_energy(j) - einitial;

          if(edelta <= 0.0){
            naccept++;
//*****            naccept_nntr++;
            return;
          }
          else if (temperature > 0.0){ 
            P = random->uniform();
	    if(P <= exp(-1*edelta*t_inverse)) {
	       naccept++;
//*****               naccept_nntr++;
	       return;
            }
          }
         else {
	  lattice[i] = i_old;
	  lattice[j] = j_old;
	  return;
	 }
       }
     }
  }  
 
//DEBUG
//else if (lattice[i] < 3) printf("Error, lattice[%i] = %i\n", i, lattice[i]);
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppBccSelfdiffusion::site_propensity(int i)
{
  error->all(FLERR,"kMC has been disabled in this version. Use rkMC\n");
  int j,k, i_old, j_old, nbor, eflag;
  double einitial,edelta,probone,proball;

  clear_events(i);

  i_old = lattice[i];
  proball = 0.0;
  probone = 0.0;
  j = 0;

  // if i is a dumbbell
  if(lattice[i] > 2) {
    i_old = lattice[i]; 

    //D can be annihilated at a sink or ...
    //proball += probone = (double) maxneigh * D_anniAtSinks; 
    proball += probone = D_anniAtSinks; 
    add_event(i,-1,probone,ANNIATSINK);
    //printf("D_anni, i = %i, probone = %f\n",i,probone);

    //... or rotate onsite or ...
    //proball += probone = (double) maxneigh * P_rot*(1.0-D_anniAtSinks);
    proball += probone = P_rot*(1.0-D_anniAtSinks);
    add_event(i,-1,probone,ONSITEROT);
    //printf("Rot, i = %i, probone = %f\n",i,probone);

    // ... or translate and rotate.  Each D can translate to 
    // 4 particular neighbors depending on orientation
    for(nbor = 0; nbor < 4; nbor++){
       if(i_old==3){ // [1-10] can rotate to [01-1]=6 & [10-1]=8
          if(nbor == 0) { j = neighbor[i][2];
          } else if(nbor == 1) { j = neighbor[i][3];
          } else if(nbor == 2) { j = neighbor[i][4];
          } else if(nbor == 3) { j = neighbor[i][5];
          } 
       } else if (i_old==4) {
         if(nbor == 0) { j = neighbor[i][0];
         } else if(nbor == 1) { j = neighbor[i][1];
         } else if(nbor == 2) { j = neighbor[i][6];
         } else if(nbor == 3) { j = neighbor[i][7];
         }
      } else if (i_old==5) {
        if(nbor == 0) { j = neighbor[i][0];
        } else if(nbor == 1) { j = neighbor[i][3];
        } else if(nbor == 2) { j = neighbor[i][4];
        } else if(nbor == 3) { j = neighbor[i][7];
        }
      } else if (i_old==6) {
        if(nbor == 0) { j = neighbor[i][1];
        } else if(nbor == 1) { j = neighbor[i][2];
        } else if(nbor == 2) { j = neighbor[i][5];
        } else if(nbor == 3) { j = neighbor[i][6];
        }
      } else if (i_old==7) {
        if(nbor == 0) { j = neighbor[i][0];
        } else if(nbor == 1) { j = neighbor[i][2];
        } else if(nbor == 2) { j = neighbor[i][5];
        } else if(nbor == 3) { j = neighbor[i][7];
        }
      } else if (i_old==8) {
        if(nbor == 0) { j = neighbor[i][1];
        } else if(nbor == 1) { j = neighbor[i][3];
        } else if(nbor == 2) { j = neighbor[i][4];
        } else if(nbor == 3) { j = neighbor[i][6];
        }
      }
      j_old = lattice[j];
      if(j_old == 1) {
	// proball += probone = VD_mutualAnni*(1.0-D_anniAtSinks);
	proball += probone = VD_mutualAnni*(1.0-D_anniAtSinks)/4.0;
	add_event(i,j,probone,ANNI);
        //printf("DV_ann, i = %i, probone = %f\n",i,probone);
      }
      else if(j_old == 2) {
        if(DD == 0 && VV == 0 && DV == 0){
          // probone = P_nntr*(1.0-D_anniAtSinks);
          probone = P_nntr*(1.0-D_anniAtSinks)/4.0;
        }
        else {
          // calculate delta E
          einitial = site_energy(i) + site_energy(j);
          lattice[i] = j_old;
          lattice[j] = i_old;
          edelta = site_energy(i) + site_energy(j) - einitial;
	  lattice[i] = i_old;
	  lattice[j] = j_old;

          // calculate probability of exchange
          //if(edelta <= 0.0) probone = P_nntr;
          if(edelta <= 0.0) probone = P_nntr/4.0;
          else if(temperature == 0.0) probone = 0.0;
          //else probone = exp(-1*edelta*t_inverse)*P_nntr*(1.0-D_anniAtSinks);
          else probone = exp(-1*edelta*t_inverse)*P_nntr*(1.0-D_anniAtSinks)/4.0;
        }
        if(probone > 0.0){
          proball += probone;
          add_event(i,j,probone,NNTRANSROT);
          //printf("NNTR, i = %i, probone = %f\n",i,probone);
	  // DEBUG
          // if(i_old < 3 || j_old != 2)
          // printf("In add_event, l[i=%i] = %i, l[j=%i]=%i, P_nntr = %f, probone = %f, proball = %f\n", i,lattice[i],j,lattice[j],P_nntr,probone,proball);
        }
     }
   }
 }

  // If vacancy, then annihilation at sink or with neighboring D
  // and exchange with neighboring M site are the possible events.
  else if(lattice[i] == 1) {
     // proball += probone = (double) maxneigh * V_anniAtSinks;
     proball += probone = V_anniAtSinks;
     if(probone > 0.0){
       add_event(i,-1,probone,ANNIATSINK);
       //printf("V_anni, i = %i, probone = %f\n",i,probone);
     } 

     // calculate probability of exchange with M neighbors
     for(nbor = 0; nbor < maxneigh; nbor++){
	i_old = lattice[i];
	j = neighbor[i][nbor];
	j_old = lattice[j];

        // if vac and dumbbell, annihilate both at rate=VD_mutualAnni
	if(j_old > 2){
	   // proball += probone = VD_mutualAnni; 
	   proball += probone = VD_mutualAnni/(double) maxneigh; 
	   add_event(i,j,probone,ANNI);
           //printf("DV_anni, i = %i, probone = %f\n",i,probone);
	}
    
        // if vac and M exchange
	else if(j_old == 2) {

	  // if bond energies = 0, no need to calculate delta E, just probability
	  if(DD == 0 && VV == 0 && DV == 0){
            proball += probone = P_vdiff/(double) maxneigh;
	    add_event(i,j,probone,EXCHANGE);
            //printf("EXCH, i = %i, probone = %f\n",i,probone);
          }
	  else{
            // calculate delta E
            einitial = site_energy(i) + site_energy(j);
	    lattice[i] = j_old;
	    lattice[j] = i_old;
            edelta = site_energy(i) + site_energy(j) - einitial;
	    lattice[i] = i_old;
	    lattice[j] = j_old;

            // calculate probablity of exchange
	    probone = 0.0;
	    if(edelta <= 0.0){ 
	      // proball += probone = P_vdiff;;
	      proball += probone = P_vdiff/(double) maxneigh;
 	      add_event(i,j,probone,EXCHANGE);
              //printf("EXCH, i = %i, probone = %f\n",i,probone);
	    }
	    else if (temperature > 0.0) {
	      //proball += probone = P_vdiff*(1.0-V_anniAtSinks)*exp(-1.0*edelta*t_inverse);
	      proball += probone = P_vdiff*(1.0-V_anniAtSinks)*exp(-1.0*edelta*t_inverse)/(double) maxneigh;
 	      add_event(i,j,probone,EXCHANGE);
              //printf("EXCH, i = %i, probone = %f\n",i,probone);
	    }
	  }
	}
     }
  }
  //if(proball > 0.0) printf("i = %i, proball = %f\n",i,proball);
  return proball;
}


/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::site_event(int i, class RandomPark *random)
{
  return site_event_linear(i,random);
}

/* ---------------------------------------------------------------------- */


void AppBccSelfdiffusion::site_event_linear(int i, class RandomPark *random)
{
  int j,k,m,isite,i_old,j_old,i_new,j_new;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;
  int updatejnbors = 0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) {
         printf("i = %i, propensity = %f, threshhold = %f, proball =%f\n", i, propensity[i2site[i]],threshhold, proball);
         error->one(FLERR,"Did not reach event propensity threshhold");
    }
  }

  if (events[ievent].style == ANNIATSINK){
     j = events[ievent].destination;
     // DEBUG
     // if(j != -1 || lattice[i] == 2) 
	// printf("Error, Annihilate at sink.\n");
     if(lattice[i] > 2) {
	naccept_danni++;
        NumD--;
     }	
     else { 
	naccept_vanni++;
        NumV--;
     }
     lattice[i] = 2;
  }

  else if (events[ievent].style == ANNI){
     j = events[ievent].destination;
     // DEBUG
     // if(lattice[i] == 2 || lattice[j] == 2){
         //printf("Error, trying to annihilate incorrectly.");
         //printf(" lattice[%i] = %i, lattice[%i] = %i\n", i,lattice[i], j,lattice[j]);}
     lattice[i] = 2;
     lattice[j] = 2;
     updatejnbors = 1;
     naccept_dvanni++;
     NumD--;
     NumV--;
  }

  else if (events[ievent].style == ONSITEROT) {
     j = events[ievent].destination;
     // DEBUG
     // if(j != -1 || lattice[i] < 3) {
	// printf("Error, trying to rotate on site incorrectly.\n");}
     while(!0){
	i_new = (int) (6*random->uniform() + 3);
	if(i_new > 8 || i_new < 3){ 
	   printf("Error, i_new for onsite rotation is out of bounds\n");}
        i_old = lattice[i];
	if (i_new != i_old) {
	   if (i_old == 3){
	     if (i_new != 4) break;
	   } else if (i_old == 4) {
	     if (i_new != 3) break;
	   } else if (i_old == 5) {
	     if (i_new != 6) break;
	   } else if (i_old == 6) {
	     if (i_new != 5) break;
	   } else if (i_old == 7) {
	     if (i_new != 8) break;
	   } else if (i_old == 8) {
	     if (i_new != 7) break;
	   }
	}
     }
     lattice[i] = i_new;
     naccept_rot++;
  }
  else if (events[ievent].style == NNTRANSROT) {
     j = events[ievent].destination;
     // DEBUG
     // if (lattice[i] < 3 || lattice[j] != 2){ 
	// printf("Error, incorrect NNTRANSROT l[i=%i]=%i, l[j=%i]=%i.\n",i,lattice[i],j,lattice[j]);}

     int i_old = lattice[i]; 

     if (i_old == 3){
        j_new = 6;
	if (random->uniform() <= 0.5) j_new = 8; 
     } else if (i_old==4) {
	j_new = 5;
        if (random->uniform() <= 0.5) j_new = 7;
     } else if (i_old==5) {
        j_new = 4;
        if (random->uniform() <= 0.5) j_new = 7;
     } else if (i_old==6) {
        j_new = 3;
        if (random->uniform() <= 0.5) j_new = 8;
     } else if (i_old==7) {
        j_new = 4;
        if (random->uniform() <= 0.5) j_new = 5;
     } else if (i_old==8) {
        j_new = 3;
        if (random->uniform() <= 0.5) j_new = 6;
     }
     j_old = lattice[j];
     lattice[j] = j_new;
     lattice[i] = j_old;
     updatejnbors = 1;
     naccept_nntr++;
  }
  else if (events[ievent].style == EXCHANGE) {
     j = events[ievent].destination;
     i_old = lattice[i];
     if (lattice[i] == 1 && lattice[j] == 2) 
	naccept_Vnn++;

     // DEBUG
     // else if((lattice[i] == 1 && lattice[j] > 2) || (lattice[i] > 2 && lattice[j] == 1))
	//printf("Error, should not be in this loop, spin[%i]=%i & spin[%i]=%i, event_style=%i\n", i,i_old,j,lattice[j],events[ievent].style);
     // else 
	// printf("Error in EXCHANGE of kMC: i = %i, j = %i\n, events_style=%i", i,j,events[ievent].style);

     lattice[i] = lattice[j];
     lattice[j] = i_old;
     updatejnbors = 1;
  }
  // compute propensity changes for self and swap site and their neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  // update site i's neighbors, this will include the exchanged site

  for (k = 0; k < maxneigh; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    // not quite sure what this does
    if (isite < 0) continue;
    // add to update list
    esites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
    echeck[isite] = 1;
  }

  // update exchanged site's neighbors
  // avoid any that have already been found
  if (updatejnbors == 1){
     for (k = 0; k < maxneigh; k++) {
       m = neighbor[j][k];
       isite = i2site[m];
       // not quite sure what this does
       if (isite < 0) continue;
       // make sure site is not already updated
       if (echeck[isite] == 1) continue;
       // add to update list
       esites[nsites++] = isite;
       propensity[isite] = site_propensity(m);
       echeck[isite] = 1;
     }
  }

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (k = 0; k < nsites; k++) echeck[esites[k]] = 0;
 // printf("Exit perform kMC\n");
}


/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::add_event(int i, int destination,
			      double propensity, int eventflag)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events =
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;
  events[freeevent].propensity = propensity;
  events[freeevent].destination = destination;
  events[freeevent].style = eventflag;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppBccSelfdiffusion::allocate_data()
{
  // for linear:
  //   make esites large enough for 1 sites and their 1,2 neighbors

  if (engstyle == LINEAR) {
    int emax = 1 + maxneigh*2;
    esites = new int[2*emax];
  }

  echeck = new int[nlocal+nghost];

  memory->create(firstevent,nlocal,"app:firstevent");
}
