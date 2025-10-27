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
#include "app_diffusion_multiphase.h"

using namespace SPPARKS_NS;
using std::map;
using std::set;

enum{LINEAR};
enum{ANNIATSINK,ANNI,ONSITEROT,TRANS,NNTRANSROT,NNNTRANSROT,EXCHANGE};

#define DELTAEVENT 100000

// This app is based on the diffusion app (for Kawasaki dynamics)
// These are the significant changes and simplifications
// 1. No Schowebel hops, only exchanges to site neighbors
// 2. Allow for an arbitrary number of phases (species for an atomic model)
// 3. Only the linear energy style is supported

/* ---------------------------------------------------------------------- */

AppDiffusionMultiphase::AppDiffusionMultiphase(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg), phase_labels()
{
  // need to double check these values

  ninteger = 2;
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 0;
  numrandom = 1;

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
}

/* ---------------------------------------------------------------------- */

AppDiffusionMultiphase::~AppDiffusionMultiphase()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->destroy(firstevent);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::input_app(char *command, int narg, char **arg)
{
  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(FLERR,str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

  if (strcmp(command,"diffusion/multiphase") == 0)
    parse_diffmultiphase(narg,arg);
  else error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppDiffusionMultiphase::parse_diffmultiphase(int narg, char **arg)
{
   // 2 args: diffusion/multiphase phase <int value>
   // 2 args: diffusion/multiphase rate <double value>
   // 5 args: diffusion/multiphase energy <double value> 

   if (narg < 2)
     error->all(FLERR,"Illegal diffusion/multiphase command");

   if (strcmp(arg[0],"phase")==0){
      if (narg != 2)
        error->all(FLERR,"Illegal diffusion/multiphase phase command: "
                   "num args != 2");
      int phase = std::atoi(arg[1]);
      phase_labels.insert(phase);
      if (phase < 1)
        error->all(FLERR,"Illegal diffusion/multiphase phase value: "
                   "must be >= 1");
   } else if (strcmp(arg[0], "rate") == 0) {
       if (narg != 3) {
            error->all(FLERR, "Illegal diffusion/multiphase rate command: "
                               "Expected 'rate <name> <value>'");
        } 
        std::string rate_name = arg[1]; // The name of the rate variable (e.g., "V_anniAtSinks")
        double rate_value = std::atof(arg[2]); // Using std::atof as requested

        // Assign the value to the corresponding member variable
        if (rate_name == "V_anniAtSinks") {
            V_anniAtSinks = rate_value;
        } else if (rate_name == "D_anniAtSinks") {
            D_anniAtSinks = rate_value;
        } else if (rate_name == "Rrot") {
            Rrot = rate_value;
        } else if (rate_name == "RnnTrans") {
            RnnTrans = rate_value;
        } else if (rate_name == "VD_mutualAnni") {
            VD_mutualAnni = rate_value;
        } else {
            error->all(FLERR, ("Unknown diffusion/multiphase rate name: " + rate_name).c_str());
        }
    } else if (strcmp(arg[0], "VtoDdiffcoeff") == 0) {
       if (narg != 3) {
            error->all(FLERR, "Illegal diffusion/multiphase rate command: "
                               "Expected 'diffcoeff <name> <value>'");
        } 
        std::string diffcoeff = arg[1]; // The name of diffusivity variable (e.g., "VtoD_DiffCoeff")
        double diffcoeff_value = std::atof(arg[2]); // Using std::atof as requested

        // Assign the value to the corresponding member variable
        if (diffcoeff == "VtoD_DiffCoeff") {
            VtoD_DiffCoeff = diffcoeff_value;
        } else {
            error->all(FLERR, ("Unknown diffusion/multiphase diffcoeff name: "+std::to_string(diffcoeff_value)).c_str());
        }
        if (diffcoeff_value > 1.0 || diffcoeff_value <= 0.0) 
           error->all(FLERR,"Illegal diffusion/multiphase diffcoeff value: "
                   "must be > 0 and <= 1");
    } 
    else if (strcmp(arg[0], "energy") == 0) {
       if (narg != 3) {
            error->all(FLERR, "Illegal diffusion/multiphase rate command: "
                               "Expected 'energy <name> <value>'");
        } 
        std::string energy = arg[1]; // The name of the energy variable (e.g., "V-V")
        double energy_value = std::atof(arg[2]); // Using std::atof as requested

        // Assign the value to the corresponding member variable
        if (energy == "VV") {
            VV = energy_value;
        } else if (energy == "DD") {
            DD = energy_value;
        } else if (energy == "DV") {
            DV = energy_value;
        } else {
            error->all(FLERR, ("Unknown diffusion/multiphase energy name: "+std::to_string(energy_value)).c_str());
        }
    }
    else {
        error->all(FLERR, "Illegal diffusion/multiphase command: "
                           "Expected 'phase', 'pin', 'energy', or 'rate'");
    }
} 


/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::grow_app()
{
  lattice = iarray[0];  //vac=1, reg=2, dumbbell=3 - 8: i1
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::init_app()
{
   if (!allocated) allocate_data();
   allocated = 1;

   dimension = domain->dimension;
   dt_sweep = 1.0/maxneigh;

   {
     // insure all site values are in set of phase labels, otherwise error

     std::set<int>::iterator not_found=phase_labels.end();
     int flag = 0;
     for (int i = 0; i < nlocal; i++) {
       if (not_found == phase_labels.find(lattice[i])) flag=1;
     }
     int flagall;
     MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
     if (flagall) error->all(FLERR,"One or more sites have invalid values");
   }
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::setup_app()
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

double AppDiffusionMultiphase::site_energy(int i)
{
  // energy of site = sum of bond weights
  // default values are JDD = 2, JVV = 1, JVD = -1,   


   double energy = 0.0;
   int ip = lattice[i];
   for (int j = 0; j < numneigh[i]; j++){
      int nj = neighbor[i][j];
      int jp = lattice[nj];
      if (ip>2 && jp>2) //if neighboring dumbbells, repel
          energy += DD;
      else if (ip==1 && jp==1) // if neighboring vacancies, repel
          energy += VV;
      else if ((ip>2 && jp==1) || (ip==1 && jp>2)) //if dumbbell and vac, attract.
          energy += -DV;
      }
   return 0.5*energy;
}

/* ----------------------------------------------------------------------
   rKMC Method
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::site_event_rejection(int i, RandomPark *random)
{

  // This app assigns spin = 1 for Vac, 2 for Metal, 3 to 8 for dumbbells in
  // 6 <110> orientations with all 8 spins on regular BCC sites.
  // It attempts diffusion of vacancies by exchanging with M and annihilates
  // both vacancy and dumbbell when they meet.
  // It also diffuses both dumbbells and M.  M can exchange with Vac
  // Dumbbells can rotate onsite, translate to nn, translate to nn & rotate,
  // or translate to nnn and rotate by 90 degrees.
  // Both dumbbell interstitials and vacancies can be annihilated at sinks
  // at a specified rate.  Note, no mass conservation is enforced, they
  // both simply become regular Metal sites.

  double einitial,edelta;
  int i_old, j_old, i_new, j_new;;
  int nbor;
  double P;

  //if i is dumbbell, randomly pick P=[0,1]
  //if P<Rrot, rotate to anther randomly selected orienation

  //else attempt translation and rotation.
  //if Rrot+Rtrn>0, Error.

  if (lattice[i] > 2){
    P = random->uniform();
    if (P <= D_anniAtSinks) { //annihilation of dumbbells at sinks
	lattice[i] = 2;
	naccept++;
	return;
    }
    if (P<Rrot) {
	while(!0){
	  i_new = (int) (6*random->uniform() + 3);
	  if(i_new > 8) printf("i_new > 8 was randomly selected\n");
	  if(i_new == 2) printf("i_new = 2 was randomly selected\n");
	  if (i_new != lattice[i]) break;
        }
          //printf("Onsite rotation, old_spin[%i] = %i, ",i,lattice[i]);
	  lattice[i] = i_new;
          naccept++;
          return;
    }

    //else attempt translation or translation + rotation
    else {
      nbor = (int) (maxneigh*random->uniform());
      if (nbor >= maxneigh) nbor = maxneigh-1;
      int j = neighbor[i][nbor];

      //if neighbor is another dumbbell, return
      if (lattice[j] > 2){
	  //printf("Dumbbell has Dumbbell neighbor\n");
	 return;
      }

      //if neighbor is a vacancy, i&j become M regular sites
      if (lattice[j] == 1){
        P = random->uniform();
	if(P <= VD_mutualAnni) {
           //printf("annihilate, old_spin_i[%i] = %i, old_spin_j[%i] = %i, ", i,lattice[i],j,lattice[j]);
           lattice[i] = 2;
	   lattice[j] = 2;
           naccept++;
           //printf("newspin_i = %i, newspin_j = %i\n", lattice[i],lattice[j]);
	   return;
	}
      }

      i_old = lattice[i];
      j_old = lattice[j];
      einitial = site_energy(i)+site_energy(j);

      //if neighbor is M on nn, then translation or translation+rotation
      if (nbor < 8) {
        //can dumbbell translate to that neighboring site, nbor
	if(i_old==3){ // orientation of the dumbbell is (1-10)
	   if(nbor==0 || nbor==1 || nbor==6 || nbor==7){
	      return;}}
	else if(i_old==4){ // (110) 
	   if(nbor==2 || nbor==3 || nbor==4 || nbor==5){
	      return;}}
	else if(i_old==5){ // (011)
	   if(nbor==1 || nbor==2 || nbor==5 || nbor==6){
	      return;}}
	else if(i_old==6){ // 01-1)
	   if(nbor==0 || nbor==3 || nbor==4 || nbor==7){
	      return;}}
	else if(i_old==7){ // (101)
	   if(nbor==1 || nbor==3 || nbor==4 || nbor==6){
	      return;}}
	else if(i_old==8){ // (10-1)
	   if(nbor==0 || nbor==2 || nbor==5 || nbor==7){
	      return;}}

        // translate with probability of 0.5.
        P = random->uniform();
   	if (P<0.5) { // Rtrn = Rtrn_rot = (1-Rrot)/2; 
           i_new = j_old;
   	   j_new = i_old;
           //printf("attempting translation, i_old[%i] = %i, neigh = %i, j_old[%i] = %i, i_new = %i, j_new = %i\n", i,i_old,nbor,j,j_old,i_new,j_new);
	}
        //translate and rotate
	else {
	  while(!0){
	    j_new = (int) (6*random->uniform() + 3);
	    if(j_new == 2) printf("j_new = 2 was randomly selected\n");
	    if(j_new > 8) printf("j_new > 8 was randomly selected\n");
	    //printf("In trns+rot, j_new = %i, i_old = %i\n", j_new, i_old);
	    if (j_new != i_old) break;
	  }
	  i_new = j_old;
          //printf("attempting trans+rot, i_old[%i] = %i, neigh = %i, j_old[%i] = %i, i_newi = %i, j_new = %i\n", i,i_old,nbor,j,j_old,i_new,j_new);
	}
      }
      else { // nbor > 7 and M neighbor is nnn, then translate and 90 degree rotation 
        // can dumbbell translate to that neighbor site, nbor
	if(i_old==3 || i_old==4) { // orientation of the dumbbell is (1-10) (110)
	   if(nbor==10 || nbor==11){
	      return;}
	   else{
	      if (i_old==3) j_new = 4;
	      else j_new = 3;
	   }
	}
	if(i_old==5 || i_old==6) { // (011) (01-1)
	   if(nbor==8 || nbor==13){
	      return;}
	   else {
	     if (i_old==5) j_new = 6;
	     else j_new = 5;
	   }
	} 
	if(i_old==7 || i_old==8) { // (101) (10-1)
	   if(nbor==9 || nbor==12){
	      return;}
	   else {
	      if (i_old==7) j_new = 8;
	      else j_new = 7;
	   }
	}
	i_new = j_old;
        //printf("nnn trans+90rot, i_old[%i] = %i, neigh = %i, j_old[%i] = %i, i_new = %i, j_new = %i\n", i,i_old,nbor,j,j_old,i_new,j_new);
      }
      lattice[i] = i_new;
      lattice[j] = j_new;
      edelta = site_energy(i) + site_energy(j) - einitial;

      if(edelta <= 0.0){
        naccept++;
	//printf("attempt successful\n");
        return;
      }
      else if (temperature > 0.0){ 
        P = random->uniform();
	if(P < exp(-1*edelta*t_inverse)) {
	   naccept++;
	   //printf("attempt successful\n");
	   return;
        }
        else{
	  lattice[i] = i_old;
	  lattice[j] = j_old;
	  //printf("attempt UNsuccessful\n");
	  return;
        }
      }
    }
  }
  // if vac, and neigh is D, annihilate
  else if (lattice[i] == 1){
    P = random->uniform();
    if (P <= V_anniAtSinks) { //annihilation of vacancies at sinks
	lattice[i] = 2;
	naccept++;
	return;
    }
      int nbor = (int) (maxneigh*random->uniform());
      if (nbor > maxneigh) nbor = maxneigh-1;
      int j = neighbor[i][nbor];
      
      if (lattice[i] == lattice[j]) return;
      else if (lattice[j] > 2){
        P = random->uniform();
	if(P <= VD_mutualAnni){
           //printf("Vac-D ANNIHILATED: i_old[%i] = %i, j_old[%i] = %i, ", i,lattice[i],j,lattice[j]);
           lattice[i] = 2;
           lattice[j] = 2;
           naccept++;
           //printf("i_new = %i, j_new = %i\n", lattice[i],lattice[j]);
           return;
	}
      }
      else { // vac & M 
  	i_old = lattice[i];
  	j_old = lattice[j];
	einitial = site_energy(i)+site_energy(j);

        lattice[i] = j_old;
   	lattice[j] = i_old;
        edelta = site_energy(i)+site_energy(j) - einitial;

	if(edelta <= 0.0){
	  naccept++;
          //printf("Vac-M exhanged: i_old[%i] = %i, j_old[%i] = %i, i_new = %i, j_new = %i\n", i,i_old,j,j_old,lattice[i],lattice[j]);
	  return;
        }
        else if (temperature > 0.0){ 
          P = random->uniform();
	  if(P <= exp(-1*edelta*t_inverse)) {
	    naccept++;
            //printf("Vac-M exhanged: i_old[%i] = %i, j_old[%i] = %i, i_new = %i, j_new = %i\n", i,i_old,j,j_old,lattice[i],lattice[j]);
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
  else { // Reg M can only exchange with vac.  This is to simplify the coding of determining 
	 // if the neighbor dumbbell's orientation allows an exchange.  Additionally, it is
	 // an uncommon event.

      int nbor = (int) (maxneigh*random->uniform());
      if (nbor > maxneigh) nbor = maxneigh-1;
      int j = neighbor[i][nbor];
      if(lattice[i] == lattice[j]) return;

      if(lattice[j] == 1){
  	i_old = lattice[i];
  	j_old = lattice[j];
	einitial = site_energy(i)+site_energy(j);
 
        lattice[i] = j_old;
   	lattice[j] = i_old;
        edelta = site_energy(i)+site_energy(j) - einitial; 

	if(edelta <= 0.0){
	  naccept++;
          //printf("Vac-M exhanged: old_spin_i[%i] = %i, old_spin_j[%i] = %i, newspin_i = %i, newspin_j = %i\n", i,i_old,j,j_old,lattice[i],lattice[j]);
	  return;
	}
        else if (temperature > 0.0){ 
          P = random->uniform();
	  if(P <= exp(-1*edelta*t_inverse)) {
	    naccept++;
            //printf("Vac-M exhanged: i_old[%i] = %i, j_old[%i] = %i, i_new = %i, j_new = %i\n", i,i_old,j,j_old,lattice[i],lattice[j]);
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

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusionMultiphase::site_propensity(int i)
{
  return site_propensity_linear(i);
}

/* ---------------------------------------------------------------------- */

double AppDiffusionMultiphase::site_propensity_linear(int i)
{
  int j,k, i_old, j_old, eflag;
  double e0,einitial,edelta,probone,proball;

  clear_events(i);

  i_old = lattice[i];
  proball = 0.0;
  probone = 0.0;

  // if i is a dumbbell
  if(lattice[i] > 2) {

     i_old = lattice[i]; 
     proball += probone = D_anniAtSinks; //D can be annihilated at a sink or ...
     add_event(i,-1,D_anniAtSinks,ANNIATSINK);

     proball += probone = Rrot; //... or rotate onsite or ...
     add_event(i,-1,Rrot,ONSITEROT);

     for(j = 0; j < numneigh[i]; j++){
           i_old = lattice[i];
	   j_old = lattice[neighbor[i][j]]; 

	// ...or annihilate w/ vacancy or ...
	if(lattice[neighbor[i][j]] == 1) { 
	   proball += probone = (1.0-D_anniAtSinks-Rrot)*VD_mutualAnni; 
	   add_event(i,neighbor[i][j],probone,ANNI);

           if(lattice[i] < 3){
	      printf("i = %i neigh = %i, and should not be in dumbbell loop\n", lattice[i], j);
	   }
           if ((lattice[i] > 2 && j_old == 1) || (lattice[i] == 1 && j_old > 2)){}
	   else {
		printf("Error in propensity, ANNI for D, lattice[%i]=%i, neigh=%i, lattice[n]=%i\n", i,lattice[i],j,j_old);
	   }
        }

	// ... or translate or translate+rotate w/ nn or translate+rot w/nnn
	else if(lattice[neighbor[i][j]] == 2) {
           einitial = site_energy(i) + site_energy(j);
	   lattice[i] = j_old;
	   lattice[neighbor[i][j]] = i_old;
           edelta = site_energy(i) + site_energy(j) - einitial;

 	   // ... or translate or translate + rotate with nn
	   if(j < 8){
	      int exchange = 1;
	      if(i_old==3){ // orientation of the dumbbell is (1-10)
           	if(j==0 || j==1 || j==6 || j==7){
                   exchange = 0;}}
       	      else if(i_old==4){ // (110)
             	if(j==2 || j==3 || j==4 || j==5){
                   exchange = 0;}}
              else if(i_old==5){ // (011)
           	if(j==1 || j==2 || j==5 || j==6){
                   exchange = 0;}}
              else if(i_old==6){ // 01-1)
           	if(j==0 || j==3 || j==4 || j==7){
                   exchange = 0;}}
              else if(i_old==7){ // (101)
           	if(j==1 || j==3 || j==4 || j==6){
                   exchange = 0;}}
              else if(i_old==8){ // (10-1)
           	if(j==0 || j==2 || j==5 || j==7){
                   exchange = 0;}}
	      if(exchange == 1){
	         if(edelta <= 0.0){
 		    //0.5 b/c it can exchage with 4/8 nn, 0.5 for trans & transrot
		    proball += probone = 0.5*RnnTrans*(1.0-D_anniAtSinks-Rrot);
	   	    add_event(i,neighbor[i][j],probone,TRANS);
		    proball += probone = 0.5*(1.0-RnnTrans)*(1.0-D_anniAtSinks-Rrot);
	   	    add_event(i,neighbor[i][j],probone,NNTRANSROT);
	      	 }
	         else if (temperature > 0.0){
		    proball += probone = 0.5*RnnTrans*(1.0-D_anniAtSinks-Rrot)*
			                 exp(-1.0*edelta*t_inverse);
	   	    add_event(i,neighbor[i][j],probone,TRANS);
		    proball += probone = 0.5*(1.0-RnnTrans)*(1.0-D_anniAtSinks-Rrot)*
					 exp(-1.0*edelta*t_inverse);
	   	    add_event(i,neighbor[i][j],probone,NNTRANSROT);
	      	 }
	      }
	   }
	   // ... or translate + rotate 90 degrees with nnn
	   else { 
	      int exchange = 1;
	      if(i_old==3 || i_old==4) { // orientation of the dumbbell is (1-10) (110)
           	if(j==10 || j==11){
              	  exchange = 0;}}
              if(i_old==5 || i_old==6) { // (011) (01-1)
           	if(j==8 || j==13){
              	  exchange = 0;}}
              if(i_old==7 || i_old==8) { // (101) (10-1)
                if(j==9 || j==12){
                  exchange = 0;}}
	      if(exchange == 1){
	         if(edelta <= 0.0){ 
		    // it can exchange with 4/6 nnn
		    proball += probone = 2.0/3.0*(1.0-D_anniAtSinks-Rrot); 
	   	    add_event(i,neighbor[i][j],probone,NNNTRANSROT);
		 }
	         else if (temperature > 0.0){
		    proball += probone = 2.0/3.0*(1.0-D_anniAtSinks-Rrot)
					 *exp(-1.0*edelta*t_inverse);
	   	    add_event(i,neighbor[i][j],probone,NNNTRANSROT);
		 }
	      }
	   }
	   lattice[i] = i_old;
	   lattice[neighbor[i][j]] = j_old;
	}
	else{
	   if(lattice[neighbor[i][j]] < 3) 
	      printf("Error, D skipped anni and trans/transrot incorrectly\n");
        }
     } 
  }
  else if(lattice[i] == 1) {
     proball += V_anniAtSinks;
     add_event(i,-1,V_anniAtSinks,ANNIATSINK);
     for(j = 0; j < numneigh[i]; j++){
	i_old = lattice[i];
	j_old = lattice[neighbor[i][j]];

        // if vac and dumbbell, annihilate both 
	if(lattice[neighbor[i][j]] > 2){
	   proball += probone = 1.0*VD_mutualAnni; 
	   add_event(i,neighbor[i][j],probone,ANNI);
           if (lattice[i] != 1 && j_old < 3){ 
	      printf("Error in propensity for ANNI. i = %i, j = %i\n", lattice[i], j_old);
	   }
	}
        // if vac and M, calculate probability of exchange
	else if(lattice[neighbor[i][j]] == 2) {
           einitial = site_energy(i) + site_energy(j);
	   lattice[i] = j_old;
	   lattice[neighbor[i][j]] = i_old;
           edelta = site_energy(i) + site_energy(j) - einitial;
	   probone = 0.0;
	   if(edelta <= 0.0){ proball += probone = 1.0 - V_anniAtSinks;}
	   else if (temperature > 0.0) {proball += probone = exp(-1.0*edelta*t_inverse);}
	   if(probone > 0.0) {
		add_event(i,neighbor[i][j],probone,EXCHANGE);
	   }
	}
	lattice[i] = i_old;
	lattice[neighbor[i][j]] = j_old;
     }
  }
  else { // lattice[i] == 2 and is M
     for(j = 0; j < numneigh[i]; j++){
	// if neighbor is vac, calculate probability of exchange
	if(lattice[neighbor[i][j]] == 1) {
	   i_old = lattice[i];
	   j_old = lattice[neighbor[i][j]];
           einitial = site_energy(i) + site_energy(j);
	   lattice[i] = j_old;
	   lattice[neighbor[i][j]] = i_old;
           edelta = site_energy(i) + site_energy(j) - einitial;
	   probone = 0.0;
	   if(edelta <= 0.0){ proball += probone = 1.0;}
	   else if (temperature > 0.0){ proball += probone = exp(-1.0*edelta*t_inverse);}
	   if(probone > 0.0) add_event(i,neighbor[i][j],probone,EXCHANGE);
	   lattice[i] = i_old;
	   lattice[neighbor[i][j]] = j_old;
	}
     }
  }
  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::site_event(int i, class RandomPark *random)
{
  return site_event_linear(i,random);
}

/* ---------------------------------------------------------------------- */

void AppDiffusionMultiphase::site_event_linear(int i, class RandomPark *random)
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
    if (ievent < 0) error->one(FLERR,"Did not reach event propensity threshhold");
  }

  if (events[ievent].style == ANNIATSINK){
     j = events[ievent].destination;
     if(j != -1 || lattice[i] == 2) printf("Error, Annihilate at sink.\n");
     lattice[i] = 2;
  }

  else if (events[ievent].style == ANNI){
     j = events[ievent].destination;
     if ((lattice[i] > 2 && lattice[j] == 1) || (lattice[i] == 1 && lattice[j] > 2)){
        i_old = lattice[i]; j_old = lattice[j];
	lattice[i] = 2;
	lattice[j] = 2;
	updatejnbors = 1;
     }
     else {printf("Error, trying to annihilate incorrectly.");
           printf(" lattice[%i] = %i, lattice[%i] = %i\n", i,lattice[i], j,lattice[j]);
     }
  }

  else if (events[ievent].style == ONSITEROT) {
     j = events[ievent].destination;
     if(j != -1 || lattice[i] < 3) {
	printf("Error, trying to rotate on site incorrectly.\n");}
     while(!0){
	i_new = (int) (6*random->uniform() + 3);
	if(i_new > 8 || i_new < 3){ 
	   printf("Error, i_new for onsite rotation is out of bounds\n");}
	if (i_new != lattice[i]) break;
     }
     lattice[i] = i_new;
  }
  else if (events[ievent].style == TRANS) {
     j = events[ievent].destination;
     // assumes i is a dumbbell and j is M
     if (lattice[i] < 3 || lattice[j] != 2) printf("Error, incorrect TRANS.\n");
     lattice[j] = lattice[i];
     lattice[i] = 2;
     updatejnbors = 1;
  }
  else if (events[ievent].style == NNTRANSROT) {
     j = events[ievent].destination;
     if (lattice[i] < 3 || lattice[j] != 2){ 
	printf("Error, incorrect NNTRANSROT.\n");}
     while(!0){
	j_new = (int) (6*random->uniform() + 3);
	if(j_new > 8 || j_new < 3) printf("Error, j_new for trans+rotation is out of bounds\n");
	if (j_new != lattice[i]) break;
     }
     lattice[j] = j_new;
     lattice[i] = 2;
     updatejnbors = 1;
  }
  else if (events[ievent].style == NNNTRANSROT) {
     j = events[ievent].destination;
     if (lattice[i] < 3 || lattice[j] != 2) {
	printf("Error, incorrect NNNTRANSROT. i = %i, j = %i\n", lattice[i], lattice[j]);}
     i_old = lattice[i];
     if (i_old == 3) j_new = 4; // (1-10)
     else if (i_old == 4) j_new = 3;
     else if (i_old == 5) j_new = 6;
     else if (i_old == 6)j_new = 5;
     else if (i_old==7) j_new = 8;
     else j_new = 7;
     lattice[j] = j_new;
     lattice[i] = 2;
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

  for (k = 0; k < numneigh[i]; k++) {
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
     for (k = 0; k < numneigh[j]; k++) {
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
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::clear_events(int i)
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

void AppDiffusionMultiphase::add_event(int i, int destination,
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

void AppDiffusionMultiphase::allocate_data()
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
