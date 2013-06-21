#include "pair_eam_alloy_spline.h"
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"

#include <sstream>
#include <fstream>
#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMAlloySpline::PairEAMAlloySpline(LAMMPS *lmp) : Pair(lmp), cutoff_max(0.0)
{
  single_enable = 0;        // 1 if single() routine exists
  restartinfo = 0;          // 1 if pair style writes restart info
  one_coeff = 1;            // 1 if allows only one coeff * * call

  comm_forward = 0;         // size of forward communication (0 if none)
  comm_reverse = 0;         // size of reverse communication (0 if none)
}

/* ---------------------------------------------------------------------- */

PairEAMAlloySpline::~PairEAMAlloySpline()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ----------------------------------------------------------------------
   Compute energy, stresses, and forces
------------------------------------------------------------------------- */

void PairEAMAlloySpline::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag, vflag);
  else evflag = vflag_fdotr = eflag_global = vflag_global = eflag_atom = vflag_atom = 0;

  double cutoff_maxsq = cutoff_max * cutoff_max;

  double **x = atom->x;
  double **F = atom->f;
  int *type = atom->type;
  int ntypes = atom->ntypes;
  int nlocal = atom->nlocal;
  bool newton_pair = force->newton_pair;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // Determine the maximum number of neighbors a single atom has
  int newMaxNeighbors = 0;
  for(int ii = 0; ii < inum; ++ii) {
    int jnum = numneigh[ilist[ii]];
    if(jnum > newMaxNeighbors) newMaxNeighbors = jnum;
  }

  // Allocate array for temporary bond info
  int maxNeighbors = two_body_pair_info.size();
  if(newMaxNeighbors > maxNeighbors) {
    two_body_pair_info.resize(newMaxNeighbors);
  }


  // Loop over full neighbor list of my atoms
  for(int ii = 0; ii < inum; ++ii) {
    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    int itype = type[i]-1;      // shifted back for reference to arrays of Splines
    double rhovalue = 0;

    int n2 = 0; // keep count of 2-body pair info for atom ii

    // Two-body interactions
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; ++jj) {    //loop over neighbors
      int j = jlist[jj];
      int jtype = type[j]-1;      // species of current neighbor. shifted back for reference to arrays
      int phi_col = itype*ntypes + jtype;
      double delx = x[j][0] - xtmp;
      double dely = x[j][1] - ytmp;
      double delz = x[j][2] - ztmp;
      double rsq = delx*delx + dely*dely + delz*delz;

      // Pair potential terms phi
      if (rsq < cutoff_maxsq) {
        double r = std::sqrt(rsq);
        double recip = 1.0/r;

        if (r < phi[phi_col].get_cutoff()) {
          // Compute phi(r_ij) and its gradient in one step
          double phigrad;
          double phival = 0.5 * phi[phi_col].splint_comb(r, phigrad);
          // note that atom types begin with 1, while vector elements begin with 0, hence jtype-1

          // Only half of the gradient contributes to the force as
          // well as half of the energy since we are double counting
          phigrad *= 0.5;

          // Compute tmp force values
          double tmp_force = phigrad * recip;

          // Add in force on atom i from atom j
          // Subtract off force on atom j from atom i
          // Newton's law: action = -reaction
          F[i][0] += tmp_force * delx;
          F[i][1] += tmp_force * dely;
          F[i][2] += tmp_force * delz;
          F[j][0] -= tmp_force * delx;
          F[j][1] -= tmp_force * dely;
          F[j][2] -= tmp_force * delz;

          // Add in piece contributed by neighbor to energy
          if (evflag) ev_tally(i,j,nlocal,newton_pair,phival,0.0,-tmp_force,delx,dely,delz);
        } // rsq < phi_xcutsq

        // Pair density terms
        if (r < rho[jtype].get_cutoff()) {
          rhovalue += rho[jtype].splint_comb(r, two_body_pair_info[n2].df); // Value of the density at atom i

          two_body_pair_info[n2].idx = j;
          two_body_pair_info[n2].r = r;
          two_body_pair_info[n2].recip = 1.0/r;
          two_body_pair_info[n2].nx = delx*recip;
          two_body_pair_info[n2].ny = dely*recip;
          two_body_pair_info[n2].nz = delz*recip;

          ++n2; // increment number of 2-body info's for this atom ii
        }

      } // rsq < max_xcutsq
    } // END LOOP OVER NEIGHBORS jj


    // Done with calculating rho[i], now we can calculate the embedding energy
    double du;
    double embedding_energy = u[itype].splint_comb(rhovalue, du);
    if(eflag) {
      if(eflag_global) eng_vdwl += embedding_energy;
      if(eflag_atom) eatom[i] += embedding_energy;
    }

    // Forces from embedding potential and pair density term
    for (int jj = 0; jj < n2; ++jj) {
      const EAM2Body& neigh_jj = two_body_pair_info[jj];
      int j = neigh_jj.idx;
      double tmp_force = neigh_jj.df * du;
      F[i][0] += neigh_jj.nx * tmp_force;
      F[i][1] += neigh_jj.ny * tmp_force;
      F[i][2] += neigh_jj.nz * tmp_force;

      F[j][0] -= neigh_jj.nx * tmp_force;
      F[j][1] -= neigh_jj.ny * tmp_force;
      F[j][2] -= neigh_jj.nz * tmp_force;
/*
      // WARNING: ADD THIS SECTION LATER WITH CORRECT PARAMETERS!
      if (vflag_either) {
        double delx = neigh_jj.nx * neigh_jj.r;
        double dely = neigh_jj.ny * neigh_jj.r;
        double delz = neigh_jj.nz * neigh_jj.r;
        ev_tally(i,j,nlocal,newton_pair,0.0,0.0,-tmp_force*neigh_jj.recip,delx,dely,delz);
      }
*/

    } // LOOP OVER NEIGHBORS JJ

  } // LOOP OVER ATOMS II

  if(vflag_fdotr) virial_fdotr_compute();

  return;
}

/* ---------------------------------------------------------------------- */

void PairEAMAlloySpline::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   Global settings
------------------------------------------------------------------------- */

void PairEAMAlloySpline::settings(int narg, char **arg)
{
  if(narg != 0) error->all(FLERR,"Illegal pair_style command");
  return;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEAMAlloySpline::coeff(int narg, char **arg)
{
  std::vector<std::string> args(narg);
  for (int i = 0; i < narg; ++i) args[i] = arg[i];

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // ensure I,J args are * *
  if (args[0] != "*" || args[1] != "*")
    error->all(FLERR,"Incorrect args for pair coefficients");

  // WARNING: This is still unnecessary
  /*
  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  map_atom_type = std::vector<std::string>(atom->ntypes, "");
  elements = std::map<std::string, int>();

  int nelements = 0;
  for (int i = 3; i < narg; ++i) {
    std::string element_name = args[i];

    // atom_type -> element_name
    map_atom_type[i-3] = element_name;

    if (element_name == "NULL") {
      elements[element_name] = -1;
      continue;
    }

    // element_name -> potential_idx
    std::map<std::string, int>::iterator it;
    it = elements.find(element_name);
    if (it == elements.end()) { // not found, count this as new element
      elements[element_name] = nelements;
      ++nelements;
    }
  }
  */

  // read potential file
  read_file(args[2]);

  // FROM PREVIOUS WARNING
  /*
  // clear setflag since coeff() called once with I,J = * *
  int n = atom->ntypes;
  for (int i = 1; i <= n; ++i)
    for (int j = i; j <= n; ++j)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  int count = 0;
  for (int i = 1; i <= n; ++i)
    for (int j = i; j <= n; ++j) {
      setflag[i][j] = 1;
      ++count;
    }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
  */
}

/* ----------------------------------------------------------------------
   Init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMAlloySpline::init_style()
{
  if(force->newton_pair == 0)
    error->all(FLERR,"Pair style eam/spline requires newton pair on");

  // Need full neighbor list.
  int irequest_full = neighbor->request(this);
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
}

/* ----------------------------------------------------------------------
   Init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAMAlloySpline::init_one(int i, int j)
{
  return cutoff_max;
}

/* ----------------------------------------------------------------------
   Set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEAMAlloySpline::read_file(std::string filename)
{
  if(comm->me == 0) {
    std::ifstream ifs;
    ifs.open(filename.c_str(), std::ifstream::in);

    if(!ifs) {
      std::stringstream oss;
      oss << "Cannot open spline EAM potential file " << filename;
      std::string error_str = oss.str();
      error->one(FLERR,error_str.c_str());
    }

    // Skip first line of file
    std::string tmp_line;
    std::getline(ifs, tmp_line);

    // Read in potentials
    int ntypes = atom->ntypes;

    // Read in phi: phi_aa, phi_ab, phi_ba, phi_bb
    int nphi = ntypes * ntypes;
    phi.resize(nphi, Spline(lmp));
    for (int i = 0; i < ntypes; ++i) {
      for (int j = i; j < ntypes; ++j) {
        ifs >> phi[i*ntypes + j];

        // keep symmetry: phi_ij = phi_ji
  if ( i != j ) phi[j*ntypes + i] = phi[i*ntypes + j];
      }
    }

    // Read in rho: rho_a, rho_b
    int nrho = ntypes;
    rho.resize(nrho, Spline(lmp));
    for (int i = 0; i < nrho; ++i) {
      ifs >> rho[i];
    }

    // Read in U: U_a, U_b
    int nu = ntypes;
    u.resize(nu, Spline(lmp));
    for (int i = 0; i < nu; ++i) {
      ifs >> u[i];
    }

    ifs.close();
  }

  // Communicate potentials
  for (int i = 0; i < phi.size(); ++i) phi[i].communicate();
  for (int i = 0; i < rho.size(); ++i) rho[i].communicate();
  for (int i = 0; i < u.size();   ++i) u[i].communicate();

  // Determine maximum cutoff radius of all relevant spline functions
  cutoff_max = 0.0;
  for (int i = 0; i < phi.size(); ++i)
    cutoff_max = std::max(cutoff_max, phi[i].get_cutoff());
  for (int i = 0; i < rho.size(); ++i)
    cutoff_max = std::max(cutoff_max, rho[i].get_cutoff());

  // Set LAMMPS pair interaction flags
  for(int i = 1; i <= atom->ntypes; ++i) {
    for(int j = 1; j <= atom->ntypes; ++j) {
      setflag[i][j] = 1;          // 0/1 = whether each i,j has been set
      cutsq[i][j] = cutoff_max;   // cutoff sq for each atom pair (neighbor.cpp thinks it is cutoff, NOT cutoff^2)
    }
  }

  return;
}

