#include "spline.h"
#include "comm.h"
#include "error_spline.h"

#include "mpi.h"
#include <string>
#include <sstream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Spline::Spline() : nknots_(0), yp0_(0), ypn_(0), xmax_shifted_(0), step_(0), invstep_(0)
{
  // ctor
}

/* ---------------------------------------------------------------------- */

Spline::Spline(const Spline& spline)
{
  nknots_ = spline.nknots_;
  resize(); // resize x,y,ypp vectors

  yp0_ = spline.yp0_;
  ypn_ = spline.ypn_;

  for (int k=0; k<nknots_; ++k) {
    x_[k] = spline.x_[k];
    y_[k] = spline.y_[k];
    ypp_[k] = spline.ypp_[k];
  }

  rehash(); // setup xmaxshift, step, and invstep
}

/* ---------------------------------------------------------------------- */

Spline::~Spline()
{
  // dtor
}

/* ----------------------------------------------------------------------
   Get radial cutoff of spline (basically last knot)
------------------------------------------------------------------------- */

double Spline::get_cutoff() const
{
  return x_[nknots_-1];
}

/* ----------------------------------------------------------------------
   Broadcasts the spline function parameters to all processors
------------------------------------------------------------------------- */

void Spline::communicate(Comm* &comm, MPI_Comm& world)
{
  MPI_Bcast(&nknots_, 1, MPI_INT, 0, world);
  MPI_Bcast(&yp0_, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ypn_, 1, MPI_DOUBLE, 0, world);
  if(comm->me != 0)
    resize(); // resize vectors that make up spline using nknots_
  MPI_Bcast(&x_[0], nknots_, MPI_DOUBLE, 0, world);
  MPI_Bcast(&y_[0], nknots_, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ypp_[0], nknots_, MPI_DOUBLE, 0, world);
  if(comm->me != 0)
    rehash(); // setup spline
  return;
}

/* ----------------------------------------------------------------------
   Read in spline from input stream
------------------------------------------------------------------------- */

namespace LAMMPS_NS {

std::istream& operator>>(std::istream& is, Spline& spline)
{
  spline.read(is);
  return is;
}

} // LAMMPS_NS

/* ----------------------------------------------------------------------
   Assignment operator for splines
------------------------------------------------------------------------- */

Spline& Spline::operator=(const Spline& rhs)
{
  nknots_ = rhs.nknots_;
  resize(); // resize x,y,ypp vectors

  yp0_ = rhs.yp0_;
  ypn_ = rhs.ypn_;

  for (int k=0; k<nknots_; ++k) {
    x_[k] = rhs.x_[k];
    y_[k] = rhs.y_[k];
    ypp_[k] = rhs.ypp_[k];
  }

  rehash(); // setup xmaxshift, step, and invstep

  return *this;
}

/* ----------------------------------------------------------------------
   Read in spline data from input stream
------------------------------------------------------------------------- */

void Spline::read(std::istream& is)
{
  // 1st line lists # knots
  is >> nknots_;
  if (nknots_ < 2)
    throw ErrorSpline("Invalid number of spline knots in potential file");

  // Resize spline using nknots_
  resize();

  // 2nd line lists first and 2nd derivs
  is >> yp0_ >> ypn_;

  // 3rd line is garbage line
  std::string tmp_line;
  std::getline(is, tmp_line);
  std::getline(is, tmp_line);

  // Read in knots
  for (int k=0; k<nknots_; ++k) {
    std::string line;
    std::getline(is, line);

    std::stringstream line_str(line);
    line_str >> x_[k] >> y_[k] >> ypp_[k];

    if (line_str.fail())
      throw ErrorSpline("Invalid knot line in potential file");
  }

  // Setup spline
  rehash();

  return;
}

/* ----------------------------------------------------------------------
   Resize spline for # knots
------------------------------------------------------------------------- */

void Spline::resize()
{
  x_.resize(nknots_);
  y_.resize(nknots_);
  ypp_.resize(nknots_);
  return;
}

/* ----------------------------------------------------------------------
   Setup spline using limited data already possessed
------------------------------------------------------------------------- */

void Spline::rehash()
{
  xmax_shifted_ = x_[nknots_-1] - x_[0];
  step_ = x_[1] - x_[0];
  invstep_ = 1.0/step_;
  return;
}
