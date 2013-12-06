#ifdef PAIR_CLASS

PairStyle(eam/alloy/spline,PairEAMAlloySpline)

#else

#ifndef LMP_PAIR_EAM_ALLOY_SPLINE_H
#define LMP_PAIR_EAM_ALLOY_SPLINE_H

#include "pair.h"
#include "spline.h"

#include <vector>
#include <map>
#include <iostream>

namespace LAMMPS_NS {

class PairEAMAlloySpline : public Pair
{
public:
  PairEAMAlloySpline(class LAMMPS *);
  virtual ~PairEAMAlloySpline();

  virtual void compute(int, int);
  void allocate();
  void settings(int, char **);
  void coeff(int, char **);

  void init_style();
  double init_one(int, int);

protected:

  std::vector<Spline> phi; 
  std::vector<Spline> rho;
  std::vector<Spline> u;

  double cutoff_max;

  virtual void read_file(std::string);    // Read in potential from file

private:
  struct EAM2Body {
    int idx;
    double r, recip, nx, ny, nz, f, df;
  };

  std::vector<EAM2Body> two_body_pair_info;

  std::vector<std::string> map_atom_type;
  std::map<std::string, int> elements;

};  // PairEAMAlloySpline

} // LAMMPS_NS

#endif  // LMP_PAIR_EAM_ALLOY_SPLINE_H
#endif  // PAIR_CLASS
