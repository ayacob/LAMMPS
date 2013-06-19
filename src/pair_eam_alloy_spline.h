#ifdef PAIR_CLASS

PairStyle(eam/spline,PairEAMSpline)

#else

#ifndef LMP_PAIR_EAM_SPLINE_H
#define LMP_PAIR_EAM_SPLINE_H

#include "pair.h"
#include "spline.h"

#include <vector>
#include <map>
#include <iostream>

namespace LAMMPS_NS {

class PairEAMSpline : public Pair
{
public:
  PairEAMSpline(class LAMMPS *);
  virtual ~PairEAMSpline();

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
  virtual void read_file_alloy(std::string); // read binary alloy potential file

private:
  struct EAM2Body {
    int idx;
    double r, recip, nx, ny, nz, f, df;
  };

  std::vector<EAM2Body> two_body_pair_info;

  std::vector<std::string> map_atom_type;
  std::map<std::string, int> elements;

};  // PairEAMSpline

} // LAMMPS_NS

#endif  // LMP_PAIR_EAM_SPLINE_H
#endif  // PAIR_CLASS
