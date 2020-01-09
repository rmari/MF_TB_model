#ifndef __MFTB__configuration__
#define __MFTB__configuration__

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include "MersenneTwister.h"
#include "global.h"

namespace TBConfig {

inline mfloat getSize(unsigned np, mfloat phi, mfloat rad) {
  return sqrt(np*M_PI*rad*rad/phi);
}

inline void printConf(std::ostream &out, std::vector<std::array<mfloat, 2>> conf) {
  out.precision(10);
  out << std::scientific;
  for (auto &pos: conf) {
    out << pos[0] << " " << pos[1] << std::endl;
  }
}

inline void printYapConf(std::ostream &out, std::vector<std::array<mfloat, 2>> conf) {
  out.precision(10);
  out << std::scientific;
  out << "r  1" << std::endl;
  out << "@  3" << std::endl;
  for (auto &pos: conf) {
    out << "c " << pos[0] << " 0 " << pos[1] << std::endl;
  }
  out << std::endl;
}

inline std::vector<std::array<mfloat, 2>> randConf(unsigned np, mfloat system_size, MTRand &r_gen) {
  std::vector<std::array<mfloat, 2>> conf(np);
  for (auto &pos: conf) {
    pos[0] = system_size*r_gen.rand();
    pos[1] = system_size*r_gen.rand();
  }
  return conf;
}


class Configuration {
public:
  Configuration(unsigned np, mfloat phi, mfloat rad):
  pos(np),
  _rad(rad),
  _phi(phi),
  _size(getSize(np, phi, rad)),
  half_size(_size/2) {};
  std::vector<std::array<mfloat, 2>> pos;
  mfloat size() const {return _size;};
  mfloat phi() const {return _phi;};
  mfloat dist_square(unsigned i, unsigned j) const {
    auto d0 = pos[j][0] - pos[i][0];
    while (d0 > half_size)
      d0 -= _size;
    while (d0 < -half_size)
      d0 += _size;
    auto d1 = pos[j][1] - pos[i][1];
    while (d1 > half_size)
      d1 -= _size;
    while (d1 < -half_size)
      d1 += _size;
    return d0*d0+d1*d1;
  }
  mfloat dist_square(std::array<mfloat, 2> pos_i, unsigned j) const {
    auto d0 = pos[j][0] - pos_i[0];
    if (d0 > half_size)
      d0 -= _size;
    if (d0 < -half_size)
      d0 += _size;
    auto d1 = pos[j][1] - pos_i[1];
    if (d1 > half_size)
      d1 -= _size;
    if (d1 < -half_size)
      d1 += _size;
    return d0*d0+d1*d1;
  }
  unsigned np() const {return pos.size();};
  mfloat rad() const {return _rad;};
private:
  mfloat _rad;
  mfloat _phi;
  mfloat _size;
  mfloat half_size;
};

inline Configuration read_conf(std::string fname) {
  std::ifstream fin (fname);
  if (!fin.good()) {
    std::cerr << " could not open " << fname << std::endl;
  }
  // auto startN = fname.find("_N")+2;
  auto startphi = fname.find("_phi");
  // unsigned N = std::stoi(fname.substr(startN, startphi-startN));

  startphi += 4;
  auto endphi = fname.find("_range");
  if (endphi == std::string::npos) {
    endphi = fname.find(".dat");
  }
  mfloat phi = std::stod(fname.substr(startphi, endphi-startphi));

  std::string dumb;
  unsigned t;
  fin >> dumb >> t;
  std::cout << "init conf time : " << t << std::endl;
  std::array<mfloat, 2> xy;
  std::vector<std::array<mfloat, 2>> positions;
  while (true) {
    fin >> xy[0] >> xy[1];
    if (fin.eof()) {
      break;
    }
    positions.push_back(xy);
  }
  mfloat rad = 1;
  Configuration conf(positions.size(), phi, rad);
  conf.pos = positions;
  return conf;
}

} // TBConfig namespace

#endif /* defined(__MFTB__configuration__) */
