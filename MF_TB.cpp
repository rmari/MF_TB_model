#include <iostream>
#include <fstream>
#include <tuple>

#include <getopt.h>

#include "MersenneTwister.h"

typedef double mfloat;

#include "boxing.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


void printConf(std::ostream &out, std::vector<std::array<mfloat, 2>> conf) {
  for (auto &pos: conf) {
    out << pos[0] << " " << pos[1] << std::endl;
  }
}

std::vector<std::array<mfloat, 2>> randConf(unsigned np, mfloat system_size, MTRand &r_gen) {
  std::vector<std::array<mfloat, 2>> conf(np);
  for (auto &pos: conf) {
    pos[0] = system_size*r_gen.rand();
    pos[1] = system_size*r_gen.rand();
  }
  return conf;
}

mfloat getSize(unsigned np, mfloat phi, mfloat rad) {
  return sqrt(np*M_PI*rad*rad/phi);
}

void checkFileExists(std::string fname) {
  std::ifstream file_test(fname.c_str());
  if (file_test.good()) {
  	file_test.close();
  	std::cerr << "The file '" << fname << "' already exists." << std::endl;
  	exit(1);
  } else {
  	file_test.close();
  }
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
    mfloat d2 = 0;
    for (auto d: {pos[j][0] - pos[i][0], pos[j][1] - pos[i][1]}) {
      while (d > half_size)
        d -= _size;
      while (d < -half_size)
        d += _size;
      d2 += d*d;
    }
    return d2;
  }
  unsigned np() const {return pos.size();};
  mfloat rad() const {return _rad;};
private:
  mfloat _rad;
  mfloat _phi;
  mfloat _size;
  mfloat half_size;
};

std::tuple<unsigned, mfloat, mfloat>  parse_args(int argc, char **argv) {
  const struct option longopts[] = {
    {"np",   required_argument, 0, 'n'},
    {"phi",  required_argument, 0, 'p'},
    {"range",  required_argument, 0, 'r'},
    {0, 0, 0, 0},
  };

  int index;
  int c;
  unsigned np = 0;
  mfloat phi = 1;
  mfloat range = -1;
  while ((c = getopt_long(argc, argv, "n:p:r:", longopts, &index)) != -1) {
    switch (c) {
      case 'n':
        np = atoi(optarg);
        break;
      case 'p':
        phi = atof(optarg);
        break;
      case 'r':
        range = atof(optarg);
        break;
      case '?':
        /* getopt already printed an error message. */
        break;
      default:
        abort();
    }
  }
  return std::make_tuple(np, phi, range);
}

void findActive(const Configuration &conf, std::vector<bool> &to_be_moved, const BoxSet &bxset, mfloat diam2) {
  to_be_moved.assign(conf.np(), false);
  for (unsigned i=0; i<conf.np(); i++) {
    auto neigh_out = bxset.neighborsOffBox(i);
    for (auto j: neigh_out) {
      if (conf.dist_square(i, j) < diam2) {
        to_be_moved[i] = true;
        to_be_moved[j] = true;
      }
    }
    auto neigh_in = bxset.neighborsInBox(i);
    for (auto j: neigh_in) {
      if (conf.dist_square(i, j) < diam2) {
        to_be_moved[i] = true;
        to_be_moved[j] = true;
      }
    }
  }
}

unsigned moveParticlesMF(Configuration &conf, const std::vector<bool> &to_be_moved, MTRand &r_gen) {
  unsigned active_nb = 0;
  for (unsigned i=0; i<conf.np(); i++) {
    if (to_be_moved[i]) {
      conf.pos[i][0] = conf.size()*r_gen.rand();
      conf.pos[i][1] = conf.size()*r_gen.rand();
      active_nb++;
    }
  }
  return active_nb;
}


unsigned moveParticles(Configuration &conf, const std::vector<bool> &to_be_moved, MTRand &r_gen, mfloat range) {
  unsigned active_nb = 0;
  for (unsigned i=0; i<conf.np(); i++) {
    if (to_be_moved[i]) {
      conf.pos[i][0] += r_gen.randNorm(0., range);
      conf.pos[i][1] += r_gen.randNorm(0., range);
      active_nb++;
    }
  }
  return active_nb;
}

int main(int argc, char **argv)
{

  auto params = parse_args(argc, argv);
  mfloat range = std::get<2>(params);
  bool mean_field = range < 0;
  mfloat rad = 1;
  Configuration conf(std::get<0>(params), std::get<1>(params), rad);
  MTRand r_gen;
  conf.pos = randConf(conf.np(), conf.size(), r_gen);

  Periodizer pbc(conf.size());

  BoxSet bxset(2, conf.np(), conf.size(), pbc);
  bxset.box(conf.pos);
  bxset.buildNeighborhoodContainers();

  std::vector<bool> to_be_moved(conf.np(), false);
  unsigned active_nb;
  double active_prop;
  unsigned tcount = 0;
  unsigned out_data_period = 10;
  unsigned out_conf_period = 1000;
  unsigned simu_stop = 50000;
  mfloat diam2 = 4*conf.rad()*conf.rad();

  std::string dfile_name, cfile_name;
  if(mean_field) {
    dfile_name = "data_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+".dat";
    cfile_name = "conf_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+".dat";
  } else {
    dfile_name = "data_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+"_range"+std::to_string(range)+".dat";
    cfile_name = "conf_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+"_range"+std::to_string(range)+".dat";
  }
  checkFileExists(dfile_name);
  std::ofstream out_data (dfile_name.c_str());
  checkFileExists(cfile_name);
  std::ofstream out_conf (cfile_name.c_str());

  do {
    findActive(conf, to_be_moved, bxset, diam2);
    if (mean_field) {
      active_nb = moveParticlesMF(conf, to_be_moved, r_gen);
    } else {
      active_nb = moveParticles(conf, to_be_moved, r_gen, range);
      pbc.periodize(conf.pos);
    }

    bxset.box(conf.pos);
    bxset.buildNeighborhoodContainers();
    active_prop = active_nb;
    active_prop /= conf.np();
    if (tcount%out_data_period == 0) {
      out_data << tcount << " " << active_prop << std::endl;
      std::cout << tcount << " " << active_prop << std::endl;
    }
    if (tcount%out_conf_period == 0) {
      out_conf.seekp(0, out_conf.beg);
      out_conf << "time: " << tcount << std::endl;
      printConf(out_conf, conf.pos);
    }
    tcount++;
  } while(active_nb&&tcount<simu_stop);

  out_data << tcount << " " << active_prop << std::endl;
  std::cout << tcount << " " << active_prop << std::endl;
  out_conf.seekp(0, out_conf.beg);
  out_conf << "time: " << tcount << std::endl;
  printConf(out_conf, conf.pos);

  out_conf.close();
  out_data.close();

  // for (unsigned i=0; i<conf.np(); i++) {
  //   auto neigh_out = bxset.neighborsOffBox(i);
  //   for (auto j: neigh_out) {
  //     std::cout << conf.dist_square(i, j) << std::endl;
  //   }
  //   auto neigh_in = bxset.neighborsInBox(i);
  //   for (auto j: neigh_in) {
  //     std::cout << conf.dist_square(i, j) << std::endl;
  //   }
  // }

}
