#include <iostream>
#include <fstream>
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
  mfloat size() {return _size;};
  mfloat phi() {return _phi;};
  mfloat dist_square(unsigned i, unsigned j) {
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
  unsigned np() {return pos.size();};
  mfloat rad() {return _rad;};
private:
  mfloat _rad;
  mfloat _phi;
  mfloat _size;
  mfloat half_size;
};

Configuration init_conf(int argc, char **argv) {
  const struct option longopts[] = {
		{"np",   required_argument, 0, 'n'},
		{"phi",  required_argument, 0, 'p'},
		{0, 0, 0, 0},
	};

	int index;
	int c;
  unsigned np = 0;
  mfloat phi = 1;
	while ((c = getopt_long(argc, argv, "n:p:", longopts, &index)) != -1) {
		switch (c) {
			case 'n':
        np = atoi(optarg);
				break;
			case 'p':
        phi = atof(optarg);
				break;
			case '?':
				/* getopt already printed an error message. */
				break;
			default:
				abort();
		}
	}


  mfloat rad = 1;
  Configuration conf(np, phi, rad);

  return conf;
}

int main(int argc, char **argv)
{

  auto conf = init_conf(argc, argv);
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
  mfloat diam2 = 4*conf.rad()*conf.rad();

  std::string dfile_name = "data_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+".dat";
  checkFileExists(dfile_name);
  std::ofstream out_data (dfile_name.c_str());
  std::string cfile_name = "conf_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+".dat";
  checkFileExists(cfile_name);
  std::ofstream out_conf (cfile_name.c_str());

  do {
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
    active_nb = 0;
    for (unsigned i=0; i<conf.np(); i++) {
      if (to_be_moved[i]) {
        conf.pos[i][0] = conf.size()*r_gen.rand();
        conf.pos[i][1] = conf.size()*r_gen.rand();
        active_nb++;
      }
    }
    pbc.periodize(conf.pos);
    bxset.box(conf.pos);
    bxset.buildNeighborhoodContainers();
    active_prop = active_nb;
    active_prop /= conf.np();
    if (tcount%100 == 0) {
      out_data << tcount << " " << active_prop << std::endl;
      std::cout << tcount << " " << active_prop << std::endl;
    }
    if (tcount%1000 == 0) {
      out_conf.seekp(0, out_conf.beg);
      out_conf << "time: " << tcount << std::endl;
      printConf(out_conf, conf.pos);
    }
    tcount++;
  } while(active_nb);

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
