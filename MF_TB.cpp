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

Configuration read_conf(std::string fname) {
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

std::tuple<unsigned, mfloat, mfloat, std::string>  parse_args(int argc, char **argv) {
  const struct option longopts[] = {
    {"np",   required_argument, 0, 'n'},
    {"phi",  required_argument, 0, 'p'},
    {"range",  required_argument, 0, 'r'},
    {"conf",  required_argument, 0, 'c'},
    {0, 0, 0, 0},
  };

  int index;
  int c;
  unsigned np = 0;
  mfloat phi = 1;
  mfloat range = -1;
  std::string conf = "";
  while ((c = getopt_long(argc, argv, "n:p:r:c:", longopts, &index)) != -1) {
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
      case 'c':
        conf = std::string(optarg);
        break;
      case '?':
        /* getopt already printed an error message. */
        break;
      default:
        abort();
    }
  }
  return std::make_tuple(np, phi, range, conf);
}

void findActive(const Configuration &conf,
                std::vector<bool> &to_be_moved,
                const ExclusiveBoxSet &bxset,
                mfloat diam2) {
  to_be_moved.assign(conf.np(), false);
  for (unsigned i=0; i<conf.np(); i++) {
    auto neigh_out = bxset.neighborsOffBox(i);
    auto pos_i = conf.pos[i];
    for (auto j: neigh_out) {
      if (conf.dist_square(pos_i, j) < diam2) {
        to_be_moved[i] = true;
        to_be_moved[j] = true;
      }
    }
    auto neigh_in = bxset.neighborsInBox(i);
    for (auto j: neigh_in) {
      if (i != j && conf.dist_square(pos_i, j) < diam2) {
        to_be_moved[i] = true;
        to_be_moved[j] = true;
      }
    }
  }
}

void findActiveFromMoved(const Configuration &conf,
                         std::set<unsigned> &just_moved,
                         const InclusiveBoxSet &bxset,
                         mfloat diam2) {
  std::set<unsigned> to_be_moved;
  for (auto i: just_moved) {
    auto neigh_out = bxset.neighbors(i);
    auto pos_i = conf.pos[i];
    for (auto j: neigh_out) {
      if (i != j && conf.dist_square(pos_i, j) < diam2) {
        to_be_moved.insert(i);
        to_be_moved.insert(j);
      }
    }
  }
  just_moved = to_be_moved;
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

unsigned moveParticlesMF(Configuration &conf, const std::set<unsigned> &to_be_moved, MTRand &r_gen) {
  for (auto i: to_be_moved) {
      conf.pos[i][0] = conf.size()*r_gen.rand();
      conf.pos[i][1] = conf.size()*r_gen.rand();
  }
  return to_be_moved.size();
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

unsigned moveParticles(Configuration &conf, const std::set<unsigned> &to_be_moved, MTRand &r_gen, mfloat range) {
  for (auto i: to_be_moved) {
    conf.pos[i][0] += r_gen.randNorm(0., range);
    conf.pos[i][1] += r_gen.randNorm(0., range);
  }
  return to_be_moved.size();
}

std::set<unsigned> nonzero(const std::vector<bool> &v) {
  std::set<unsigned> nnz;
  for (unsigned i=0; i<v.size(); i++) {
    nnz.insert(i);
  }
  return nnz;
}

int main(int argc, char **argv)
{
  MTRand r_gen;

  auto params = parse_args(argc, argv);
  mfloat range = std::get<2>(params);
  bool mean_field = range < 0;

  auto conffile = std::get<3>(params);
  mfloat rad = 1;
  Configuration conf(std::get<0>(params), std::get<1>(params), rad);

  if (conffile != "") {
    conf = read_conf(conffile);
    if (std::get<0>(params) != conf.np()) {
      std::cerr << " conf size inconsistant with CL args" << std::endl;
    }
    if (std::get<1>(params) != conf.phi()) {
      std::cerr << " conf phi inconsistant with CL args"<< std::endl;
    }
  } else {
    conf.pos = randConf(conf.np(), conf.size(), r_gen);
  }
  Periodizer pbc(conf.size());

  ExclusiveBoxSet bxset(10*rad, conf.np(), conf.size(), pbc);
  InclusiveBoxSet bxset_in(4*rad, conf.np(), conf.size(), pbc);

  bxset.box(conf.pos);
  bxset.buildNeighborhoodContainers();

  std::vector<bool> to_be_moved(conf.np(), false);
  unsigned active_nb;
  double active_prop;
  unsigned tcount = 0;
  unsigned out_data_period = 10;
  unsigned out_conf_period = 1000;
  unsigned simu_stop = 1000000;
  mfloat diam2 = 4*conf.rad()*conf.rad();

  std::string dfile_name, cfile_name;
  if (!conffile.empty()) {
    cfile_name = conffile;
    cfile_name.replace(cfile_name.find("conf_mftb"), cfile_name.find("conf_mftb")+9, "conf2_mftb");
    dfile_name = cfile_name;
    dfile_name.replace(cfile_name.find("conf2"), cfile_name.find("conf2")+5, "data2");
  } else {
    if (mean_field) {
      dfile_name = "data_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+".dat";
      cfile_name = "conf_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+".dat";
    } else {
      dfile_name = "data_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+"_range"+std::to_string(range)+".dat";
      cfile_name = "conf_mftb_N"+std::to_string(conf.np())+"_phi"+std::to_string(conf.phi())+"_range"+std::to_string(range)+".dat";
    }
  }
  checkFileExists(dfile_name);
  std::ofstream out_data (dfile_name.c_str());
  checkFileExists(cfile_name);
  std::ofstream out_conf (cfile_name.c_str());
  bool targeted_search = false;
  std::set<unsigned> to_be_moved_label;
  do {
    if (!targeted_search) {
      findActive(conf, to_be_moved, bxset, diam2);
      if (mean_field) {
        active_nb = moveParticlesMF(conf, to_be_moved, r_gen);
      } else {
        active_nb = moveParticles(conf, to_be_moved, r_gen, range);
        pbc.periodize(conf.pos);
      }
    } else {
      findActiveFromMoved(conf, to_be_moved_label, bxset_in, diam2);
      if (mean_field) {
        active_nb = moveParticlesMF(conf, to_be_moved_label, r_gen);
      } else {
        active_nb = moveParticles(conf, to_be_moved_label, r_gen, range);
        pbc.periodize(conf.pos);
      }
    }
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

    if (!targeted_search && conf.np() >= 1e5 && active_prop < 0.02) {
      to_be_moved_label = nonzero(to_be_moved);
      to_be_moved.clear();
      targeted_search = true;
    }
    if (targeted_search) {
      bxset_in.box(conf.pos);
    } else {
      bxset.box(conf.pos);
      bxset.buildNeighborhoodContainers();
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

  return 0;
}
