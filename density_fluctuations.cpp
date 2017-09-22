#include <iostream>
#include <fstream>
#include <string>

#include "configuration.h"
#include "boxing.h"
#include "MersenneTwister.h"


std::vector<double> getLs(double system_size, unsigned L_nb) {
  double L_step = pow(0.5*system_size, 1./L_nb);
  std::vector<double> Ls (L_nb, 1);


  for (int i = 1; i<L_nb; i++) {
    Ls[i] = Ls[i-1]*L_step;
  }
  return Ls;
}

int main(int argc, char **argv)
{
  std::string conffile(argv[1]);
  unsigned loc_nb = std::stoi(argv[2]);
  unsigned L_nb = std::stoi(argv[3]);


  auto conf = TBConfig::read_conf(conffile);
  TBPBC::Periodizer pbc(conf.size());

  std::string ofname = conffile;
  ofname.replace(ofname.find("conf"), ofname.find("conf")+4, "dens_fluc");
  std::ofstream outfile(ofname);

  MTRand r_gen;

  double rho_total = conf.np()/(conf.size()*conf.size());
  auto Ls = getLs(conf.size(), L_nb);

  for (auto L: Ls) {
    double L2 = L*L;
    double bxsize;
    if (conf.np() < 100000) {
      bxsize = L;
    } else {
      double min_size = 0.02*sqrt(conf.np());
      bxsize = L > min_size ? L : min_size;
    }
    TBBoxing::InclusiveBoxSet bxset(bxsize, conf.np(), conf.size(), pbc);
    bxset.box(conf.pos);

    double delta_rho_2 = 0;
    for (int _k = 0; _k<loc_nb; _k++) {
      double rho = 0;
      std::array<double, 2> pos = {conf.size()*r_gen.rand(), conf.size()*r_gen.rand()};
      auto neighbors = bxset.neighbors(pos);
      for (auto j: neighbors) {
        if (conf.dist_square(pos, j) < L2) {
          rho += 1;
        }
      }
      rho /= M_PI*L2;
      delta_rho_2 += (rho - rho_total)*(rho - rho_total);
    }
    delta_rho_2 /= loc_nb;
    std::cout << L << " " << L/conf.size() << " " << delta_rho_2 << std::endl;
    outfile << L << " " << L/conf.size() << " " << delta_rho_2 << std::endl;
  }
  outfile.close();
  return 0;
}
