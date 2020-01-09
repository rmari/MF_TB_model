#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>


#include "global.h"
#include "PeriodicBoundary.h"
#include "boxing.h"
#include "configuration.h"


namespace TBDynamics {

std::set<unsigned> nonzero(const std::vector<bool> &v) {
  std::set<unsigned> nnz;
  for (unsigned i=0; i<v.size(); i++) {
    nnz.insert(i);
  }
  return nnz;
}

void findActive(const TBConfig::Configuration &conf,
                std::vector<bool> &to_be_moved,
                const TBBoxing::ExclusiveBoxSet &bxset,
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

void findPassiveDisp(const TBConfig::Configuration &conf,
		                const std::vector<bool> &to_be_moved,
		                std::vector<mfloat> &passive_disp,
		                MTRand &r_gen) {
	passive_disp.assign(passive_disp.size(), 0);
	for (unsigned i=0; i<conf.np(); i++) {
	  	if (to_be_moved[i]) {
		    auto pos_i = conf.pos[i];
			for (unsigned j=0; j<conf.np(); j++) {
				if (!to_be_moved[j]) {
			  		auto d2 = conf.dist_square(pos_i, j);
			  		passive_disp[j] += 2*(r_gen.rand()-0.5)/d2;
			  	}
			}
		}
	}
}

void findActiveFromMoved(const TBConfig::Configuration &conf,
                         std::set<unsigned> &just_moved,
                         const TBBoxing::InclusiveBoxSet &bxset,
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


unsigned moveParticles(TBConfig::Configuration &conf, const std::vector<bool> &to_be_moved, 
						const std::vector<mfloat> &passive_disp, 
						MTRand &r_gen, mfloat range) {
  unsigned active_nb = 0;
  for (unsigned i=0; i<conf.np(); i++) {
    if (to_be_moved[i]) {
      conf.pos[i][0] += r_gen.randNorm(0., range);
      conf.pos[i][1] += r_gen.randNorm(0., range);
      active_nb++;
    } else {
	    conf.pos[i][0] += passive_disp[i]*r_gen.randNorm(0., range);
    	conf.pos[i][1] += passive_disp[i]*r_gen.randNorm(0., range);
    }
  }
  return active_nb;
}

TBBoxing::ExclusiveBoxSet setupExclusiveBoxing (const TBConfig::Configuration &conf,
                                                const TBPBC::Periodizer &pbc) {
  double exclusive_bx_size;
  if (conf.np() > 10000000) {
    exclusive_bx_size = 0.01*sqrt(conf.np())*conf.rad();
  } else {
    exclusive_bx_size = 10*conf.rad();
  }

  return TBBoxing::ExclusiveBoxSet(exclusive_bx_size, conf.np(), conf.size(), pbc);
}


void runTB(TBConfig::Configuration &conf,
           MTRand &r_gen,
           std::string dfile_name,
           std::string cfile_name,
           mfloat range)
{
  TBPBC::Periodizer pbc(conf.size());

  auto bxset = setupExclusiveBoxing(conf, pbc);

  bxset.box(conf.pos);
  bxset.buildNeighborhoodContainers();

  std::vector<bool> to_be_moved(conf.np(), false);
  std::vector<mfloat> passive_disp(conf.np(), 0);
  unsigned active_nb;
  double active_prop;
  unsigned tcount = 0;
  unsigned out_data_period = 10;
  unsigned out_conf_period = 1000;
  unsigned simu_stop = 1000000;
  mfloat diam2 = 4*conf.rad()*conf.rad();


  checkFileExists(dfile_name);
  std::ofstream out_data (dfile_name.c_str());
  checkFileExists(cfile_name);
  std::ofstream out_conf (cfile_name.c_str());
  do {
	findActive(conf, to_be_moved, bxset, diam2);
	findPassiveDisp(conf, to_be_moved, passive_disp, r_gen);
	active_nb = moveParticles(conf, to_be_moved, passive_disp, r_gen, range);
	pbc.periodize(conf.pos);
    active_prop = active_nb;
    active_prop /= conf.np();

    if (tcount%out_data_period == 0) {
      out_data << tcount << " " << active_prop << std::endl;
      std::cout << tcount << " " << active_prop << std::endl;
    }
    if (tcount%out_conf_period == 0) {
      std::ofstream out_conf (cfile_name.c_str());
      out_conf << "time: " << tcount << std::endl;
      TBConfig::printConf(out_conf, conf.pos);
      out_conf.close();
    }

	bxset.box(conf.pos);
	bxset.buildNeighborhoodContainers();
    tcount++;
	// TBConfig::printYapConf(out_conf, conf.pos);
  } while(active_nb&&tcount<simu_stop);

  out_data << tcount << " " << active_prop << std::endl;
  out_data.close();
  std::cout << tcount << " " << active_prop << std::endl;

  TBConfig::printYapConf(out_conf, conf.pos);
  out_conf.close();
}

} // TBDynamics namespace
