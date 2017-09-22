#include <iostream>
#include <fstream>
#include <tuple>

#include <getopt.h>

#include "MersenneTwister.h"
#include "global.h"
#include "configuration.h"


#define no_argument 0
#define required_argument 1
#define optional_argument 2

namespace TBDynamics {
extern void runTB(TBConfig::Configuration &conf,
           MTRand &r_gen,
           std::string dfile_name,
           std::string cfile_name,
           mfloat range);
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




int main(int argc, char **argv)
{
  MTRand r_gen;

  auto params = parse_args(argc, argv);
  mfloat range = std::get<2>(params);
  bool mean_field = range < 0;

  auto conffile = std::get<3>(params);
  mfloat rad = 1;
  TBConfig::Configuration conf(std::get<0>(params), std::get<1>(params), rad);

  if (!conffile.empty()) {
    conf = TBConfig::read_conf(conffile);
    if (std::get<0>(params) != conf.np()) {
      std::cerr << " conf size inconsistant with CL args" << std::endl;
    }
    if (std::get<1>(params) != conf.phi()) {
      std::cerr << " conf phi inconsistant with CL args"<< std::endl;
    }
  } else {
    conf.pos = TBConfig::randConf(conf.np(), conf.size(), r_gen);
  }
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

  TBDynamics::runTB(conf, r_gen, dfile_name, cfile_name, range);

  return 0;
}
