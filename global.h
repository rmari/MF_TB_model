#ifndef __MFTB__global__
#define __MFTB__global__

typedef double mfloat;

inline void checkFileExists(std::string fname) {
  std::ifstream file_test(fname.c_str());
  if (file_test.good()) {
  	file_test.close();
  	std::cerr << "The file '" << fname << "' already exists." << std::endl;
  	exit(1);
  } else {
  	file_test.close();
  }
}
#endif
