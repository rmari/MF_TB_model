#ifndef __MFTB__PBC__
#define __MFTB__PBC__

#include <vector>
#include <array>
#include "global.h"

namespace TBPBC {

class Periodizer{
public:
	Periodizer(mfloat system_size): size(system_size) {};
	void periodize(std::array<mfloat,2> &pos);
	std::array<mfloat,2> periodized(const std::array<mfloat,2> &pos);
	void periodize(std::vector<std::array<mfloat,2>> &pos);
	std::vector<std::array<mfloat,2>> periodized(const std::vector<std::array<mfloat,2>> &pos);
private:
	mfloat size;
};

inline void Periodizer::periodize(std::array<mfloat,2> &pos) {
	for (unsigned i=0; i<2; i++) {
		while (pos[i] > size)
			pos[i] -= size;
		while (pos[i] < 0)
			pos[i] += size;
	}
}

inline std::array<mfloat,2> Periodizer::periodized(const std::array<mfloat,2> &pos) {
	auto p = pos;
	periodize(p);
	return p;
}

inline void Periodizer::periodize(std::vector<std::array<mfloat,2>> &pos) {
	for (auto &p: pos) {
		for (unsigned i=0; i<2; i++) {
			while (p[i] > size)
				p[i] -= size;
			while (p[i] < 0)
				p[i] += size;
		}
	}
}

inline std::vector<std::array<mfloat,2>> Periodizer::periodized(const std::vector<std::array<mfloat,2>> &pos) {
	auto p = pos;
	periodize(p);
	return p;
}
} // namespace TBPBC

#endif /* defined(__MFTB__PBC__) */
