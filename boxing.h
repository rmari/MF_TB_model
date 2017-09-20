//
//  Box.h
//  initially from LF_DEM (Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.)
//
//
//

#ifndef __MFTB__Box__
#define __MFTB__Box__
#include <set>
#include <unordered_set>
#include <vector>
#include <array>
#include <sstream>
#include <algorithm>

#include <stdexcept>

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

class Box{
private:
	std::vector <Box*> neighbors;
	std::set <unsigned> container; // ordered, important for no double counting
	std::vector <unsigned>	neighborhood_container;

public:
	void addNeighborBox(Box* neigh_box);
	void resetNeighborBoxes();
	void add(unsigned);
	void remove(unsigned);
	void buildNeighborhoodContainer();

	const std::vector <unsigned> & neighborsOffBox() const;
	std::vector <unsigned> neighborsInBox(unsigned i) const;
	const std::set <unsigned> & getContainer() const {return container;}
};



inline void Box::addNeighborBox(Box* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	neighbors.push_back(neigh_box);
}

inline void Box::resetNeighborBoxes() {
	neighbors.clear();
	neighborhood_container.clear();
}

inline void Box::add(unsigned i)
{
	container.insert(i);
}

inline void Box::remove(unsigned i)
{
	container.erase(i);
}

inline void Box::buildNeighborhoodContainer()
{
	neighborhood_container.clear();
	std::size_t size = 0;
	for (const auto& box : neighbors) {
		size += box->getContainer().size();
	}
	neighborhood_container.resize(size);
	int j = 0;
	for (const auto& box : neighbors) {
		for (const int& k : box->container) {
			neighborhood_container[j] = k;
			j++;
		}
	}
}

const std::vector <unsigned> & Box::neighborsOffBox() const {
	return neighborhood_container;
}

std::vector <unsigned> Box::neighborsInBox(unsigned i) const {
	return std::vector <unsigned>(++container.find(i), container.end());
}

class ExtendedBox{
private:
	std::vector <ExtendedBox*> neighbors;
	std::vector <unsigned> container;

public:
	void addNeighborBox(ExtendedBox* neigh_box);
	void addFromNeighbor(unsigned);
	void add(unsigned);
	void removeFromNeighbor(unsigned);
	void remove(unsigned);

	std::vector <unsigned> getContainer() const {return
	  																						std::vector<unsigned>(container.begin(), container.end());}

};

inline void ExtendedBox::addNeighborBox(ExtendedBox* neigh_box)
{
	if (neigh_box == this) {
		return;
	}
	neighbors.push_back(neigh_box);
}

inline void ExtendedBox::add(unsigned i)
{
	container.push_back(i);
	for (auto &bx: neighbors) {
		bx->addFromNeighbor(i);
	}
}

inline void ExtendedBox::addFromNeighbor(unsigned i)
{
	container.push_back(i);
}

inline void ExtendedBox::remove(unsigned i)
{
	removeFromNeighbor(i);
	for (auto &bx: neighbors) {
		bx->removeFromNeighbor(i);
	}
}

inline void ExtendedBox::removeFromNeighbor(unsigned i)
{
	auto elm = std::find(container.begin(), container.end(), i);
	*elm = *(--container.end());
	container.pop_back();
}

template<class BxClass>
class BoxSet{
protected:
	unsigned box_nb;
	bool _is_boxed;
	mfloat box_size;
	std::vector<BxClass> boxes;
	std::vector<BxClass*> boxMap;
	Periodizer pbc;
	mfloat total_size;
	BxClass* whichBox(const std::array<mfloat,2> &pos);
	virtual void assignNeighbors()=0;
public:
	BoxSet(mfloat box_min_size,
				std::size_t np,
				mfloat system_size,
				Periodizer PBC);
	void box(unsigned i, std::array<mfloat, 2> position_i);
	void box(std::vector<std::array<mfloat, 2>> position);
};

template<class BxClass>
inline BoxSet<BxClass>::BoxSet(mfloat box_min_size,
                      std::size_t np,
											mfloat system_size,
											Periodizer PBC):
_is_boxed(false),
pbc(PBC),
total_size(system_size)
{
	std::cout << "Setting up Cell List System ... ";
	if (box_min_size > system_size) {
		throw std::runtime_error(" system smaller than boxes ");
	}
	boxMap.resize(np);
	for (unsigned i=0; i<boxMap.size(); i++) {
		boxMap[i] = NULL;
	}
	box_nb = (unsigned)(system_size/box_min_size);
	box_size = system_size/box_nb;

	if (box_nb > 3) {
		_is_boxed = true;
	} else {
		box_nb = 1;
		box_size = system_size;
	}
	boxes.resize(box_nb*box_nb);
	std::cout << " [ok]" << std::endl;
}

template<class BxClass>
inline BxClass* BoxSet<BxClass>::whichBox(const std::array<mfloat, 2> &pos)
{
	unsigned label = (unsigned)(box_nb*(pos[0]/total_size))+box_nb*(unsigned)(box_nb*(pos[1]/total_size));
	if (label >= boxes.size()) {
		std::ostringstream error_str;
		error_str  << " BoxSet: trying to box position out of boundaries \"" << pos[0] << ", " << pos[1]	<< "\""\
		           << " (system size "<< box_nb*box_size << ")" << std::endl;
		throw std::runtime_error(error_str.str());
	}
	return &boxes[label];
}

template<class BxClass>
inline void BoxSet<BxClass>::box(unsigned i, std::array<mfloat, 2> position_i)
{
	BxClass* b = whichBox(position_i);
	if (b != boxMap[i]) {
		if (boxMap[i] != NULL) {
			boxMap[i]->remove(i);
		}
		b->add(i);
		boxMap[i] = b;
	}
}

template<class BxClass>
inline void BoxSet<BxClass>::box(std::vector<std::array<mfloat, 2>> position)
{
	for (unsigned i=0; i<position.size(); i++) {
		BxClass* b = whichBox(position[i]);
		if (b != boxMap[i]) {
			if (boxMap[i] != NULL) {
				boxMap[i]->remove(i);
			}
			b->add(i);
			boxMap[i] = b;
		}
	}
}

class ExclusiveBoxSet : public BoxSet<Box> {
private:
	void assignNeighbors();
public:
	ExclusiveBoxSet(mfloat box_min_size,
									std::size_t np,
									mfloat system_size,
									Periodizer PBC)
									: BoxSet(box_min_size, np, system_size, PBC) {assignNeighbors();};
	void buildNeighborhoodContainers();
	const std::vector <unsigned>& neighborsOffBox(unsigned i) const;
	std::vector <unsigned> neighborsInBox(unsigned i) const;
};


inline void ExclusiveBoxSet::assignNeighbors()
{
	for (unsigned i=0; i<box_nb; i++) {
		for (unsigned j=0; j<box_nb; j++) {
			std::array<mfloat, 2> pos = {(i+0.5)*box_size, (j+0.5)*box_size};
			unsigned box_label = i + j*box_nb;
			auto &bx = boxes[box_label];
			// rhs 3 neighbors
			auto delta0 = box_size;
			for (const auto& b : {-1, 0, 1}) {
				auto delta1 = b*box_size;
				auto pos_neigh = pos;
				pos_neigh[0] += delta0;
				pos_neigh[1] += delta1;
				bx.addNeighborBox(whichBox(pbc.periodized(pos_neigh)));
			}
			// right-on-top neighbor
			auto delta1 = box_size;
			auto pos_neigh = pos;
			pos_neigh[1] += delta1;
			bx.addNeighborBox(whichBox(pbc.periodized(pos_neigh)));
		}
	}
}

inline void ExclusiveBoxSet::buildNeighborhoodContainers()
{
	for (auto& bx : boxes) {
		bx.buildNeighborhoodContainer();
	}
}

inline const std::vector<unsigned>& ExclusiveBoxSet::neighborsOffBox(unsigned i) const
{
	return (boxMap[i])->neighborsOffBox();
}

inline std::vector<unsigned> ExclusiveBoxSet::neighborsInBox(unsigned i) const
{
	return (boxMap[i])->neighborsInBox(i);
}


class InclusiveBoxSet: public BoxSet<ExtendedBox> {
private:
	void assignNeighbors();
public:
	InclusiveBoxSet(mfloat box_min_size,
									std::size_t np,
									mfloat system_size,
									Periodizer PBC)
									: BoxSet(box_min_size, np, system_size, PBC) {assignNeighbors();};
	std::vector<unsigned> neighbors(unsigned i) const;

};

inline void InclusiveBoxSet::assignNeighbors()
{
	// becomes interesting if you don't check for all particles' neighborhoods,
	// but instead wants to know the neighborhood of a small subset
	for (unsigned i=0; i<box_nb; i++) {
		for (unsigned j=0; j<box_nb; j++) {
			std::array<mfloat, 2> pos = {(i+0.5)*box_size, (j+0.5)*box_size};
			unsigned box_label = i + j*box_nb;
			auto &bx = boxes[box_label];
			// rhs 3 neighbors
			auto delta0 = box_size;
			for (const auto& b : {-1, 0, 1}) {
				auto delta1 = b*box_size;
				auto pos_neigh = pos;
				pos_neigh[0] += delta0;
				pos_neigh[1] += delta1;
				bx.addNeighborBox(whichBox(pbc.periodized(pos_neigh)));
			}
			delta0 = -box_size;
			for (const auto& b : {-1, 0, 1}) {
				auto delta1 = b*box_size;
				auto pos_neigh = pos;
				pos_neigh[0] += delta0;
				pos_neigh[1] += delta1;
				bx.addNeighborBox(whichBox(pbc.periodized(pos_neigh)));
			}
			// right-on-top neighbor
			auto pos_neigh = pos;
			pos_neigh[1] += box_size;
			bx.addNeighborBox(whichBox(pbc.periodized(pos_neigh)));

			pos_neigh = pos;
			pos_neigh[1] -= box_size;
			bx.addNeighborBox(whichBox(pbc.periodized(pos_neigh)));
		}
	}
}

inline std::vector<unsigned> InclusiveBoxSet::neighbors(unsigned i) const
{
	return (boxMap[i])->getContainer();
}
#endif /* defined(__MFTB__Box__) */
