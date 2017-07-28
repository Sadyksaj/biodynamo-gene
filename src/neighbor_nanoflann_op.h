#ifndef NEIGHBOR_NANOFLANN_OP_H_
#define NEIGHBOR_NANOFLANN_OP_H_

#include <utility>
#include <vector>
#include "inline_vector.h"
#include "nanoflann.h"

namespace bdm {

using nanoflann::KDTreeSingleIndexAdaptorParams;
using nanoflann::L2_Simple_Adaptor;
using nanoflann::KDTreeSingleIndexAdaptor;

// https://github.com/jlblancoc/nanoflann/blob/master/examples/pointcloud_adaptor_example.cpp
template <typename Derived>
struct NanoFlannAdapter {
  using coord_t = float;

  const Derived& obj;  //!< A const ref to the data set origin

  /// The constructor that sets the data set source
  explicit NanoFlannAdapter(const Derived& obj_) : obj(obj_) {}

  /// CRTP helper method
  inline const Derived& derived() const { return obj; }  // NOLINT

  /// Must return the number of data points
  inline size_t kdtree_get_point_count() const {  // NOLINT
    return derived().size();
  }

  /// Returns the distance between the vector "p1[0:size-1]" and the data point
  /// with index "idx_p2" stored in the class:
  inline coord_t kdtree_distance(const coord_t* p1,  // NOLINT
                                 const size_t idx_p2, size_t size) const {
    const auto& position = derived()[idx_p2].GetPosition();
    const coord_t d0 = p1[0] - position[0];
    const coord_t d1 = p1[1] - position[1];
    const coord_t d2 = p1[2] - position[2];
    return d0 * d0 + d1 * d1 + d2 * d2;
  }

  /// Returns the dim'th component of the idx'th point in the class:
  /// Since this is inlined and the "dim" argument is typically an immediate
  /// value, the "if/else's" are actually solved at compile time.
  inline coord_t kdtree_get_pt(const size_t idx, int dim) const {  // NOLINT
    return derived()[idx].GetPosition()[dim];
  }

  /// Optional bounding-box computation: return false to default to a standard
  /// bbox computation loop.
  ///   Return true if the BBOX was already computed by the class and returned
  ///   in
  ///   "bb" so it can be avoided to redo it again.
  ///   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3
  ///   for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const {  // NOLINT
    return false;
  }
};

class NeighborNanoflannOp {
 public:
  NeighborNanoflannOp() {}
  explicit NeighborNanoflannOp(float distance) : distance_(distance) {}
  ~NeighborNanoflannOp() {}

  template <typename TContainer>
  void Compute(TContainer* cells) const {
    typedef NanoFlannAdapter<TContainer> Adapter;
    const Adapter nf_cells(*cells);  // The adaptor

    // construct a 3D kd-tree index:
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, Adapter>,
                                     Adapter, 3>
        MyKdTree;

    // three dimensions; max leafs: 10
    MyKdTree index(3, nf_cells, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    std::vector<std::pair<size_t, float>> ret_matches;
    ret_matches.reserve(8);
    nanoflann::SearchParams params;
    params.sorted = false;
    float search_radius = distance_;
    InlineVector<int, 8> neighbors;

// calc neighbors
#pragma omp parallel for firstprivate(ret_matches, params, search_radius, \
                                      neighbors)
    for (size_t i = 0; i < cells->size(); i++) {
      // fixme make param
      // according to roman 50 - 100 micron
      ret_matches.clear();
      neighbors.clear();

      const auto& position = (*cells)[i].GetPosition();

      // calculate neighbors
      const size_t n_matches =
          index.radiusSearch(&position[0], search_radius, ret_matches, params);

      // transform result (change data structure - remove self from list)
      neighbors.reserve(n_matches - 1);
      for (size_t j = 0; j < n_matches; j++) {
        if (ret_matches[j].first != i) {
          neighbors.push_back(ret_matches[j].first);
        }
      }
      (*cells)[i].SetNeighbors(neighbors);
    }
  }

 private:
  float distance_ = 3000;
};

}  // namespace bdm

#endif  // NEIGHBOR_NANOFLANN_OP_H_
