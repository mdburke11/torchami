#include "tami_base.hpp"

// Small functions compiled in c++ to see if there's a speed improvement

at::Tensor epsilon_2D(at::Tensor K) {
  return -2 *
         at::cos(K).sum(/*dim=*/1,
                        /*keepdim=*/true); // sum along dim=1, and keepdim=True
}