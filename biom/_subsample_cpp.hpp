// -----------------------------------------------------------------------------
// Copyright (c) 2023-2023, The BIOM Format Development Team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
// -----------------------------------------------------------------------------

#ifndef _SUBSAMPLE_HPP
#define _SUBSAMPLE_HPP

#include <random>
#include <vector>

class WeightedSample
  {
  public:
   WeightedSample(uint32_t _max_count, uint32_t _n, uint32_t random_seed);
   void do_sample(double* data_arr, int start, int end);

private:
    uint32_t max_count;
    uint32_t n;
    std::mt19937 generator;

    // use persistent buffer to minimize allocation costs
    std::vector<uint64_t> data;  // original values
    std::vector<uint32_t> sample_out;     // random output buffer
    std::vector<uint32_t> data_out; // computed values
};

#endif
