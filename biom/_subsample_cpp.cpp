// -----------------------------------------------------------------------------
// Copyright (c) 2023-2023, The BIOM Format Development Team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
// -----------------------------------------------------------------------------

#include "_subsample_cpp.hpp"
#include <algorithm>

// Adapted from unifrac-binaries code:
// https://github.com/biocore/unifrac-binaries/blob/ba11b1b80c56ae13dff6b6e364352abb2f2b0faa/src/biom_subsampled.cpp#L115

// Equivalent to iterator over np.repeat
// https://github.com/biocore/biom-format/blob/b0e71a00ecb349a6f5f1ca64a23d71f380ddc19c/biom/_subsample.pyx#LL64C24-L64C55
class WeightedSampleIterator
  {
  public:
    // While we do not implememnt the whole random_access_iterator interface
    // we want the implementations to use operator- and that requires random
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = int64_t;
    using value_type        = uint32_t;
    using pointer           = const uint32_t*;
    using reference         = const uint32_t&;

    WeightedSampleIterator(uint64_t *_data_in, uint32_t _idx, uint64_t _cnt)
    : data_in(_data_in)
    , idx(_idx)
    , cnt(_cnt)
    {}

    reference operator*() const { return idx; }
    pointer operator->() const { return &idx; }

    WeightedSampleIterator& operator++()
    {  
       cnt++;
       if (cnt>=data_in[idx]) {
         cnt = 0;
         idx++;
       }
       return *this;
    }

    WeightedSampleIterator operator++(int) { WeightedSampleIterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== (const WeightedSampleIterator& a, const WeightedSampleIterator& b)
    {
       return (a.data_in == b.data_in) && (a.idx == b.idx) && (a.cnt==b.cnt);
    };

    friend bool operator!= (const WeightedSampleIterator& a, const WeightedSampleIterator& b)
    {
       return !((a.data_in == b.data_in) && (a.idx == b.idx) && (a.cnt==b.cnt));
    };

    friend int64_t operator-(const WeightedSampleIterator& b, const WeightedSampleIterator& a)
    {
       int64_t diff = 0;
       //assert(a.data_in == b.data_in);
       //assert(a.idx <= b.idx);
       //assert((a.idx > b.idx) || (a.cnt<=b.cnt));

       for (uint32_t i = a.idx; i<b.idx; i++) {
          diff += a.data_in[i];
       }

       return diff + b.cnt - a.cnt;
    };

  private:

    uint64_t *data_in;
    uint32_t idx; // index of data_in
    uint64_t cnt; // how deep in data_in[idx] are we (must be < data_in[idx])
};


WeightedSample::WeightedSample(uint32_t _max_count, uint32_t _n, uint32_t random_seed)
    : max_count(_max_count)
    , n(_n)
    , generator(random_seed)
    , data(max_count)
    , sample_out(n)
    , data_out(max_count)
{}

void WeightedSample::do_sample(double* data_base, int start, int end) {
        double* data_arr = data_base+start;
        unsigned int length = end-start;
        for (unsigned int j=0; j<length; j++) data_out[j] = 0;

        // note: We are assuming length>=n
        //      Enforced by the caller (via filtering)
        for (uint32_t j=0; j<length; j++) data[j] = data_arr[j];
        std::sample(WeightedSampleIterator(data.data(),0,0),
                    WeightedSampleIterator(data.data(),length,0),
                    sample_out.begin(), n,
                    generator);

        for (uint32_t j=0; j<n; j++) data_out[sample_out[j]]++;

        for (unsigned int j=0; j<length; j++) data_arr[j] = data_out[j];
}

