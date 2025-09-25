/* This file is part of COVISE.

   You can use it under the terms of the GNU Lesser General Public License
   version 2.1 or later, see lgpl-2.1.txt.

 * License: LGPL 2+ */

#include "Map1D.h"
#include <algorithm>
#include <stdexcept>

Map1D::Map1D(int l, const int *list): mapping_(nullptr), min_(0), max_(0), length_(0), method_(TRIVIAL)
{
    setMap(l, list);
}

Map1D::~Map1D()
{
    delete[] mapping_;
}

Map1D &Map1D::operator=(const Map1D &rhs)
{
    if (this != &rhs) {
        delete[] mapping_;

        min_ = rhs.min_;
        max_ = rhs.max_;
        length_ = rhs.length_;
        method_ = rhs.method_;
        labels_ = rhs.labels_;

        if (rhs.mapping_ && method_ == TRIVIAL) {
            int size = max_ - min_ + 1;
            mapping_ = new int[size];
            std::copy(rhs.mapping_, rhs.mapping_ + size, mapping_);
        } else {
            mapping_ = nullptr;
        }
    }
    return *this;
}

void Map1D::setMap(int l, const int *list)
{
    delete[] mapping_;
    mapping_ = nullptr;
    labels_.clear();

    if (l <= 0 || !list) {
        length_ = 0;
        method_ = TRIVIAL;
        return;
    }

    length_ = l;

    // Find min and max values
    min_ = *std::min_element(list, list + l);
    max_ = *std::max_element(list, list + l);

    int span = max_ - min_ + 1;

    if (span <= TRIVIAL_LIMIT) {
        // Use trivial mapping (array-based)
        method_ = TRIVIAL;
        mapping_ = new int[span];

        // Initialize with invalid indices
        std::fill(mapping_, mapping_ + span, -1);

        // Set up mapping
        for (int i = 0; i < l; ++i) {
            mapping_[list[i] - min_] = i;
        }
    } else {
        // Use hash mapping
        method_ = HASHING;
        for (int i = 0; i < l; ++i) {
            labels_[list[i]] = i;
        }
    }
}

const int &Map1D::operator[](int i) const
{
    static int invalid_index = -1;

    if (method_ == TRIVIAL) {
        if (mapping_ && i >= min_ && i <= max_) {
            return mapping_[i - min_];
        }
    } else { // HASHING
        auto it = labels_.find(i);
        if (it != labels_.end()) {
            return it->second;
        }
    }

    return invalid_index;
}
