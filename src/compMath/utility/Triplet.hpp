//
// Created by petrov on 12.02.2022.
//

#ifndef SLAE_TRIPLET_HPP
#define SLAE_TRIPLET_HPP
#include <cstdio>

template<typename T>
struct Triplet{
    std::size_t i;
    std::size_t j;
    T value;
    bool operator<(Triplet<T> const & rgh) const{
        return this->i < rgh.i || (this->i == rgh.i && this->j < rgh.j);
    };
};
#endif//SLAE_TRIPLET_HPP
