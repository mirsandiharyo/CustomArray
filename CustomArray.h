/*
 * CustomArray.h
 *
 *  Created on: Jun 11, 2020
 *      Author: mirsandiharyo
 */

#ifndef CUSTOMARRAY_H_
#define CUSTOMARRAY_H_

#include <stdint.h>
#include <cmath>
#include <vector>
#include <string.h>
#include "AllignedAllocator.h"

namespace csar {

typedef uint32_t Index_t;

struct Index3 {
    Index_t i;
    Index_t j;
    Index_t k;

    inline friend std::ostream& operator<<(std::ostream &os, const Index3 &other) {
        os << "(" << other.i << "," << other.j << "," << other.k << ")";
        return os;
    }

    inline void SetZero() {
        i = j = k = 0;
    }
};

template<class T, std::size_t AlignmentN = 32>
class variables {
public:
    typedef std::vector<T, AlignmentAllocator<T, AlignmentN> > InternalVec;
    typedef Index_t IndexType;

    /*
     * \brief Default constructor set all data to zero initially
     */
    variables() :
            NI(0), NJ(0), NK(0), NI_m1(0), NJ_m1(0), NK_m1(0), NJNK(0), i_size(0) {
    }

    /*
     * \brief Default destructor
     */
    ~variables() {
        var.clear();
    }

    /*
     * \brief Resizes pseudo 'i', 'j' and 'k' indices
     */
    inline void init_rise(Index_t _NI, Index_t _NJ, Index_t _NK) {
        NI = _NI;
        NJ = _NJ;
        NK = _NK;
        if (_NI >= 1)
            NI_m1 = _NI - 1;
        else
            NI_m1 = 0;
        if (_NJ >= 1)
            NJ_m1 = _NJ - 1;
        else
            NJ_m1 = 0;
        if (_NK >= 1)
            NK_m1 = _NK - 1;
        else
            NK_m1 = 0;
        NJNK = NJ * NK;
    }

    /*
     * \brief Resizes pseudo 'i' and 'j' indices
     */
    inline void init_rise(Index_t _NI, Index_t _NJ) {
        NI = _NI;
        NJ = _NJ;
        NK = 1;
        if (_NI >= 1)
            NI_m1 = _NI - 1;
        else
            NI_m1 = 0;
        if (_NJ >= 1)
            NJ_m1 = _NJ - 1;
        else
            NJ_m1 = 0;
        NK_m1 = 0;
        NJNK = NJ * NK;
    }

    /*
     * \brief Resizes pseudo 'i' index
     */
    inline void init_rise(Index_t _NI) {
        NI = _NI;
        NJ = 1;
        NK = 1;
        if (_NI >= 1)
            NI_m1 = _NI - 1;
        else
            NI_m1 = 0;
        NJ_m1 = 0;
        NK_m1 = 0;
        NJNK = NJ * NK;
    }

    /*
     * \brief Allocates memory for a vector (the same as .allocate() method)
     * \note If new size is greater, the appended elements will not be initialized!
     */
    inline void resize(Index_t new_size) {

        /* Do not allocate if internal size is equal to the provided one */
        if (new_size == i_size) return;

        /*
         * Check if data will be allocated as 16-bytes aligned. If no - increase
         * vector size so that all internal points will be 16-bytes aligned
         */
        size_t size_of_T = sizeof(T);
        int checker = 0;
        checker = (new_size * size_of_T) % AlignmentN;
        if (checker != 0) {
            double real_block = new_size * size_of_T / (double)AlignmentN; // Calculate size of real block it it would be aligned
            int new_block = static_cast<int>((double)round(real_block + 0.5)); // Round up size of real block, such vector of new size
            Index_t new_size = new_block * AlignmentN / size_of_T;
            var.resize(new_size);                           // will be "aligned" by AlignmentN bytes
        }
        else
            var.resize(new_size);
        __it_beg = var.begin();
        __it_end = var.end();

//        if (new_size > i_size)
//            memset(var.data() + i_size, 0x00, sizeof(T) * (new_size - i_size));

        i_size = new_size;
    }

    /*
     * \brief Allocates memory for a pseudo 3d vector
     */
    inline void allocate(Index_t _NI, Index_t _NJ, Index_t _NK) {

        Index_t new_size = _NI * _NJ * _NK;

        if (new_size == i_size)
            return;

        this->init_rise(_NI, _NJ, _NK);
        this->resize(new_size);

        __it_beg = var.begin();
        __it_end = var.end();
    }

    /*
     * \brief Allocates memory for a 1d vector
     */
    inline void allocate(Index_t size) {
        this->init_rise(size);
        this->resize(size);
    }

    /*
     * \brief Adds element to the end of a vector
     */
    inline void add_elemet(T NewElement) {
        var.push_back(NewElement);
        __it_end = var.end();
        ++i_size;
    }

    /*
     * \brief Adds N elements from a given vector to the end of a current vector
     * TODO: add checker for the case if N > NewElements.size() and preallocate memory
     * of var
     */
    inline void add_elements(Index_t N, variables<T, AlignmentN> NewElements) {
        for(Index_t i = 0; i < N; ++i)
            var.push_back(NewElements(i));
        __it_end = var.end();
        i_size += N;
    }

    /*
     * \brief Adds N number of zero elements to the end of a vector
     * TODO: preallocate memory of var
     */
    inline void add_elements(Index_t N) {
        for(Index_t i = 0; i < N; ++i)
            var.push_back(0);
        __it_end = var.end();
        i_size += N;
    }

    /*
     * \brief Clears a vector
     * \note Clears the data but not necessarily frees the memory (for optimization purpose)!
     */
    inline void clear() {
        var.clear();
        init_rise(0);
        i_size = 0;
        NI = 0;
        NJ = 0;
        NK = 0;
        NJNK = 0;
    }

    /*
     * \brief Erases a vector
     * \note Completely frees the memory!
     */
    inline void erase() {
        var.clear();
        InternalVec(var).swap(var);
        init_rise(0);
        i_size = 0;
        NI = 0;
        NJ = 0;
        NK = 0;
        NJNK = 0;
    }

    /*
     * \brief Returns global id calculated from triplet of indices
     */
    inline Index_t pos(Index3 ind) {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return ind.k;
    }

    inline Index_t pos(Index3 ind) const {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return ind.k;
    }

    /*
     * \brief Returns global id calculated from three indices
     */
    inline Index_t pos(Index_t i, Index_t j, Index_t k) {
        k += NK * j;
        k += NJNK * i;
        return k;
    }

    inline Index_t pos(Index_t i, Index_t j, Index_t k) const {
        k += NK * j;
        k += NJNK * i;
        return k;
    }

    /*
     * \brief Returns global id calculated from two indices
     */
    inline Index_t pos(Index_t i, Index_t j) {
        j += NJ * i;
        return j;
    }

    inline Index_t pos(Index_t i, Index_t j) const {
        j += NJ * i;
        return j;
    }

    /*
     * \brief Returns global id calculated from triplet of indices
     */
    inline Index_t GetId(Index3 ind) {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return ind.k;
    }

    inline Index_t GetId(Index3 ind) const {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return ind.k;
    }

    inline Index_t GetId(Index_t i, Index_t j, Index_t k) {
        k += NK * j;
        k += NJNK * i;
        return k;
    }

    /*
     * \brief Returns element from the given by global index position
     * \note This method has different behavior than one from STL. It has no range check.
     */
    inline T &at(Index3 ind) {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return *(var.begin() + ind.k);
    }

    inline T at(Index3 ind) const {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return *(var.begin() + ind.k);
    }

    /*
     * \brief Returns element and global id from the given by two indices position
     * \note This method has different behavior than one from STL. It has no range check.
     */
    inline T &at(Index3 ind, Index_t &pos) {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        pos = ind.k;
        return *(var.begin() + pos);
    }

    inline T at(Index3 ind, Index_t &pos) const {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        pos = ind.k;
        return *(var.begin() + pos);
    }

    /*
     * \brief Returns element from the given by three indices position
     * \note This method has different behavior than one from STL. It has no range check.
     */
    inline T &at(Index_t i, Index_t j, Index_t k) {
        k += NK * j;
        k += NJNK * i;
        return *(var.begin() + k);
    }

    inline T at(Index_t i, Index_t j, Index_t k) const {
        k += NK * j;
        k += NJNK * i;
        return *(var.begin() + k);
    }

    /*!
     * \brief Returns element from the given by two indices position
     * \note This method has different behavior than one from STL. It has no range check.
     */
    inline T &at(Index_t i, Index_t j) {
        j += NJ * i;
        return *(var.begin() + j);
    }

    inline T at(Index_t i, Index_t j) const {
        j += NJ * i;
        return *(var.begin() + j);
    }

    /*!
     * \brief Returns element from the given position
     * \note This method has different behavior than one from STL. It has no range check.
     */
    inline T &at(Index_t id) {
        return *(var.begin() + id);
    }

    inline T at(Index_t id) const {
        return *(var.begin() + id);
    }

    /*
     * \brief Returns element from the given by global index position
     */
    inline T &operator [](Index_t n) {
        return *(var.begin() + n);
    }

    inline T operator [](Index_t n) const {
        return *(var.begin() + n);
    }

    /*!
     * \brief Returns element from the given by triplet of indices position
     */
    inline T &operator ()(Index3 ind) {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return *(var.begin() + ind.k);
    }

    inline T operator ()(Index3 ind) const {
        ind.k += NK * ind.j;
        ind.k += NJNK * ind.i;
        return *(var.begin() + ind.k);
    }

    /*!
     * \brief Returns element from the given by i-j-k triplet position
     */
    inline T &operator ()(Index_t i, Index_t j, Index_t k) {
        k += NK * j;
        k += NJNK * i;
        return *(var.begin() + k);
    }

    inline T operator ()(Index_t i, Index_t j, Index_t k) const {
        k += NK * j;
        k += NJNK * i;
        return *(var.begin() + k);
    }

    /*!
     * \brief Returns element and global id from the given by i-j-k triplet position
     */
    inline T &operator ()(Index_t i, Index_t j, Index_t k, Index_t &id) {
        k += NK * j;
        k += NJNK * i;
        id = k;
        return *(var.begin() + k);
    }

    /*!
     * \brief Returns element from the given by i-j pair position
     */
    inline T &operator ()(Index_t i, Index_t j) {  // only for 2d k = 0
        j += NJ * i;
        return *(var.begin() + j);
    }

    inline T operator ()(Index_t i, Index_t j) const {
        j += NJ * i;
        return *(var.begin() + j);
    }

    /*!
     * \brief Returns element from the given by global index position
     */
    inline T &operator ()(Index_t n) {
        return *(var.begin() + n);
    }

    inline T operator ()(Index_t n) const {
        return *(var.begin() + n);
    }

    /*!
     * \brief Back calculates triplet of indices based on provided \e id
     */
    inline Index3 GetTriplet(const Index_t id) {
        Index_t i = id / NJNK;
        Index_t j = (id - i * NJNK) / NK;
        Index_t k = id - NK * (j + NJ * i);

        return {i, j, k};
    }

    inline Index3 GetTriplet(const Index_t id) const {
        Index_t i = id / NJNK;
        Index_t j = (id - i * NJNK) / NK;
        Index_t k = id - NK * (j + NJ * i);

        return {i, j, k};
    }

    /*!
     * \brief Back calculates triplet of indices based on provided \e id
     */
    inline void GetTriplet(const Index_t id, Index_t &i, Index_t &j, Index_t &k) {
        i = id / NJNK;
        j = (id - i * NJNK) / NK;
        k = id - NK * (j + NJ * i);
    }

    inline void GetTriplet(const Index_t id, Index_t &i, Index_t &j, Index_t &k) const {
        i = id / NJNK;
        j = (id - i * NJNK) / NK;
        k = id - NK * (j + NJ * i);
    }

    /*!
     * \brief Back calculates pair of indices based on provided \e id
     */
    inline Index3 GetPair(const Index_t id) {
        Index_t i = id / NJ;
        Index_t j = id - i * NJNK;

        return {i, j, 0};
    }

    inline Index3 GetPair(const Index_t id) const {
        Index_t i = id / NJ;
        Index_t j = id - i * NJNK;

        return {i, j, 0};
    }

    /*!
     * \brief Back calculates pair of indices based on provided \e id
     */
    inline void GetPair(const Index_t id, Index_t &i, Index_t &j) {
        i = id / NJ;
        j = id - i * NJNK;
    }

    inline void GetPair(const Index_t id, Index_t &i, Index_t &j) const {
        i = id / NJ;
        j = id - i * NJNK;
    }

    inline Index3 GetExtent() {
        return {NI, NJ, NK};
    }

    inline Index3 GetExtent() const {
        return {NI, NJ, NK};
    }

    /*!
     * \brief Returns number of elements
     * \note Method returns number of elements based on provided dimensions, not true size of the
     * internal vector
     */
    inline Index_t size() {
        return i_size;
    }

    /*!
     * \brief Returns number of elements
     * \note Method returns number of elements based on provided dimensions, not true size of the
     * internal vector
     */
    inline Index_t size() const {
        return i_size;
    }

    /*!
     * \brief Returns size of this class in Bytes
     */
    inline size_t ClassSize() {
        return sizeof(variables);
    }

    inline size_t ClassSize() const {
        return sizeof(variables);
    }

    /*!
     * \brief Returns number of elements along pseudo 'k' direction
     * \note It will return 1 for 2d case
     */
    inline Index_t NK_init() {
        return NK;
    }

    inline Index_t NK_init() const {
        return NK;
    }

    /*!
     * \brief Returns number of elements along pseudo 'j' direction
     */
    inline Index_t NJ_init() {
        return NJ;
    }

    inline Index_t NJ_init() const {
        return NJ;
    }

    /*!
     * \brief Returns number of elements along pseudo 'i' direction
     */
    inline Index_t NI_init() {
        return NI;
    }

    inline Index_t NI_init() const {
        return NI;
    }

    /*!
     * \brief Returns reference to the internal stl vector of data
     */
    InternalVec &ref() {
        return var;
    }

    /*!
     * \brief Returns a pointer to the specified position
     */
    T *ref(Index_t shift) {
        return var.data() + shift;
    }

    /*!
     * \brief Returns pointer to the internal data
     */
    T *data() {
        return var.data();
    }

    const T *data() const {
        return var.data();
    }

    /*!
     * \brief Returns true if pseudo 3d space is used
     */
    inline bool Is3d() {
        if (NK == 0 || NK == 1)
            return false;
        else
            return true;
    }

    inline bool Is3d() const {
        if (NK == 0 || NK == 1)
            return false;
        else
            return true;
    }

    /*!
     * \brief Sets internal vector to zero
     */
    inline void setZero() {
        for(Index_t n = 0; n < i_size; ++n) {
            var[n] = (T)0;
        }
    }

    /*!
     * \brief Assigns a constant value to all elements of the vector
     */
    inline void setConstant(T value) {
        for(Index_t n = 0; n < i_size; ++n) {
            var[n] = value;
        }
    }

    /*!
     * \brief Sorts all elements in a vector
     */
    void Sort() {
        sort(var.begin(), var.end());
    }

    /*!
     * \brief Returns maximum element in a vector and its position
     */
    template<typename P>
    inline T FindMax(P &Position) {
        __it = __it_beg;
        ++__it_beg;

        for(; __it_beg != __it_end; ++__it_beg) {
            if (*__it < *__it_beg) {
                __it = __it_beg;
            }
        }

        __it_beg = var.begin();
        Position = distance(__it_beg, __it);
        return *__it;
    }

    inline variables<T> &operator = (const variables<T> &other) {
        if (&other == this)
            return *this;

        if (i_size == other.size()) {
            for(Index_t n = 0; n < i_size; ++n)
                var[n] = other[n];
            return *this;
        }

        this->clear();
        this->allocate(other.NI_init(), other.NJ_init(), other.NK_init());
        for(Index_t n = 0; n < i_size; ++n)
            var[n] = other[n];
        return *this;
    }
private:

    /*
     * Returns iterator to the maximum element within specified range
     */
    template<typename _ForwardIterator>
    _ForwardIterator my_max_element_changed(_ForwardIterator __first, _ForwardIterator __last) {
        if (__first == __last) return __first;
        _ForwardIterator __result = __first;
        ++__first;

        for(; __first != __last; ++__first)
            if (*__result < *__first) __result = __first;

        return __result;
    }

private:
    InternalVec var;
    typename InternalVec::iterator __it;
    typename InternalVec::iterator __it_beg;
    typename InternalVec::iterator __it_end;
    Index_t NI;
    Index_t NJ;
    Index_t NK;
    Index_t NI_m1;
    Index_t NJ_m1;
    Index_t NK_m1;
    Index_t NJNK;
    Index_t i_size;     // internal size of the vector, based on provided dimensions
};
}



#endif /* CUSTOMARRAY_H_ */
