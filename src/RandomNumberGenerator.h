#ifndef _SPHAEROSIM_RANDOM_NUMBER_GENERATOR_H_
#define _SPHAEROSIM_RANDOM_NUMBER_GENERATOR_H_

#ifndef _INCLUDED_VECTOR_H_

#include <vector>

#define _INCLUDED_VECTOR_H_
#endif

#define R   32      // number of words in the state register

#define M1   3
#define M2  24
#define M3  10

#define MAT0POS(t, v)  (v ^ (v >>  (t)))
#define MAT0NEG(t, v)  (v ^ (v << -(t)))
#define Identity(v)   (v)

#define V0            states_[ current_index_        ]
#define VM1           states_[(current_index_+M1) % R]
#define VM2           states_[(current_index_+M2) % R]
#define VM3           states_[(current_index_+M3) % R]
#define VRm1          states_[(current_index_+31) % R]

#define newV0         states_[(current_index_+31) % R]
#define newV1         states_[ current_index_        ]

// pseudo random number generator with the 'Well Equidistributed Long-period Linear' (WELL) algorithm
// taken from ACM Transactions on Mathematical Software, Vol. V, No. N, Month 20YY, Pages 1�14.
class RandomNumberGenerator {
public:
    // constructor & destructor
    RandomNumberGenerator();

    // generate a random number
    unsigned int rand() {
        unsigned int z0;
        unsigned int z1;
        unsigned int z2;

        z0 = VRm1;
        z1 = Identity(V0) ^ MAT0POS (+8, VM1);
        z2 = MAT0NEG (-19, VM2) ^ MAT0NEG (-14, VM3);
        newV1 = z1 ^ z2;
        newV0 = MAT0NEG (-11, z0) ^ MAT0NEG (-7, z1) ^ MAT0NEG (-13, z2);
        current_index_ = (current_index_ + R - 1) % R;

        return states_[current_index_];
    }

    // set the state registers
    void SetStateRegisters(const std::vector<unsigned int> &states);

    // get the state registers
    void GetStateRegisters(std::vector<unsigned int> *dst) const;

private:
    // state registers
    unsigned int states_[R];
    // current index
    int current_index_;
};

#endif  // _SPHAEROSIM_RANDOM_NUMBER_GENERATOR_H_
