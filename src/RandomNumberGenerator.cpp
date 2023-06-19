#include "RandomNumberGenerator.h"

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_TIME_H_

#include <time.h>

#define _INCLUDED_TIME_H_
#endif

#ifndef _INCLUDED_STDLIB_H_

#include <stdlib.h>

#define _INCLUDED_STDLIB_H_
#endif

RandomNumberGenerator::RandomNumberGenerator() {
    srand((unsigned int) time(NULL));

    union value {
        unsigned short two_byte_values[2];
        unsigned int four_byte_value;
    };
    for (int j = 0; j < R; j++) {
        value new_val;
        new_val.two_byte_values[0] = (unsigned short) ::rand() + (unsigned short) ::rand();
        new_val.two_byte_values[1] = (unsigned short) ::rand() + (unsigned short) ::rand();
        states_[j] = new_val.four_byte_value;
    }
    current_index_ = 0;
}

// set the state registers
void RandomNumberGenerator::
SetStateRegisters(const std::vector<unsigned int> &states) {
    if (states.size() != R) {
        throw std::runtime_error("The number of statevalues must equal 32");
    }

    for (std::size_t cur = 0; cur < states.size(); ++cur)
        states_[cur] = states[cur];
}

// get the state registers
void RandomNumberGenerator::
GetStateRegisters(std::vector<unsigned int> *dst) const {
    for (int j = 0; j < R; ++j)
        dst->push_back(states_[j]);
}

