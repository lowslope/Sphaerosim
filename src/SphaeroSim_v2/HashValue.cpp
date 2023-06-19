// Copyright 2014 Marcel Spekowius
#include "HashValue.h"

#include <stdint.h>
#include <vector>

#include "RandomNumberGenerator.h"

#ifndef NULL
#define NULL 0
#endif

namespace SphaeroSim {

    HashValue *HashValue::global_instance_ = NULL;

    HashValue *HashValue::I() {
        if (global_instance_ == NULL) {
            global_instance_ = new HashValue();
        }
        return global_instance_;
    }

    void HashValue::FreeSingleton() {
        if (global_instance_ != NULL) {
            delete global_instance_;
        }
    }

    HashValue::HashValue() {
        InitHashTable();
    }

    void HashValue::InitHashTable() {
        SphaeroSim::RandomNumberGenerator *rng =
                new SphaeroSim::RandomNumberGenerator();

        union {
            uint32_t dwords[2];
            uint64_t qword;
        } conversion;

        for (uint32_t counter = 0; counter < 0x10003; ++counter) {
            conversion.dwords[0] = rng->rand();
            conversion.dwords[1] = rng->rand();
            hash_keys_[counter] = conversion.qword;
        }

        delete rng;
    }

    uint64_t HashValue::
    Recalculate(const uint64_t old_hash, const uint64_t new_data) {
        hash_value_union h;
        h.value = new_data ^ old_hash;

        uint64_t result;
        result = old_hash ^ hash_keys_[h.words[0]];
        result = result ^ hash_keys_[((h.words[1] + 1) & 0xFFFF)];
        result = result ^ hash_keys_[((h.words[2] + 2) & 0xFFFF)];
        result = result ^ hash_keys_[((h.words[3] + 3) & 0xFFFF)];
        return result;
    }

}