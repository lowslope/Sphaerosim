// Copyright 2014 Marcel Spekowius
#ifndef _SPHAEROSIM_HASH_VALUE_H_
#define _SPHAEROSIM_HASH_VALUE_H_

#include <stdint.h>

namespace SphaeroSim {

    class HashValue {
    public:
        // singleton
        static HashValue *I();

        static void FreeSingleton();

        uint64_t Recalculate(const uint64_t old_hash, const uint64_t new_data);

    private:
        HashValue();

        union hash_value_union {
            uint64_t value;
            uint32_t dwords[2];
            uint16_t words[4];
        };

        void InitHashTable();

        uint64_t hash_keys_[0x10000 + 3];

        static HashValue *global_instance_;
    };

}
#endif  // _SPHAEROSIM_HASH_VALUE_H_
