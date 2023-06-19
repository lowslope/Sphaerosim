#ifndef _SPHAEROSIM_CACHE_H_
#define _SPHAEROSIM_CACHE_H_

#ifndef _INCLUDED_MAP_H_
#include <map>
#define _INCLUDED_MAP_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_OMP_H_

#include <omp.h>

#define _INCLUDED_OMP_H_
#endif


namespace SphaeroSim {

// simple cache which uses a map for the key/value assignment
    template<typename key_type, typename value_type>
    class Cache {
    public:
        explicit Cache(const int32_t number_entries) {
            int max_threads = omp_get_max_threads();
            cache_data_.resize(max_threads);

            for (int cur_thread = 0; cur_thread < max_threads; ++cur_thread) {
                cache_data_[cur_thread].number_entries_ = number_entries / max_threads;
                cache_data_[cur_thread].erase_counter_ = 0;
                cache_data_[cur_thread].key_counter_ = 0;
                cache_data_[cur_thread].keys_.
                        resize(cache_data_[cur_thread].number_entries_);
            }
        }

        inline bool CachedValue(const key_type &key,
                                value_type *dst) const {
            int thread_id = omp_get_thread_num();
            typename std::map<key_type, value_type>::const_iterator it =
                    cache_data_[thread_id].data_.find(key);
            if (it == cache_data_[thread_id].data_.end()) {
                return false;
            }
            *dst = it->second;
            return true;
        }

        inline void StoreValue(const key_type &key, const value_type &value) {
            if (cache_data()->data_.size() == cache_data()->number_entries_) {
                cache_data()->data_.erase(cache_data()->keys_[cache_data()->erase_counter_]);
                ++cache_data()->erase_counter_;
                if (cache_data()->erase_counter_ >= cache_data()->number_entries_) {
                    cache_data()->erase_counter_ = 0;
                }
            }
            cache_data()->data_[key] = value;

            cache_data()->keys_[cache_data()->key_counter_] = key;
            ++cache_data()->key_counter_;
            if (cache_data()->key_counter_ >= cache_data()->number_entries_) {
                cache_data()->key_counter_ = 0;
            }
        }

    private:
        struct CacheData {
            int32_t number_entries_;
            std::vector<key_type> keys_;
            int32_t erase_counter_;
            int32_t key_counter_;
            std::map<key_type, value_type> data_;
        };

        inline CacheData *cache_data() {
            return &cache_data_[omp_get_thread_num()];
        }

        // one cachedataentry for each thread
        std::vector<CacheData> cache_data_;

        Cache &operator=(const Cache &);
    };

}

#endif  // _SPHAEROSIM_CACHE_H_
