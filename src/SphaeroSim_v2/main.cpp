#ifndef _INCLUDED_STDIO_H_

#include <stdio.h>

#define _INCLUDED_STDIO_H_
#endif

#ifndef _INCLUDED_STRING_H_

#include <string>

#define _INCLUDED_STRING_H_
#endif

#ifndef _SPHAEROSIM_APPLICATION_H_

#include "Application.h"

#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_
#include "Exception.h"
#endif

#ifndef _UTILITY_HASH_VALUE_H_

#include "HashValue.h"

#endif

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

//#ifndef INSPECTOR
//#define INTEL_NO_ITTNOTIFY_API
//#endif

#include <set>

#include "RandomNumberGenerator.h"

#include "ComsolFile.h"

/*
SphaeroSim::RandomNumberGenerator* rng_;
bool MonteCarlo(const double prob) {
  double comparison;
#pragma omp critical
  {
  comparison = rng_->rand() /
           static_cast<double>(std::numeric_limits<unsigned int>::max());
  }
  if (comparison < prob)
    return true;
  return false;
}*/
inline float sqrt7(float x) {
    unsigned int i = *(unsigned int *) &x;
    // adjust bias
    i += 127 << 23;
    // approximation of square root
    i >>= 1;
    return *(float *) &i;
}


inline float sqrt3(const float x) {
    union {
        int i;
        float x;
    } u;

    u.x = x;
    u.i = (1 << 29) + (u.i >> 1) - (1 << 22);
    return u.x;
}

#define SQRT_MAGIC_F 0x5f3759df

inline float sqrt2(const float x) {
    const float xhalf = 0.5f * x;

    union // get bits for floating value
    {
        float x;
        int i;
    } u;
    u.x = x;
    u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
    return x * u.x * (1.5f - xhalf * u.x * u.x);// Newton step, repeating increases accuracy
}

int main(int num_arguments, char *arguments[]) {

    /*for (int64_t i = 0; i < 10000; ++i) {
      const double the_val = 0.1231 * (i + 100);
      const double res = sqrt(the_val);
      const double res_2 = sqrt2(static_cast<float>(the_val));
      if (std::abs(res - res_2) > 1e-6)
        printf("Error\n");
    }*/
/*  rng_ = new SphaeroSim::RandomNumberGenerator();
  const double rate = 1;
  const double timestep = 1e-4;
  const double time = 10;
  const int num_experiments = 1000;

  std::vector<int64_t> num_events;
  num_events.resize(num_experiments, 0);

  
  const int64_t num_steps = (int64_t)(time / timestep + 0.5);
#pragma omp parallel for
  for (int64_t cur_exp = 0; cur_exp < num_experiments; ++cur_exp) {
    for (int64_t cur_step = 0; cur_step < num_steps; ++cur_step) {
      const double prob = timestep * rate;
      if (MonteCarlo(prob) == true)
        num_events[cur_exp] += 1;
    }
  }
  */
    /*
    std::multiset<int> myset;
    std::multiset<int>::iterator itlow,itup;

    for (int i=1; i<10; i++) myset.insert(10); // 10 20 30 40 50 60 70 80 90

    itlow=myset.lower_bound (5);                //       ^
    itup=myset.upper_bound (60);

    myset.erase(myset.begin(), myset.lower_bound(11));
    */

    // get the parameterfile
    LIKWID_MARKER_INIT;
    #pragma omp parallel
    {
        LIKWID_MARKER_THREADINIT;
        //__itt_suppress_push(__itt_suppress_memory_errors | __itt_suppress_threading_errors);

    };
    srand(66);
    LIKWID_MARKER_REGISTER("ChangeCellToSolid");
    std::string parameter_file = "";
    for (int cur_arg = 0; cur_arg < num_arguments; ++cur_arg) {
        std::string arg = arguments[cur_arg];
        if (arg == "-f") {
            if (cur_arg != num_arguments - 1) {
                parameter_file = arguments[cur_arg + 1];
                break;
            }
        }
    }

    try {
        // initialize the singletons
        SphaeroSim::HashValue::I();

        SphaeroSim::Application app(parameter_file);
        app.Run();

        SphaeroSim::HashValue::FreeSingleton();
    } catch (const SphaeroSim::Exception &exc) {
        printf("\nERROR:\n");
        printf("  %s\n", exc.GetTypeDescription().c_str());
        printf("  %s\n\n", exc.GetErrorText().c_str());
        system("pause");
    }
        #pragma omp taskwait
    LIKWID_MARKER_CLOSE;
    return 0;
}