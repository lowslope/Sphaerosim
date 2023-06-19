#include "Exception.h"

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

namespace SphaeroSim {

    Exception::Exception(const std::string &type_description,
                         const std::string &error_text) :
            std::runtime_error(error_text) {
        _type_description = type_description;
        _error_text = error_text;
    }

}