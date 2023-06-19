#ifndef _SPHAEROSIM_EXCEPTION_H_
#define _SPHAEROSIM_EXCEPTION_H_

#ifndef _INCLUDED_STRING_H_

#include <string>

#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_STDEXCEPT_H_

#include <stdexcept>

#define _INCLUDED_STDEXCEPT_H_
#endif

namespace SphaeroSim {

    class Exception : public std::runtime_error {
    public:
        Exception(const std::string &type_description,
                  const std::string &error_text);

        // Getter
        inline std::string GetTypeDescription() const {
            return _type_description;
        }

        inline std::string GetErrorText() const {
            return _error_text;
        }

    private:
        std::string _type_description;
        std::string _error_text;
    };

}  // namespace SphaeroSim
#endif  // _SPHAEROSIM_EXCEPTION_H_
