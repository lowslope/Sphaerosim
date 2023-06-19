#include_directories(Eigen)
add_subdirectory(Eigen)
add_subdirectory(SphaeroSim_v2)

STRING(REGEX REPLACE "\n" ";" SOURCES_Eigen "${SOURCES_Eigen}")
STRING(REGEX REPLACE "\n" ";" SOURCES_Sphaerosim "${SOURCES_Sphaerosim}")
STRING(REGEX REPLACE "\n" ";" DIR_Eigen "${DIR_Eigen}")

foreach (SOURCE ${SOURCES_Eigen})
    set(SOURCES "${SOURCES}${CMAKE_CURRENT_SOURCE_DIR}/Eigen/${SOURCE};")
endforeach (SOURCE)
foreach (SOURCE ${SOURCES_Sphaerosim})
    #set(SOURCES "${SOURCES};SphaeroSim_v2/${SOURCE}")
    set(SOURCES "${SOURCES}${SOURCE};")
endforeach (SOURCE)

foreach (SOURCE ${DIR_Eigen})
    #set(SOURCES "${SOURCES};SphaeroSim_v2/${SOURCE}")
    set(INCLUDE "${INCLUDE}${CMAKE_CURRENT_SOURCE_DIR}/Eigen/${SOURCE};")
endforeach (SOURCE)

STRING(REGEX REPLACE "/./" "/" SOURCES "${SOURCES}")
STRING(REGEX REPLACE "/./" "/" INCLUDE "${INCLUDE}")

set(SOURCES "${SOURCES}" PARENT_SCOPE)
set(INCLUDE "${INCLUDE}" PARENT_SCOPE)