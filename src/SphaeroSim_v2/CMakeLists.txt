file(GLOB SOURCES_CPP "*.cpp")
file(GLOB SOURCES_H "*.h")
set(SOURCES_Sphaerosim "${SOURCES_CPP}\n${SOURCES_H}" PARENT_SCOPE)
#add_executable(Sphaerosim ${SOURCES_Sphaerosim})