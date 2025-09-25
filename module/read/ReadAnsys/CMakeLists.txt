set(HEADERS ReadAnsys.h)

set(SOURCES ReadAnsys.cpp)

add_module(ReadAnsys "read ANSYS data" ${HEADERS} ${SOURCES})
