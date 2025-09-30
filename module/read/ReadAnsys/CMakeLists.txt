set(HEADERS
    ReadAnsys.h
    ReadRST.h
    DerivedData.h
    DOFData.h
    DOFOptions.h
    Element.h
    GeometryHeader.h
    Map1D.h
    Node.h
    RstHeader.h
    SolutionHeader.h
    ANSYS.h)

set(SOURCES ReadAnsys.cpp ReadRST.cpp Map1D.cpp Node.cpp ReadRSTDerived.cpp ANSYS.cpp)

add_module(ReadAnsys "read ANSYS data" ${HEADERS} ${SOURCES})
