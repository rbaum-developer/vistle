add_module(Sample "sample data on points, unstructured and uniform grids to a uniform grid" Sample.cpp)

add_vtkm_module(SampleVtkm "maps the fields of one dataset to another" Sample.cpp Sample.h)
target_compile_definitions(SampleVtkm PRIVATE Sample=SampleVtkm)
target_compile_definitions(SampleVtkm PRIVATE SAMPLEVTKM)