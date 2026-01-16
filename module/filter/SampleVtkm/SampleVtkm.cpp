//#include <viskores/include/viskores/filter/resampling/viskores_filter_resampling_export.h>
#include <viskores/filter/resampling/Probe.h>

#include "SampleVtkm.h"

MODULE_MAIN(SampleVtkm)

useing namespace vistle;

SampleVtkm::SampleVtkm(const std::string &name, int moduleID, mpi::communicator comm, int numPorts,
                       bool requireMappedData)
: vistle::vtkm::VtkmModule(name, moduleID, comm, numPorts, requireMappedData)
{
    createInputPort("data_in", "object with data to be sampled");
    createInputPort("ref_in", "target grid", Port::Flags::NOCOMPUTE);

    m_out = createOutputPort("data_out", "data sampled to target grid");

    m_mode = addIntParameter("mode", "interpolation mode", Linear, Parameter::Choice);
    m_createCelltree = addIntParameter("create_celltree", "create celltree", 0, Parameter::Boolean);
    m_valOutside = addIntParameter("value_outside", "value to be used if target is outside source domain", Linear,
                                   Parameter::Choice);
    m_userDef = addFloatParameter("user_defined_value", "user defined value if target outside source domain", 1.0);
    m_hits = addIntParameter("mulit_hits", "handling if multiple interpolatied values found for one target point ",
                             Linear, Parameter::Choice);
    V_ENUM_SET_CHOICES(m_mode, InterpolationMode);
    V_ENUM_SET_CHOICES(m_valOutside, ValueOutside);
    V_ENUM_SET_CHOICES(m_hits, MultiHits);
}

std::unique_ptr<vtkm::filter::Filter> SampleVtkm::setUpFilter() const
{
    auto filter = std::make_unique<vtkm::filter::resampling::Probe>();
    return filter;
}

ModuleStatusPtr SampleVtkm::prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                              const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                              viskores::cont::DataSet &dataset) const
{
    // Implement any specific input field preparation if needed
    return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);
}

vistle::Object::const_ptr SampleVtkm::prepareOutputGrid(const viskores::cont::DataSet &dataset,
                                                        const vistle::Object::const_ptr &inputGrid) const
{
    // Implement output grid preparation if needed
    return nullptr;
}

vistle::DataBase::ptr SampleVtkm::prepareOutputField(const viskores::cont::DataSet &dataset,
                                                     const vistle::Object::const_ptr &inputGrid,
                                                     const vistle::DataBase::const_ptr &inputField,
                                                     const std::string &fieldName,
                                                     const vistle::Object::const_ptr &outputGrid) const
{
    // Implement output field preparation if needed
    return nullptr;
}
