//#include <viskores/include/viskores/filter/resampling/viskores_filter_resampling_export.h>
#include <viskores/filter/resampling/Probe.h>
#include <vistle/core/object.h>

#include "SampleVtkm.h"

using namespace vistle;

MODULE_MAIN(SampleVtkm)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(InterpolationMode,
                                    (First) // value of first vertex
                                    (Mean) // mean value of all vertices
                                    (Nearest) // value of nearest vertex
                                    (Linear) // barycentric/multilinear interpolation
)
DEFINE_ENUM_WITH_STRING_CONVERSIONS(ValueOutside, (NaN)(Zero)(userDefined))
DEFINE_ENUM_WITH_STRING_CONVERSIONS(MultiHits, (Any)(Average))

SampleVtkm::SampleVtkm(
    const std::string &name, int moduleID, mpi::communicator comm, int numPorts,
    bool requireMappedData) //TODO: numPorts 2 input and 1 output ports - just one input port needs mapped data
: VtkmModule(name, moduleID, comm)
{
    createInputPort("data_in", "object with data to be sampled"); // TODO: input of the Execute function
    createInputPort("ref_in", "target grid", Port::Flags::NOCOMPUTE); // TODO: Use get geometry to set as target drif

    m_out = createOutputPort("data_out", "data sampled to target grid");

    m_mode = addIntParameter("mode", "interpolation mode", Linear, Parameter::Choice);
    m_createCelltree = addIntParameter("create_celltree", "create celltree", 0, Parameter::Boolean);
    m_valOutside = addIntParameter("value_outside", "value to be used if target is outside source domain", Linear,
                                   Parameter::Choice); //<todo: use SetInvalidValue of the Probe filter to specify this
    m_userDef = addFloatParameter("user_defined_value", "user defined value if target outside source domain", 1.0);
    m_hits = addIntParameter("mulit_hits", "handling if multiple interpolatied values found for one target point ",
                             Linear, Parameter::Choice);

    V_ENUM_SET_CHOICES_SCOPE(m_mode, InterpolationMode, );
    V_ENUM_SET_CHOICES_SCOPE(m_valOutside, ValueOutside, );
    V_ENUM_SET_CHOICES_SCOPE(m_hits, MultiHits, );
}

std::unique_ptr<viskores::filter::Filter> SampleVtkm::setUpFilter() const
{
    auto filter = std::make_unique<viskores::filter::resampling::Probe>();
    // TODO: check if mapped data exists in second input port
    // TODO: set up the target geometry for the Probe filter:
    //filter->SetGeometry(ref_in->getObject()->getInterface<viskores::cont::DataSet>());
    filter->SetInvalidValue(m_valOutside->getValue());
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
