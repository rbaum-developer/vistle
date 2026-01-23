//#include <viskores/include/viskores/filter/resampling/viskores_filter_resampling_export.h>
#include <viskores/filter/resampling/Probe.h>
#include <vistle/core/object.h>
#include <vistle/vtkm/convert.h>

#include "SampleVtkm.h"

using namespace vistle;

MODULE_MAIN(SampleVtkm)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(ValueOutside, (NaN)(Zero)(userDefined))

SampleVtkm::SampleVtkm(
    const std::string &name, int moduleID, mpi::communicator comm, int numPorts,
    bool requireMappedData) //TODO: numPorts 2 input and 1 output ports - just one input port needs mapped data
: VtkmModule(name, moduleID, comm)
{
    createInputPort("ref_in", "target grid", Port::Flags::NOCOMPUTE);

    m_valOutside = addIntParameter("value_outside", "value to be used if target is outside source domain", NaN,
                                   Parameter::Choice);
    m_userDef = addFloatParameter("user_defined_value", "user defined value if target outside source domain", 1.0);

    V_ENUM_SET_CHOICES_SCOPE(m_valOutside, ValueOutside, );
}

std::unique_ptr<viskores::filter::Filter> SampleVtkm::setUpFilter() const
{
    auto filter = std::make_unique<viskores::filter::resampling::Probe>();
    // set up the target geometry for the Probe filter:
    filter->SetGeometry(m_ref_in);

    // set up handling of values outside the input domain:
    float valOut = getInvalidValue();
    filter->SetInvalidValue(valOut);
    return filter;
}

ModuleStatusPtr SampleVtkm::prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                              const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                              viskores::cont::DataSet &dataset) const
{
    if (port->getName() == "ref_in") {
        auto object = port->objects().back();
        auto split = splitContainerObject(object);
        auto coords = Coords::as(split.geometry);
        if (object->getTimestep() < 1 and coords) {
            // convert to viskores data set
            ModuleStatusPtr status = vistle::vtkmSetGrid(m_ref_in, coords);
            if (!status) {
                return Error("Failed to set grid in dataset");
            }
        } else {
            return Error("No valid grid object received on ref_in port");
        }
    }
    return Success();
}

vistle::Object::const_ptr SampleVtkm::prepareOutputGrid(const viskores::cont::DataSet &dataset,
                                                        const vistle::Object::const_ptr &inputGrid) const
{
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

float SampleVtkm::getInvalidValue() const
{
    if (m_valOutside->getValue() == NaN) {
        return NAN;
    }
    if (m_valOutside->getValue() == Zero) {
        return 0.0f;
    }
    if (m_valOutside->getValue() == userDefined) {
        const float v = static_cast<float>(m_userDef->getValue());
        if (v)
            return v;
    }
    return 0.0f;
}
