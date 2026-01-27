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

    m_valOutside =
        addIntParameter("value_outside", "value to be used if target is outside source domain", NaN, Parameter::Choice);
    m_userDef = addFloatParameter("user_defined_value", "user defined value if target outside source domain", 1.0);

    V_ENUM_SET_CHOICES_SCOPE(m_valOutside, ValueOutside, );
}

std::unique_ptr<viskores::filter::Filter> SampleVtkm::setUpFilter() const {
    viskores::cont::DataSet refCopy;
    bool refValidCopy = false;
    {
        std::lock_guard<std::mutex> lock(m_refMutex);
        refCopy = m_ref_in;
        refValidCopy = m_refValid;
    }

    if (!refValidCopy || refCopy.GetNumberOfCoordinateSystems() == 0) {
        sendError("SampleVtkm::setUpFilter: No valid reference geometry set");
        return nullptr;
    }

    auto filter = std::make_unique<viskores::filter::resampling::Probe>();
    filter->SetGeometry(refCopy);
    filter->SetInvalidValue(getInvalidValue());
    return filter;
}

bool SampleVtkm::objectAdded(int sender, const std::string &senderPort, const Port *port)
{
    if (port->getName() == "ref_in") {
        viskores::cont::DataSet local;
        bool localValid = false;

        m_ref_in = viskores::cont::DataSet{};
        m_refValid = false;

        if (port->objects().empty()) {
            sendError("No object on ref_in port");
            return true;
        }
        
        auto object = port->objects().back();
        if (!object) {
            sendError("Invalid object on ref_in port");
            return true;
        }
        
        if (object->getTimestep() < 1) {
            // Try to get it as a GeometryInterface directly (like a Coords object)
            const GeometryInterface *geo = object->getInterface<GeometryInterface>();
            if (geo) {
                ModuleStatusPtr status = vistle::vtkmSetGrid(local, object);
                if (status) {
                    localValid = local.GetNumberOfCoordinateSystems() > 0;
                    if (localValid) {
                        // Check cell set type
                        auto cellSet = local.GetCellSet();
                        if (cellSet.IsType<viskores::cont::CellSetSingleType<>>()) {
                            auto singleType = cellSet.AsCellSet<viskores::cont::CellSetSingleType<>>();
                            if (singleType.GetCellShapeAsId() == 1) { // VTK_VERTEX
                                sendError("Reference grid contains vertex cells which are not supported by Probe filter");
                                localValid = false;
                            }
                        }
                    }
                    if (!localValid && local.GetNumberOfCoordinateSystems() > 0)
                        sendError("Reference grid has no coordinate system or unsupported cell type");
                } else {
                    sendError("Failed to convert grid to viskores format");
                }
            } else {
                // Try to extract grid from DataBase wrapper
                DataBase::const_ptr data = DataBase::as(object);
                if (data && data->grid()) {
                    ModuleStatusPtr status = vistle::vtkmSetGrid(local, data->grid());
                    if (status) {
                        localValid = local.GetNumberOfCoordinateSystems() > 0;
                        if (!localValid)
                            sendError("Reference grid has no coordinate system");
                    } else {
                        sendError("Failed to convert grid to viskores format");
                    }
                } else {
                    sendError("Could not extract grid from object on ref_in port");
                }
            }
        }

        {
            std::lock_guard<std::mutex> lock(m_refMutex);
            m_ref_in = std::move(local);
            m_refValid = localValid;
        }
    }
    return true;
}

ModuleStatusPtr SampleVtkm::prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                              const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                              viskores::cont::DataSet &dataset) const
{
    // Test why field is not valid:
    if (!field) {
        sendError("Input field is null");
        return Error("Input field is null");
    }
    return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);
}

/* vistle::Object::const_ptr SampleVtkm::prepareOutputGrid(const viskores::cont::DataSet &dataset,
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
 */
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
