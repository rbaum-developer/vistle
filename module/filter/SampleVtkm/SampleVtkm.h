#ifndef VISTLE_SAMPLEVTKM_SAMPLEVTKM_H
#define VISTLE_SAMPLEVTKM_SAMPLEVTKM_H

#include <vistle/vtkm/vtkm_module.h>

#include <vistle/module/module.h>
#include <vistle/core/object.h>
#include <vistle/core/grid.h>
#include <vistle/core/uniformgrid.h>
#include <vistle/core/structuredgrid.h>
#include <vistle/core/rectilineargrid.h>
#include <vistle/core/unstr.h>

class SampleVtkm: public vistle::VtkmModule {
public:
    SampleVtkm(const std::string &name, int moduleID, mpi::communicator comm, int numPorts = 1,
               bool requireMappedData = false);

private:
    ModuleStatusPtr prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                      const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                      viskores::cont::DataSet &dataset) const override;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;

    vistle::Object::const_ptr prepareOutputGrid(const viskores::cont::DataSet &dataset,
                                                const vistle::Object::const_ptr &inputGrid) const override;

    vistle::DataBase::ptr prepareOutputField(const viskores::cont::DataSet &dataset,
                                             const vistle::Object::const_ptr &inputGrid,
                                             const vistle::DataBase::const_ptr &inputField,
                                             const std::string &fieldName,
                                             const vistle::Object::const_ptr &outputGrid) const override;
    float getInvalidValue() const;

    vistle::IntParameter *m_valOutside;
    vistle::FloatParameter *m_userDef;

    mutable viskores::cont::DataSet m_ref_in;

    bool m_useCelltree = false;
};

#endif // VISTLE_SAMPLEVTKM_SAMPLEVTKM_H
