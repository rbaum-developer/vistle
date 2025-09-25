#ifndef VISTLE_READANSYS_READANSYS_H
#define VISTLE_READANSYS_READANSYS_H

#include <vistle/module/reader.h>

#include <vector>
#include <string>

#include "ReadRST.h"
#include "Map1D.h"
#include "DOFOptions.h"

class ReadAnsys: public vistle::Reader {
public:
    ReadAnsys(const std::string &name, int moduleID, mpi::communicator comm);
    ~ReadAnsys() override;

    // reader interface
    bool examine(const vistle::Parameter *param) override;
    bool read(vistle::Reader::Token &token, int timestep = -1, int block = -1) override;
    bool prepareRead() override;
    bool finishRead() override;
    static const int V_OFFSET = 500;
    static const int EX_OFFSET = 1000;

private:
    DOFOptions m_DOFOptions;
    bool m_inMapLoading;
    ReadRST m_readRST;
    int m_open_err;
    int SetNodeChoices();

    vistle::FloatParameter *m_scale;
    vistle::StringParameter *m_filename;
    vistle::StringParameter *m_solutionChoice, *m_NodeChoices, *m_ElementChoices, *m_SolidShellComponents,
        *m_BeamComponents, *m_AxiShellComponents, *m_TopBottomOpts, *m_ThermalFluxOpts, *m_VolEnergyOpts,
        *m_MagFluxDensOpts;
    vistle::IntParameter *m_vertex_based, *m_outputDecode;
    vistle::Port *m_grid_out, *m_field, *m_materials;
};

#endif
