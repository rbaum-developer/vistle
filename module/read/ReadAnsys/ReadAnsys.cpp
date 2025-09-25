#include "ReadAnsys.h"
#include <boost/algorithm/string/predicate.hpp>
#include <vistle/core/lines.h>
#include <vistle/core/object.h>
#include <vistle/core/points.h>
#include <vistle/core/polygons.h>
#include <vistle/core/rectilineargrid.h>
#include <vistle/core/structuredgrid.h>
#include <vistle/core/uniformgrid.h>
#include <vistle/core/unstr.h>
#include <vistle/core/vec.h>
#include <boost/mpi/communicator.hpp>

#include <vistle/util/filesystem.h>

#include "ReadRST.h"
#include "ReadRST.cpp"
#include "DOFOptions.h"

using namespace vistle;

/* const char *dofname[] = {
    "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ", 
    "AX", "AY", "AZ", "VX", "VY", "VZ", 
    "PRES", "TEMP", "VOLT", "MAG", "ENKE", "ENDS", 
    "EMF", "CURR"
};

const char *exdofname[] = {
    "HDSP", "CONC", "WARP", "TIME", "PHASE", 
    "FLUX", "HEAT", "RADI", "MASS", "VISC"
}; */

const char *dofname[32] = {"UX",      "UY",      "UZ",   "ROTX", "ROTY",    "ROTZ",    "AX",      "AY",
                           "AZ",      "VX",      "VY",   "VZ",   "unused1", "unused2", "unused3", "unused4",
                           "unused5", "unused6", "PRES", "TEMP", "VOLT",    "MAG",     "ENKE",    "ENDS",
                           "EMF",     "CURR",    "SP01", "SP02", "SP03",    "SP04",    "SP05",    "SP06"};

const char *exdofname[28] = {"DENS", "VISC", "EVIS", "COND", "ECON", "LMD1", "LMD2", "LMD3", "LMD4", "LMD5",
                             "LMD6", "EMD1", "EMD2", "EMD3", "EMD4", "EMD5", "EMD6", "PTOT", "TTOT", "PCOE",
                             "MACH", "STRM", "HFLU", "HFLM", "YPLU", "TAUW", "SPHT", "CMUV"};

MODULE_MAIN(ReadAnsys)

bool isCollectionFile(const std::string &fn)
{
    constexpr const char *collectionEndings[] = {".rst", ".rfl", ".rth", ".rmg"};
    for (const auto ending: collectionEndings)
        if (boost::algorithm::ends_with(fn, ending))
            return true;

    return false;
}

ReadAnsys::ReadAnsys(const std::string &name, int moduleID, mpi::communicator comm): Reader(name, moduleID, comm)
{
    m_filename = addStringParameter("filename", "name of .rst file", "", Parameter::ExistingFilename);
    setParameterFilters(m_filename, "ANSYS Result Files (*.rst;*.rfl;*.rth;*.rmg)");

    // Define Parameters:
    m_scale = addFloatParameter("ScaleGridDisplacement", "scale grid displacement", 1.0);

    const std::vector<std::string> NodeChoices = {"none"};
    const std::vector<std::string> ElementChoices = {"none",
                                                     "Stresses",
                                                     "Elastic strains",
                                                     "Plastic strains",
                                                     "Creep strains",
                                                     "Thermal strains",
                                                     "Field fluxes",
                                                     "Volume and energies",
                                                     "Magnetic flux density"};

    const std::vector<std::string> SolidComponents = {"none", "XX", "YY", "ZZ", "XY", "YZ",
                                                      "ZX",   "T1", "T2", "T3", "TI", "TIGE"};
    const std::vector<std::string> BeamComponents = {"none", "Axial", "Yp", "Ym", "Zp", "Zm"};
    const std::vector<std::string> AxiShellComponents = {"none", "Meridional", "ThroughThickness", "Hoop",
                                                         "Meridional-hoop"};
    const std::vector<std::string> TopBottomOpts = {"Top", "Bottom", "Average"};
    const std::vector<std::string> ThermalFluxOpts = {"none", "QX", "QY", "QZ", "Q"};
    const std::vector<std::string> VolEnergyOpts = {"Volume", "SENE", "KENE"};
    const std::vector<std::string> MagFluxDensOpts = {"B", "BX", "BY", "BZ", "BSUM"};
    const std::vector<std::string> SolutionChoices = {"OnlyGeometry", "NoteData", "ElementData"};

    m_solutionChoice = addStringParameter("SolutionChoice", "Please enter your choice", "", Parameter::Choice);
    setParameterChoices(m_solutionChoice, SolutionChoices);
    m_NodeChoices = addStringParameter("DOF_Solution", "Degrees of freedom of the solution", "", Parameter::Choice);
    setParameterChoices(m_NodeChoices, NodeChoices);
    m_ElementChoices = addStringParameter("Derived_Solution", "Derived variables", "", Parameter::Choice);
    setParameterChoices(m_ElementChoices, ElementChoices);
    m_SolidShellComponents = addStringParameter("SolidShellComponents", "Stress components", "", Parameter::Choice);
    setParameterChoices(m_SolidShellComponents, SolidComponents);
    m_BeamComponents = addStringParameter("BeamComponents", "Beam stress components", "", Parameter::Choice);
    setParameterChoices(m_BeamComponents, BeamComponents);
    m_AxiShellComponents =
        addStringParameter("AxiShellComponents", "Axisymmetric-shell stress components", "", Parameter::Choice);
    setParameterChoices(m_AxiShellComponents, AxiShellComponents);
    m_TopBottomOpts = addStringParameter("TopBottomOpts", "Top/bottom/average options", "", Parameter::Choice);
    setParameterChoices(m_TopBottomOpts, TopBottomOpts);
    m_ThermalFluxOpts = addStringParameter("ThermalFlux", "Thermal flux", "", Parameter::Choice);
    setParameterChoices(m_ThermalFluxOpts, ThermalFluxOpts);
    m_VolEnergyOpts = addStringParameter("VolEnergyOpts", "Volume and energy", "", Parameter::Choice);
    setParameterChoices(m_VolEnergyOpts, VolEnergyOpts);
    m_MagFluxDensOpts = addStringParameter("MagFluxDensity", "Magnetic flux density", "", Parameter::Choice);
    setParameterChoices(m_MagFluxDensOpts, MagFluxDensOpts);

    m_vertex_based = addIntParameter("AlwaysVertexBased", "vertex based", true, Parameter::Boolean);
    m_outputDecode = addIntParameter("OutputNodeDecode", "output node decode", false, Parameter::Boolean);

    // Define Output Ports:
    m_grid_out = createOutputPort("grid", "grid data"); // outputs unstructured grid
    m_field = createOutputPort("data", "output data");
    m_materials = createOutputPort("materials", "output material labels");

    setParallelizationMode(ParallelizeTimeAndBlocks);
    observeParameter(m_filename);
}

ReadAnsys::~ReadAnsys()
{}

bool ReadAnsys::examine(const vistle::Parameter *param)
{
    if (param == nullptr || param == m_filename) {
        sendError("File %s does not exist or is not valid", m_filename->getValue().c_str());
        //m_files.clear();
        return false;
    }

    const std::string filename = m_filename->getValue();
    if (isCollectionFile(filename)) {
        return true;
    } else {
        sendError("File %s is not of a supported kind, none of \".rst\", \".rfl\", \".rth\", \".rmg\"",
                  filename.c_str());
    }
    return true;
}

bool ReadAnsys::prepareRead()
{
    return true;
}

bool ReadAnsys::finishRead()
{
    return true;
}

bool ReadAnsys::read(Token &token, int timestep, int block)
{
    const std::string filename = m_filename->getValue();
    if (filename.empty()) {
        sendError("No filename specified");
        return false;
    }


    return true;
}

int compar(const void *l, const void *r)
{
    int *left = (int *)l;
    int *right = (int *)r;
    if (*left < *right)
        return -1;
    if (*left == *right)
        return 0;
    return 1;
}

int ReadAnsys::SetNodeChoices()
{
    std::vector<std::string> Choices;
    std::vector<int> variable_codes;
    std::string filename = m_filename->getValue();
    m_open_err = m_readRST.OpenFile(filename);
    // loop over rstheader_.numsets_ and accumulate
    // in variable_codes non repeated codes
    int numset;
    for (numset = 1; numset <= m_readRST.getNumTimeSteps(); ++numset) {
        int ErrorReadSHDR = m_readRST.ReadSHDR(numset);
        if (ErrorReadSHDR) {
            sendError("SetNodeChoices: Error in ReadSHDR");
            return ErrorReadSHDR;
        }
        // the important info is in solheader_.numdofs_,
        // solheader_.dof_, solheader_.numexdofs_, solheader_.exdof_
        int new_dof;
        for (new_dof = 0; new_dof < m_readRST.solheader_.numdofs_; ++new_dof) {
            size_t old_dof;
            int repeated = 0;
            for (old_dof = 0; old_dof < variable_codes.size(); ++old_dof) {
                if (m_readRST.solheader_.dof_[new_dof] == variable_codes[old_dof]) {
                    // dof is repeated
                    repeated = 1;
                    break;
                }
            }
            if (!repeated) {
                variable_codes.push_back(m_readRST.solheader_.dof_[new_dof]);
            }
        }
        // now add extra dofs
        int extra_new_dof;
        for (extra_new_dof = 0; extra_new_dof < m_readRST.solheader_.numexdofs_; ++extra_new_dof) {
            size_t old_dof;
            int repeated = 0;
            for (old_dof = 0; old_dof < variable_codes.size(); ++old_dof) {
                if (m_readRST.solheader_.exdof_[extra_new_dof] == variable_codes[old_dof] - EX_OFFSET) {
                    // dof is repeated
                    repeated = 1;
                    break;
                }
            }
            if (!repeated) {
                variable_codes.push_back(m_readRST.solheader_.exdof_[extra_new_dof] + EX_OFFSET);
            }
        }
    }
    // now order variable_codes;
    int *var_codes = &variable_codes[0];
    qsort(var_codes, variable_codes.size(), sizeof(int), compar);

    std::vector<int> my_variable_codes;
    size_t code;
    for (code = 0; code < variable_codes.size(); ++code) {
        my_variable_codes.push_back(variable_codes[code]);
        size_t my_size = my_variable_codes.size();
        // extra dofs are always scalar fields
        switch (variable_codes[code]) {
        case 3:
            if (my_size >= 3 && my_variable_codes[my_size - 2] == 2 && my_variable_codes[my_size - 3] == 1) {
                my_variable_codes.push_back(variable_codes[code] + V_OFFSET);
            }
            break;
        case 6:
            if (my_size >= 3 && my_variable_codes[my_size - 2] == 5 && my_variable_codes[my_size - 3] == 4) {
                my_variable_codes.push_back(variable_codes[code] + V_OFFSET);
            }
            break;
        case 9:
            if (my_size >= 3 && my_variable_codes[my_size - 2] == 8 && my_variable_codes[my_size - 3] == 7) {
                my_variable_codes.push_back(variable_codes[code] + V_OFFSET);
            }
            break;
        case 12:
            if (my_size >= 3 && my_variable_codes[my_size - 2] == 11 && my_variable_codes[my_size - 3] == 10) {
                my_variable_codes.push_back(variable_codes[code] + V_OFFSET);
            }
            break;
        }
    }

    // Initialize Choices vector with "none" as first element
    Choices.reserve(1 + my_variable_codes.size());
    Choices.push_back("none");

    size_t choice;
    for (choice = 0; choice < my_variable_codes.size(); ++choice) {
        if (my_variable_codes[choice] < V_OFFSET) {
            Choices.push_back(std::string(dofname[my_variable_codes[choice] - 1]));
        } else if (my_variable_codes[choice] < EX_OFFSET) {
            int substract = my_variable_codes[choice] - V_OFFSET;
            switch (substract) {
            case 3:
                Choices.push_back("U");
                break;
            case 6:
                Choices.push_back("ROT");
                break;
            case 9:
                Choices.push_back("A");
                break;
            case 12:
                Choices.push_back("V");
                break;
            }
        } else {
            int substract = my_variable_codes[choice] - EX_OFFSET;
            Choices.push_back(std::string(exdofname[substract - 1]));
        }
    }

    if (!m_inMapLoading && !m_filename) {
        setParameterChoices(m_NodeChoices, Choices);
    }

    // set the correct state for m_DOFOptions
    m_DOFOptions.num_options_ = (int)my_variable_codes.size() + 1;
    delete[] m_DOFOptions.options_;
    delete[] m_DOFOptions.codes_;
    m_DOFOptions.options_ = new std::string[m_DOFOptions.num_options_];
    m_DOFOptions.codes_ = new std::vector<int>[m_DOFOptions.num_options_];
    m_DOFOptions.options_[0] = "none";
    m_DOFOptions.codes_[0].push_back(0);
    for (choice = 0; choice < my_variable_codes.size(); ++choice) {
        m_DOFOptions.options_[choice + 1] = Choices[choice + 1];
        if (my_variable_codes[choice] < V_OFFSET) {
            m_DOFOptions.codes_[choice + 1].push_back(my_variable_codes[choice]);
        } else if (my_variable_codes[choice] < EX_OFFSET) {
            int substract = my_variable_codes[choice] - V_OFFSET;
            m_DOFOptions.codes_[choice + 1].push_back(substract - 2);
            m_DOFOptions.codes_[choice + 1].push_back(substract - 1);
            m_DOFOptions.codes_[choice + 1].push_back(substract - 0);
        } else {
            m_DOFOptions.codes_[choice + 1].push_back(my_variable_codes[choice]);
        }
    }

    return 0;
}
