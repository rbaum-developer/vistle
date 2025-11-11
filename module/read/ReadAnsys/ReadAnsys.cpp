#include "ReadAnsys.h"
#include "ReadRST.h"
#include "ANSYS.h"
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

int ReadAnsys::onlyGeometry(UnstructuredGrid::ptr entityGrid)
{
    ANSYS &elem_db_ = ANSYS::get_handle();

    std::vector<int> dummy;
    std::string filename = m_filename->getValue();
    int problems = m_readRST->Read(filename, 1, dummy);
    if (problems) {
        std::cout << "Problems occurred while reading the file " << filename << std::endl;
        return problems;
    }

    // now use nodeindex_, elemindex_, ety_, node_, element_
    int numVertices = 0;
    int elem;
    std::vector<int> e_l;
    std::vector<int> v_ansys_l;
    std::vector<int> t_l;
    int num_supp_elems = 0;

    const auto numElems = m_readRST->getNumElement();
    std::cout << "Number of elements in RST file: " << numElems << std::endl;
    const auto elements = m_readRST->getElements();
    std::cout << "Number of element types in RST file: " << m_readRST->getNumETypes() << std::endl;
    const auto etypes = m_readRST->getETypes();
    const auto numEtypes = m_readRST->getNumETypes();

    for (elem = 0; elem < m_readRST->getNumElement(); ++elem) {
        if (!m_readRST) {
            std::cerr << "ERROR: m_readRST is null\n";
            return -1;
        }
        if (!elements) {
            std::cerr << "ERROR: elements array is null\n";
            return -1;
        }
        int elemIndex = elements[elem].type_; //TODO: invalid element index on position 0
        std::cout << "Element " << elem << " with index " << elemIndex << std::endl;
        if (!etypes || elemIndex < 0 || static_cast<size_t>(elemIndex) >= numEtypes) {
            std::cerr << "ERROR: invalid element-type index " << elemIndex << " (numEtypes=" << numEtypes
                      << "), skipping element " << elem << "\n";
            continue;
        }
        const EType *etype = &etypes[elemIndex];
        std::cout << "Element type: " << etype->id_ << std::endl;
        std::cout << "  routine: " << etype->routine_ << std::endl;
        int routine = etype->routine_;

        int noNodes = getNumberOfNodes(elem, routine);

        if (noNodes <= 0)
            continue; // non-supported element
        ++num_supp_elems;

        t_l.push_back(elem_db_.ElementType(routine, noNodes));
        e_l.push_back(numVertices);

        int vert;

        if (elem_db_.getVistleType(routine) == ANSYS::TYPE_4_NODE_PLANE ||
            elem_db_.getVistleType(routine) == ANSYS::TYPE_8_NODE_PLANE) {
            switch (noNodes) {
            case 3:
                for (vert = 0; vert < 3; ++vert)
                    v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
                break;
            case 4:
                for (vert = 0; vert < 4; ++vert)
                    v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
                break;
            }
        } else if (elem_db_.getVistleType(routine) == ANSYS::TYPE_10_NODE_SOLID) {
            for (vert = 0; vert < noNodes; ++vert)
                v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
        } else if (elem_db_.getVistleType(routine) == ANSYS::TYPE_20_NODE_SOLID) {
            switch (noNodes) {
            case 4:
                for (vert = 0; vert < 4; ++vert) {
                    if (vert != 3)
                        v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
                    else
                        v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert + 1]);
                }
                break;
            case 5:
                for (vert = 0; vert < 5; ++vert)
                    v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
                break;
            case 6:
                for (vert = 0; vert < 7; ++vert) {
                    if (vert != 3)
                        v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
                }
                break;
            case 8:
                for (vert = 0; vert < 8; ++vert)
                    v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
                break;
            }
        } else {
            for (vert = 0; vert < noNodes; ++vert)
                v_ansys_l.push_back(m_readRST->getElements()[elem].nodes_[vert]);
        }
        numVertices += noNodes;
    }

    // decode nodes
    std::vector<int> nodeCodes;
    std::vector<float> x_l, y_l, z_l;
    for (int node = 0; node < m_readRST->getNumNodes(); ++node) {
        nodeCodes.push_back(int(m_readRST->getNodes()[node].id_));
        x_l.push_back(float(m_readRST->getNodes()[node].x_));
        y_l.push_back(float(m_readRST->getNodes()[node].y_));
        z_l.push_back(float(m_readRST->getNodes()[node].z_));
    }
    Map1D nodeDecode(m_readRST->getNumNodes(), &nodeCodes[0]);

    std::vector<int> v_l;
    for (size_t vi = 0; vi < v_ansys_l.size(); ++vi) {
        v_l.push_back(nodeDecode[v_ansys_l[vi]]);
    }

    // Create Vistle UnstructuredGrid using modern pattern
    const Index numElements = (Index)e_l.size();
    const Index numCorners = (Index)v_l.size();
    const Index numVerts = (Index)x_l.size();

    //UnstructuredGrid::ptr entityGrid = std::make_shared<UnstructuredGrid>(numElements, numCorners, numVerts);

    // Copy vertex coordinates
    auto x_coords = entityGrid->x().data();
    auto y_coords = entityGrid->y().data();
    auto z_coords = entityGrid->z().data();
    for (Index i = 0; i < numVerts; ++i) {
        x_coords[i] = x_l[i];
        y_coords[i] = y_l[i];
        z_coords[i] = z_l[i];
    }

    // Copy element list
    auto el = entityGrid->el().data();
    for (Index i = 0; i < numElements; ++i) {
        el[i] = e_l[i];
    }

    // Copy element types
    auto tl = entityGrid->tl().data();
    for (Index i = 0; i < numElements; ++i) {
        tl[i] = t_l[i];
    }

    // Copy connectivity list
    auto cl = entityGrid->cl().data();
    for (Index i = 0; i < numCorners; ++i) {
        cl[i] = v_l[i];
    }

    // make grid:
    /* grid = MakeGridAndObjects(&e_l, std::vector<int> &v_l,
                                                    std::vector<float> &x_l, std::vector<float> &y_l,
                                                    std::vector<float> &z_l, std::vector<int> &t_l,
                                                    const float *const *field,
                                                    int ftype, // 0 = scalar, 1 = vector
                                                    const int *materials);;
  */
    // Note: updateMeta and addObject are expected to be handled by the caller
    return 0;
}

int ReadAnsys::getNumberOfNodes(int elem, int routine)
{
    // Count unique, non-zero node ids for the element. ANSYS elements can have
    // up to 20 nodes (20-node solids), so iterate a safe upper bound.
    const int MAX_ANSYS_NODES = 20;
    int noNodes = 0;
    const int *nodes = m_readRST->getElements()[elem].nodes_;
    std::vector<int> seen;
    for (int i = 0; i < MAX_ANSYS_NODES; ++i) {
        int nid = nodes[i];
        if (nid == 0)
            break;
        if (std::find(seen.begin(), seen.end(), nid) == seen.end())
            seen.push_back(nid);
    }
    noNodes = static_cast<int>(seen.size());
    /* 
    switch (elem_db_.getVistleType(routine)) {
    case ANSYS::TYPE_TARGET:
    case ANSYS::TYPE_TARGET_2D:
        // contact/target elements: skip for now
        noNodes = 0;
        break;

    case ANSYS::TYPE_4_NODE_PLANE:
    case ANSYS::TYPE_8_NODE_PLANE:
    {
        // Count unique, non-zero node entries for the element
        const int maxnodes = ANSYS::TYPE_8_NODE_PLANE + 8; // conservative upper bound
        const auto &nodes = m_readRST->getElements()[elem].nodes_;
        std::vector<int> seen;
        for (int i = 0; i < maxnodes; ++i) {
            int nid = nodes[i];
            if (nid == 0)
                break;
            if (std::find(seen.begin(), seen.end(), nid) == seen.end())
                seen.push_back(nid);
        }
        noNodes = (int)seen.size();
        break;
    }

    case ANSYS::TYPE_10_NODE_SOLID:
        // represented as tetrahedral
        noNodes = 4;
        break;

    case ANSYS::TYPE_20_NODE_SOLID:
    {
        // Count unique node ids (handles degenerate solids)
        const int maxnodes20 = ANSYS::TYPE_20_NODE_SOLID;
        const auto &nodes20 = m_readRST->getElements()[elem].nodes_;
        std::vector<int> seen20;
        for (int i = 0; i < maxnodes20; ++i) {
            int nid = nodes20[i];
            if (nid == 0)
                break;
            if (std::find(seen20.begin(), seen20.end(), nid) == seen20.end())
                seen20.push_back(nid);
        }
        noNodes = (int)seen20.size();
        break;
    }

    case ANSYS::TYPE_HEXAEDER:
    {
        // Handle hexahedra and degenerate variants by counting unique non-zero nodes
        int maxnodes = UnstructuredGrid_Num_Nodes[elem_db_.getVistleType(routine)];
        const auto &hnodes = m_readRST->getElements()[elem].nodes_;
        std::vector<int> seen;
        for (int i = 0; i < maxnodes; ++i) {
            int nid = hnodes[i];
            if (nid == 0)
                break;
            if (std::find(seen.begin(), seen.end(), nid) == seen.end())
                seen.push_back(nid);
        }
        noNodes = (int)seen.size();
        break;
    }

    default:
        // Fallback: use the standard number for the covariant element type
        noNodes = UnstructuredGrid_Num_Nodes[elem_db_.getVistleType(routine)];
        break;
    }
 */
    return noNodes;
}

ReadAnsys::ReadAnsys(const std::string &name, int moduleID, mpi::communicator comm)
: Reader(name, moduleID, comm), m_readRST(new ReadRST())
{
    m_filename = addStringParameter("filename", "name of .rst file", "", Parameter::ExistingFilename);
    setParameterFilters(m_filename, "ANSYS Result Files (*.rst *.rstp *.rth *.rmg *.lnn)");

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

    //TODO: Fix that solution choice is onlyGeometry on default
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

    setParallelizationMode(Serial);
    observeParameter(m_filename);

    new ANSYS(); // ensure ANSYS element database is initialized
}

ReadAnsys::~ReadAnsys()
{
    delete m_readRST;
}

bool ReadAnsys::examine(const vistle::Parameter *param)
{
    if (!param || param == m_filename) {
        const std::string filename = m_filename->getValue();
        if (isCollectionFile(filename)) {
            return true;
        } else {
            sendError("File %s is not of a supported kind, none of \".rst\", \".rfl\", \".rth\", \".rmg\"",
                      filename.c_str());
        }
    } else {
        sendError("File %s does not exist or is not valid", m_filename->getValue().c_str());
        //m_files.clear();
        return false;
    }
    //std::cout << "ReadAnsys::examine returning true with file: " << m_filename->getValue() << std::endl;
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

    if (m_solutionChoice->getValue() == "NodeData") // Nodal solution
    {
        m_open_err = SetNodeChoices();
        if (m_open_err != 0) {
            sendWarning("Problems when setting node choices");
            return false;
        }
    }

    std::cout << "ReadAnsys in read method: " << filename << std::endl;
    std::string solutionChoice = m_solutionChoice->getValue();
    std::cout << "Solution choice: " << solutionChoice << std::endl;
    // only geometry, nodal or element data?
    UnstructuredGrid::ptr grid;
    int problems = 0;
    if (solutionChoice == "OnlyGeometry") {
        try {
            problems = onlyGeometry(grid);
            updateMeta(grid);
            addObject(m_grid_out, grid);
        } catch (const std::exception &e) {
            sendError("Problems occurred while reading grid information of the file");
        }
    }
    /*
    else if (solutionChoice == "NodeData") {
        problems = nodalData();
    }
    else if (solutionChoice == "ElementData") {
        problems = derivedData();
    }
    */
    // readRST_.Reset(ReadRST::RADIKAL); // TODO: Fix when ReadRST is properly integrated
    if (problems) {
        sendError("Problems occurred while reading the file");
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
    m_open_err = m_readRST->OpenFile(filename);
    // loop over rstheader_.numsets_ and accumulate
    // in variable_codes non repeated codes
    int numset;
    for (numset = 1; numset <= m_readRST->getNumTimeSteps(); ++numset) {
        int ErrorReadSHDR = m_readRST->ReadSHDR(numset);
        if (ErrorReadSHDR) {
            sendError("SetNodeChoices: Error in ReadSHDR");
            return ErrorReadSHDR;
        }
        // the important info is in solheader_.numdofs_,
        // solheader_.dof_, solheader_.numexdofs_, solheader_.exdof_
        int new_dof;
        for (new_dof = 0; new_dof < m_readRST->solheader_.numdofs_; ++new_dof) {
            size_t old_dof;
            int repeated = 0;
            for (old_dof = 0; old_dof < variable_codes.size(); ++old_dof) {
                if (m_readRST->solheader_.dof_[new_dof] == variable_codes[old_dof]) {
                    // dof is repeated
                    repeated = 1;
                    break;
                }
            }
            if (!repeated) {
                variable_codes.push_back(m_readRST->solheader_.dof_[new_dof]);
            }
        }
        // now add extra dofs
        int extra_new_dof;
        for (extra_new_dof = 0; extra_new_dof < m_readRST->solheader_.numexdofs_; ++extra_new_dof) {
            size_t old_dof;
            int repeated = 0;
            for (old_dof = 0; old_dof < variable_codes.size(); ++old_dof) {
                if (m_readRST->solheader_.exdof_[extra_new_dof] == variable_codes[old_dof] - EX_OFFSET) {
                    // dof is repeated
                    repeated = 1;
                    break;
                }
            }
            if (!repeated) {
                variable_codes.push_back(m_readRST->solheader_.exdof_[extra_new_dof] + EX_OFFSET);
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

// Adapted MakeGridAndObjects for Vistle
// builds an UnstructuredGrid and a Data Vec (scalar or vector) from ANSYS-style arrays
UnstructuredGrid::ptr ReadAnsys::MakeGridAndObjects(std::vector<int> &e_l, std::vector<int> &v_l,
                                                    std::vector<float> &x_l, std::vector<float> &y_l,
                                                    std::vector<float> &z_l, std::vector<int> &t_l,
                                                    const float *const *field,
                                                    int ftype, // 0 = scalar, 1 = vector
                                                    const int *materials //, UnstructuredGrid::ptr &grid_out,
                                                    //DataBase::ptr &data_out
)
{
    // basic sizes
    size_t numElements = (size_t)e_l.size();
    size_t numCorners = (size_t)v_l.size();
    size_t numVertices = (size_t)x_l.size();

    if (numElements <= 0 || numVertices == 0) {
        sendError("MakeGridAndObjects: empty mesh data    using namespace vistle;");
        //return;
    }

    // create grid
    //UnstructuredGrid::ptr grid_out = std::make_shared<UnstructuredGrid>(numElements, numCorners, numVertices);
    UnstructuredGrid::ptr grid_out(new UnstructuredGrid(numElements, numCorners, numVertices));


    // fill coordinates
    auto gx = grid_out->x().data();
    auto gy = grid_out->y().data();
    auto gz = grid_out->z().data();
    for (size_t i = 0; i < numVertices; ++i) {
        gx[i] = x_l[i];
        gy[i] = y_l[i];
        gz[i] = z_l[i];
    }

    // fill element lists
    auto gel = grid_out->el().data();
    auto gcl = grid_out->cl().data();
    auto gtl = grid_out->tl().data();

    for (size_t i = 0; i < numElements; ++i) {
        gel[i] = e_l[i];
        gtl[i] = t_l[i];
    }
    // copy corner/vertex indices
    for (size_t i = 0; i < numCorners; ++i) {
        gcl[i] = v_l[i];
    }
    // final element entry: total number of corners
    if (numElements > 0)
        gel[numElements] = (int)numCorners;

    return grid_out;
    /* 
    // create data if provided
    data_out.reset();
    if (field && field[0]) {
        if (ftype == 0) { // scalar per-vertex
            Vec<Scalar>::ptr s = std::make_shared<Vec<Scalar>>(numVertices);
            auto sp = s->x().data();
            for (index_t i = 0; i < numVertices; ++i) {
                if (field[0][i] != ReadRST::FImpossible_)
                    sp[i] = (Scalar)field[0][i];
                else
                    sp[i] = (Scalar)0.0;
            }
            s->setGrid(grid_out);
            s->setMapping(DataBase::Mapping::Vertex);
            data_out = s;Run a quick build/compile to verify
        } else if (ftype == 1) { // vector per-vertex
            Vec<Scalar, 3>::ptr v = std::make_shared<Vec<Scalar, 3>>(numVertices);
            auto vx = v->x().data();
            auto vy = v->y().data();
            auto vz = v->z().data();
            for (index_t i = 0; i < numVertices; ++i) {/home/rosba/covise/src/module/general/ReadANSYS/ReadANSYS1.cpp
                if (field[0][i] != ReadRST::FImpossible_ && field[1][i] != ReadRST::FImpossible_ &&
                    field[2][i] != ReadRST::FImpossible_) {
                    vx[i] = (Scalar)field[0][i];
                    vy[i] = (Scalar)field[1][i];
                    vz[i] = (Scalar)field[2][i];
                } else {
                    vx[i] = vy[i] = vz[i] = (Scalar)0.0;
                }
            }
            v->setGrid(grid_out);
            v->setMapping(DataBase::Mapping::Vertex);
            data_out = v;
        }
    } */
}
