#include "Sample.h"
#include <vistle/core/object.h>

#ifdef SAMPLEVTKM
#include <viskores/cont/DataSet.h>
#include <viskores/filter/resampling/Probe.h>

#include <vistle/vtkm/convert.h>
#endif

using namespace vistle;

MODULE_MAIN(Sample)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(InterpolationMode,
                                    (First) // value of first vertex
                                    (Mean) // mean value of all vertices
                                    (Nearest) // value of nearest vertex
                                    (Linear) // barycentric/multilinear interpolation
)
DEFINE_ENUM_WITH_STRING_CONVERSIONS(ValueOutside, (NaN)(Zero)(userDefined))
DEFINE_ENUM_WITH_STRING_CONVERSIONS(MultiHits, (Any)(Average))

Sample::Sample(const std::string &name, int moduleID, mpi::communicator comm): Module(name, moduleID, comm)
{
    setSchedulingPolicy(message::SchedulingPolicy::Single);
    setReducePolicy(message::ReducePolicy::PerTimestepOrdered);

    createInputPort("data_in", "object with data to be sampled");
    createInputPort("ref_in", "target grid", Port::Flags::NOCOMPUTE);

    m_out = createOutputPort("data_out", "data sampled to target grid");

    m_mode = addIntParameter("mode", "interpolation mode", Linear, Parameter::Choice);
    m_createCelltree = addIntParameter("create_celltree", "create celltree", 0, Parameter::Boolean);
    m_valOutside =
        addIntParameter("value_outside", "value to be used if target is outside source domain", NaN, Parameter::Choice);
    m_userDef = addFloatParameter("user_defined_value", "user defined value if target outside source domain", 1.0);
    m_hits = addIntParameter("mulit_hits", "handling if multiple interpolatied values found for one target point ",
                             Linear, Parameter::Choice);
    V_ENUM_SET_CHOICES(m_mode, InterpolationMode);
    V_ENUM_SET_CHOICES(m_valOutside, ValueOutside);
    V_ENUM_SET_CHOICES(m_hits, MultiHits);
}

int Sample::SampleToGrid(const vistle::GeometryInterface *target, vistle::DataBase::const_ptr inData,
                         const vistle::Object::const_ptr targetObj, vistle::DataBase::ptr &sampled)
{
    int found = 0;
    auto inObj = inData->grid();
    const float INVALID_VALUE = getInvalidValue();
    const GridInterface *inGrid = inObj->getInterface<GridInterface>();
    if (!inGrid) {
        std::cerr << "Failed to pass grid" << std::endl;
        return found;
    }
    Vec<Scalar>::const_ptr scal = Vec<Scalar>::as(inData);
    if (!scal) {
        std::cerr << "no scalar data received" << std::endl;
        return found;
    }
#ifdef SAMPLEVTKM
    viskores::filter::resampling::Probe filter;

    // Set up target geometry
    viskores::cont::DataSet refVTKMs;
    vistle::vtkmSetGrid(refVTKMs, targetObj);
    filter.SetGeometry(refVTKMs);
    filter.SetInvalidValue(INVALID_VALUE);

    // Set up source data and field
    viskores::cont::DataSet inVTKMs;
    vistle::vtkmSetGrid(inVTKMs, inData->grid());
    std::string fieldName = "data_in";
    vistle::vtkmAddField(inVTKMs, inData, fieldName);
    filter.SetActiveField(fieldName);

    std::string activeFieldName = filter.GetActiveFieldName();
    std::cout << "Rank " << rank() << " executing filter" << std::endl;
    std::cout << "Active field name: " << activeFieldName << std::endl;

    if (inVTKMs.GetNumberOfFields() == 0) {
        std::cerr << "Dataset has no fields" << std::endl;
    }
    // Execute filter
    auto outVTKMs = filter.Execute(inVTKMs);

    // Convert output grid back to Vistle
    vistle::Object::const_ptr outGrid = vistle::vtkmGetGeometry(outVTKMs);
    if (!outGrid) {
        sendError("Failed to extract output grid from VTK-m result");
        return found;
    }

    // Convert output field back to Vistle
    std::string outFieldName = fieldName; // or the filter's output field name
    sampled = vistle::vtkmGetField(outVTKMs, outFieldName);
    if (!sampled) {
        sendError("Failed to extract output field %s from VTK-m result", outFieldName.c_str());
        return found;
    } else {
        sampled->setGrid(outGrid);
        found = 1;
    }

    std::cout << "Output dataset has " << outVTKMs.GetNumberOfFields() << " fields" << std::endl;
    outVTKMs.PrintSummary(std::cout);

    for (int i = 0; i < outVTKMs.GetNumberOfFields(); ++i) {
        auto field = outVTKMs.GetField(i);
        std::cout << "Field " << i << ": " << field.GetName() << std::endl;
    }

//outGrid = vistle::vtkmGetGeometry(outVTKM); // convert back to vistle dataout
//sampled= ?
#else
    const Scalar *data = scal->x().data();

    Index numVert = target->getNumVertices();
    Vec<Scalar>::ptr dataOut(new Vec<Scalar>(numVert));
    Scalar *ptrOnData = dataOut->x().data();

    for (Index i = 0; i < numVert; ++i) {
        Vector3 v = target->getVertex(i);
        Index cellIdxIn =
            inGrid->findCell(v, InvalidIndex, m_useCelltree ? GridInterface::NoFlags : GridInterface::NoCelltree);
        if (cellIdxIn != InvalidIndex) {
            GridInterface::Interpolator interp = inGrid->getInterpolator(cellIdxIn, v, DataBase::Vertex, mode);
            Scalar p = interp(data);
            ptrOnData[i] = p;
            found = 1;
        } else {
            ptrOnData[i] = NO_VALUE;
        }
    }
    sampled = DataBase::as(Object::ptr(dataOut));
#endif

    return found;
}

bool Sample::reduce(int timestep)
{
    mode = (GridInterface::InterpolationMode)m_mode->getValue();

    const float valOut = getInvalidValue();

    if (m_createCelltree->getValue())
        m_useCelltree = true;
    int nProcs = 1;
    if (comm().size() > 0)
        nProcs = comm().size();
    int numLocalObjs = 0;
    numLocalObjs = objListLocal.size();
    std::vector<int> objPerRank;

    mpi::all_gather(comm(), numLocalObjs, objPerRank);

    for (int r = 0; r < nProcs; ++r) {
        if (objPerRank[r] < 1)
            continue;
        for (int n = 0; n < objPerRank[r]; ++n) {
            Object::const_ptr objAtIdx;
            if (r == rank()) {
                objAtIdx = objListLocal.at(n);
            }
            broadcastObject(objAtIdx, r);

            const GeometryInterface *geo = objAtIdx->getInterface<GeometryInterface>();
            DataBase::ptr sampledData;
            std::vector<unsigned> foundList;
            std::vector<Object::const_ptr> sampledDataList;

            unsigned found = 0;
            if (!dataList.empty()) {
                DataBase::const_ptr data = dataList.at(0);

                for (unsigned long dIdx = 0; dIdx < dataList.size(); ++dIdx) {
                    data = dataList.at(dIdx);
                    if (data->getTimestep() != timestep)
                        continue;
                    int foundSample = 0;
                    if (geo)
                        foundSample += SampleToGrid(geo, data, objAtIdx, sampledData);

                    if (foundSample > 0)
                        sampledDataList.push_back(sampledData);
                    found += foundSample;
                }
            }
            std::cout << "Rank " << rank() << " found " << found << " samples for timestep " << timestep << std::endl;

            mpi::gather(comm(), found, foundList, r);

            if (rank() == r) {
                for (unsigned long fr = 0; fr < foundList.size(); ++fr) {
                    unsigned nFound = foundList.at(fr);
                    if ((nFound > 0) && (fr != unsigned(r))) {
                        for (unsigned fIdx = 0; fIdx < nFound; ++fIdx) {
                            auto sampledElem = receiveObject(comm(), fr);
                            sampledDataList.push_back(sampledElem);
                        }
                    }
                }

                Vec<Scalar>::ptr outData(new Vec<Scalar>(geo->getNumVertices()));
                auto globDatVec = outData->x().data();
                if (m_hits->getValue() == Average) {
                    if (!sampledDataList.empty()) {
                        for (Index bIdx = 0; bIdx < geo->getNumVertices(); ++bIdx) {
                            globDatVec[bIdx] = 0;
                        }
                        std::vector<int> numHits(geo->getNumVertices(), 0);
                        for (auto elem: sampledDataList) {
                            auto locDat = Vec<Scalar, 1>::as(elem);
                            auto locDatVec = locDat->x().data();

                            for (Index bIdx = 0; bIdx < locDat->getSize(); ++bIdx) {
                                if (locDatVec[bIdx] != NO_VALUE) {
                                    numHits[bIdx] += 1;
                                    globDatVec[bIdx] += locDatVec[bIdx];
                                }
                            }
                        }
                        for (Index bIdx = 0; bIdx < geo->getNumVertices(); ++bIdx) {
                            if (numHits[bIdx] > 0)
                                globDatVec[bIdx] /= numHits[bIdx];
                            else
                                globDatVec[bIdx] = valOut;
                        }
                        numHits.clear();
                    } else {
                        for (Index bIdx = 0; bIdx < geo->getNumVertices(); ++bIdx) {
                            globDatVec[bIdx] = valOut;
                        }
                    }
                } else {
                    for (Index bIdx = 0; bIdx < geo->getNumVertices(); ++bIdx) {
                        globDatVec[bIdx] = valOut;
                    }

                    if (!sampledDataList.empty()) {
                        for (auto elem: sampledDataList) {
                            auto locDat = Vec<Scalar, 1>::as(elem);
                            auto locDatVec = locDat->x().data();

                            for (Index bIdx = 0; bIdx < locDat->getSize(); ++bIdx) {
                                if (locDatVec[bIdx] != NO_VALUE) {
                                    globDatVec[bIdx] = locDatVec[bIdx];
                                }
                            }
                        }
                    }
                }

                Object::const_ptr outGrid = objAtIdx;
                outData->setTimestep(timestep);
                outData->setBlock(blockIdx.at(n));
                outData->setMapping(DataBase::Vertex);
                outData->setGrid(outGrid);
                outData->describe("scalar", id());
                updateMeta(outData);
                addObject(m_out, outData);
            } else {
                if (found > 0) {
                    for (auto elem: sampledDataList) {
                        sendObject(comm(), elem, r);
                    }
                }
            }
            sampledDataList.clear();
            foundList.clear();
        }
    }

    comm().barrier();
    objPerRank.clear();

    if (dataList.empty() || (timestep == dataList.at(0)->getNumTimesteps() - 1) ||
        (dataList.at(0)->getNumTimesteps() < 2)) {
        dataList.clear();
        blockIdx.clear();
        objListLocal.clear();
    }
    return true;
}

bool Sample::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    #ifndef SAMPLEVTKM
    if (m_createCelltree->getValue()) {
        Object::const_ptr obj = task->expect<Object>("data_in");
        if (!obj)
            return true;
        auto unstr = UnstructuredGrid::as(obj);
        auto str = StructuredGrid::as(obj);

        if (!str && !unstr) {
            DataBase::const_ptr data = DataBase::as(obj);
            if (data && data->grid()) {
                unstr = UnstructuredGrid::as(data->grid());
                str = StructuredGrid::as(data->grid());
            }
        }
        if (str) {
            if (!str->hasCelltree()) {
                str->getCelltree();
            }
        } else if (unstr) {
            if (!unstr->hasCelltree()) {
                unstr->getCelltree();
            }
        } else {
            sendInfo("Source grid type not supported");
        }
    }
    #endif
    return true;
}

bool Sample::objectAdded(int sender, const std::string &senderPort, const Port *port)
{
    auto object = port->objects().back();

    if (port->getName() == "data_in") {
        auto data = DataBase::as(object);
        if (data && data->grid()) {
            dataList.push_back(data);
        }
    } else if (port->getName() == "ref_in") {
        if (object->getTimestep() < 1) {
            const GeometryInterface *geo = object->getInterface<GeometryInterface>();
            if (geo) {
                objListLocal.push_back(object);
                blockIdx.push_back(object->getBlock());
            } else {
                DataBase::const_ptr data = DataBase::as(object);
                if (data && data->grid()) {
                    objListLocal.push_back(data->grid());
                    blockIdx.push_back(object->getBlock());
                }
            }
        }
    }
    return true;
}

bool Sample::changeParameter(const Parameter *p)
{
    return true;
}

float Sample::getInvalidValue() const
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
