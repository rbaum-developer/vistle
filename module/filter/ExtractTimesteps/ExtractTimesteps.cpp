#include <vistle/alg/objalg.h>
#include <vistle/core/database.h>
#include <vistle/core/object.h>
#include <vistle/util/coRestraint.h>

#include "ExtractTimesteps.h"

#include <algorithm>

using namespace vistle;

ExtractTimesteps::ExtractTimesteps(const std::string &name, int moduleID, mpi::communicator comm)
: Module(name, moduleID, comm)
{
    for (unsigned i = 0; i < NUMPORTS; ++i) {
        p_in[i] = createInputPort(i == 0 ? "data_in" : "data_in" + std::to_string(i), "input data");
        p_out[i] = createOutputPort(i == 0 ? "data_out" : "data_out" + std::to_string(i), "filtered data");
        linkPorts(p_in[i], p_out[i]);
        if (i > 0)
            setPortOptional(p_in[i], true);
    }

    p_timesteps =
        addStringParameter("timesteps", "timestep indices to keep, e.g. single timesteps 0,2,4 or ranges 1-5 or every nth timestep in a range 1-10/2", "all", Parameter::Restraint);
}

bool ExtractTimesteps::prepare()
{
    m_timesteps.clear();
    m_timesteps.add(p_timesteps->getValue());
    {
        std::lock_guard<std::mutex> lock(m_timestepWarningMutex);
        m_checkedTimestepSelection = false;
    }
    return true;
}

bool ExtractTimesteps::compute(const std::shared_ptr<BlockTask> &task) const
{
    Object::const_ptr inputs[NUMPORTS];
    DataComponents splits[NUMPORTS];
    Object::const_ptr grid;
    int first = -1;

    for (unsigned i = 0; i < NUMPORTS; ++i) {
        if (!p_out[i]->isConnected()){
            continue;
        }
        if (!p_in[i]->isConnected()) {
            sendError("require input on %s for %s", p_in[i]->getName().c_str(), p_out[i]->getName().c_str());
            return true;
        }

        inputs[i] = task->expect<Object>(p_in[i]);
        if (!inputs[i]) {
            sendError("no input on %s for %s", p_in[i]->getName().c_str(), p_out[i]->getName().c_str());
            return true;
        }
        splits[i] = splitContainerObject(inputs[i]);
        if (first < 0)
            first = i;

        if (splits[i].geometry) {
            if (!grid) {
                grid = splits[i].geometry;
            } else if (grid->getHandle() != splits[i].geometry->getHandle()) {
                sendError("input ports have different grids");
                return true;
            }
        }
    }

    if (!grid) {
        sendError("no grid on any input port");
        return true;
    }

    int numTimesteps = inputs[first]->getNumTimesteps();
    if (numTimesteps <= 0 && splits[first].mapped)
        numTimesteps = splits[first].mapped->getNumTimesteps();
    if (numTimesteps <= 0)
        numTimesteps = grid->getNumTimesteps();

    if (numTimesteps > 0) {
        if (!m_checkedTimestepSelection) {
            std::lock_guard<std::mutex> lock(m_timestepWarningMutex);
            m_checkedTimestepSelection = true;
            const auto &selectedTimesteps = m_timesteps.getValues();
            const auto outOfRange = std::find_if(selectedTimesteps.begin(), selectedTimesteps.end(),
                                                 [numTimesteps](ssize_t selected) {
                                                     return selected >= numTimesteps;
                                                 });
            if (outOfRange != selectedTimesteps.end()) {
                sendWarning("selected timestep %zd is outside input data with %d timesteps", *outOfRange,
                            numTimesteps);
            }
        }
    }

    const auto timeSel = selectTimesteps(splits[first].timestep, numTimesteps);
    if (!timeSel.selected){        
        return true; // skip this timestep
    }
    auto outputGrid = grid->clone();
    outputGrid->setTimestep(timeSel.outputTimestep);
    if (timeSel.isTimeSeries)
        outputGrid->setNumTimesteps(timeSel.outputNumTimesteps);
    updateMeta(outputGrid);

    for (unsigned i = 0; i < NUMPORTS; ++i) {
        if (!p_out[i]->isConnected())
            continue;

        Object::ptr out;
        if (splits[i].mapped) {
            out = splits[i].mapped->clone();
            if (auto data = DataBase::as(out))
                data->setGrid(outputGrid);
            out->setTimestep(timeSel.outputTimestep);
            if (timeSel.isTimeSeries)
                out->setNumTimesteps(timeSel.outputNumTimesteps);
            updateMeta(out);
        } else if (splits[i].geometry) {
            out = outputGrid;
        } else {
            out = inputs[i]->clone();
            out->setTimestep(timeSel.outputTimestep);
            if (timeSel.isTimeSeries)
                out->setNumTimesteps(timeSel.outputNumTimesteps);
            updateMeta(out);
        }
        task->addObject(p_out[i], out);
    }
    return true;
}

ExtractTimesteps::TimestepSelection ExtractTimesteps::selectTimesteps(int timestep, int numTimesteps) const
{
    TimestepSelection timeSel;
    timeSel.outputTimestep = timestep;
    timeSel.outputNumTimesteps = numTimesteps;

    if (timestep >= 0 && !m_timesteps(timestep))
        return timeSel;
    timeSel.selected = true;

    timeSel.isTimeSeries = numTimesteps > 1 || (numTimesteps < 0 && timestep >= 0);
    if (timeSel.isTimeSeries && timestep >= 0) {
        auto countSelected = [&](int begin, int end) {
            int count = 0;
            for (int t = begin; t < end; ++t) {
                if (m_timesteps(t))
                    ++count;
            }
            return count;
        };
        timeSel.outputTimestep = countSelected(0, timestep);
        if (numTimesteps > 0)
            timeSel.outputNumTimesteps = countSelected(0, numTimesteps);
    }
    return timeSel;
}

bool ExtractTimesteps::changeParameter(const Parameter *param)
{
    if (param == p_timesteps) {
        sendInfo("Changed parameter %s to %s", p_timesteps->getName().c_str(), p_timesteps->getValue().c_str());
    }
    return Module::changeParameter(param);
}

MODULE_MAIN(ExtractTimesteps)
