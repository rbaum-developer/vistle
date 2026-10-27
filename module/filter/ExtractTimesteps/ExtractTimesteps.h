#ifndef VISTLE_EXTRACTTIMESTEPS_EXTRACTTIMESTEPS_H
#define VISTLE_EXTRACTTIMESTEPS_EXTRACTTIMESTEPS_H

#include <vistle/module/module.h>

#include <mutex>

class ExtractTimesteps: public vistle::Module {
public:
    ExtractTimesteps(const std::string &name, int moduleID, mpi::communicator comm);

private:
    static const unsigned NUMPORTS = 3;
    bool compute(const std::shared_ptr<vistle::BlockTask> &task) const override;
    bool prepare() override;
    bool changeParameter(const vistle::Parameter *param) override;

    struct TimestepSelection {
        bool selected = false;
        bool isTimeSeries = false;
        int outputTimestep = -1;
        int outputNumTimesteps = -1;
    };
    TimestepSelection selectTimesteps(int timestep, int numTimesteps) const;
    vistle::coRestraint m_timesteps;
    mutable std::mutex m_timestepWarningMutex;
    mutable bool m_checkedTimestepSelection = false;

    vistle::StringParameter *p_timesteps = nullptr;
    vistle::Port *p_in[NUMPORTS] = {};
    vistle::Port *p_out[NUMPORTS] = {};
};

#endif
