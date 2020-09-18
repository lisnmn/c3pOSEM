//
// Created by svf on 2020/9/14.
//

#ifndef C3POSEM_OSEM_H
#define C3POSEM_OSEM_H

#include <thread>
#include <memory>

#include "plugin.h"

#include "inifile.h"
#include "UReconPETOSEMAlgor.h"
#include "UProgressNotify.h"
#include "UPETScanner.h"

class OSEM : public AbstractPlugin{
public:
    OSEM() = default;
    virtual ~OSEM() = default;

    bool start(const char* path) override;

    bool stop() override;

    double progress() override;


private:
    inifile::IniFile mIniFile;
    UPETScanner mScanner;
    UReconPETOSEMPara mOSEMPara;
    UProgressNotify mNotify;
    UReconPETOSEMRayTracingAlgor mOSEM;
    std::shared_ptr<std::thread> m_pThreadOSEM;

    void initScanner();
};

PLUGIN(OSEM, "OSEM", "0.1")

#endif //C3POSEM_OSEM_H
