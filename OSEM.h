//
// Created by svf on 2020/9/14.
//

#ifndef C3POSEM_OSEM_H
#define C3POSEM_OSEM_H

#include "plugin.h"

class OSEM : public AbstractPlugin{
public:
    OSEM() = default;
    virtual ~OSEM() = default;

    bool start(const char* path) override;

    bool stop() override;

    double progress() override;
};


#endif //C3POSEM_OSEM_H
