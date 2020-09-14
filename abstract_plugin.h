//
// Created by svf on 2020/9/14.
//

#ifndef C3POSEM_ABSTRACT_PLUGIN_H
#define C3POSEM_ABSTRACT_PLUGIN_H

class AbstractPlugin {
public:
    AbstractPlugin() = default;
    virtual ~AbstractPlugin() = default;
    virtual bool start(const char*) = 0;
    virtual bool stop() = 0;
    virtual double progress() = 0;
};

#endif //C3POSEM_ABSTRACT_PLUGIN_H
