//
// Created by svf on 2020/9/14.
//

#ifndef C3POSEM_PLUGIN_H
#define C3POSEM_PLUGIN_H

#include "abstract_plugin.h"

typedef AbstractPlugin* (*CreatePluginFunc)();

struct Plugin {
    const char* mPluginName;
    const char* mPluginVersion;
    CreatePluginFunc mCreateFunc;
};

#define PLUGIN(classType, pluginName, pluginVersion) \
    extern "C" {                                     \
        __attribute__ ((visibility ("default"))) AbstractPlugin* createPlugin() {             \
            return new classType();                  \
        }                                            \
        __attribute__ ((visibility ("default"))) Plugin exports = {                           \
            pluginName,                              \
            pluginVersion,                           \
            createPlugin,                            \
        };                                           \
    }
#endif //C3POSEM_PLUGIN_H
