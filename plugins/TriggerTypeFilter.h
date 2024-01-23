#ifndef __TRIGGERTYPE_FILTER__
#define __TRIGGERTYPE_FILTER__

#include <iostream>
#include "TObjString.h"

#include "interface/PluginBase.h"
#include "interface/TrgTree.h"

class TriggerTypeFilter: public PluginBase
{
public:
    //---ctors---
    TriggerTypeFilter(){};

    //---dtor---
    ~TriggerTypeFilter(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    bool             filterEvents_;
    string           filterName_;
    map<int, string> maskToName_;
    unsigned int     triggerBoard_;
    TObjString       trg_;
    TrgTree          trgTree_;
};

DEFINE_PLUGIN(TriggerTypeFilter);

#endif
