#ifndef __WF_ANALYZER__
#define __WF_ANALYZER__

#include <iostream>
#include <typeinfo>

#include "TObjString.h"

#include "interface/PluginBase.h"
#include "interface/DigiTree.h"
#include "interface/WFTree.h"
#include "interface/WFClass.h"
#include "interface/WFClassLiTeDTU.h"
#include "interface/WFClassNINO.h"
#include "interface/WFClassClock.h"
#include "interface/WFViewer.h"
#include "interface/RecoEventAnalyzer.h"

class WFAnalyzer: public PluginBase
{
public:
    //---ctors---
    WFAnalyzer(){};

    //---dtor---
    ~WFAnalyzer(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    string                      srcInstance_;
    string                      trgInstance_;
    string                      trg_;
    vector<string>              channelsNames_;
    vector<string>              timeRecoTypes_;
    map<string, vector<float> > timeOpts_;
    DigiTree                    digiTree_;
    RecoEventAnalyzer           eventAnalyzer_;
    WFTree                      outWFTree_;
    map<string, WFClass*>       WFs_;
    map<string, map<string, TH1*> > templates_;
    map<string, TH1*>           spikeTemplates_;
};

DEFINE_PLUGIN(WFAnalyzer);

#endif
