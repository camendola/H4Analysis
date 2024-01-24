#ifndef __DIGITIZER_RECO__
#define __DIGITIZER_RECO__

#include <iostream>
#include <algorithm>

#include "interface/PluginBase.h"
#include "interface/DigiTree.h"
#include "interface/WFTree.h"
#include "interface/WFClass.h"
#include "interface/WFClassLiTeDTU.h"
#include "interface/WFClassNINO.h"
#include "interface/WFClassClock.h"
#include "interface/WFViewer.h"

/**
   The DigitizerReco plugin loads the samples from the digiSample branches
   and creates the WFClass accordingly to the options specified in the config file.

   - Global config file options:
     + `pluginType`: set to `DigitizerReco`
     + `channelsNames`: the list of channel to be read.
     + `calibration`: calibration file (the ROOT file should contain a DigitizerCalibration
       object named "dVdTCalibration").

   - Per-channel options:
     + `type`: WFClass type. Available classes:
       1. WFClass
       2. WFClassClock
       3. WFClassNINO
     + `nSamples` [optional]: maximum number of sample to be loaded. If `nSamples` is larger than the
       maximum number of samples available in raw data (nMax) then `nSamples = nMax`.
     + `tUnit`: time unit (i.e. sampling period).
     + `polarity` [optional, default=1]: pulse polarity. Possible values: `+1`, `-1`. Each sample amplitude is
       multiplied by the polarity after baseline subtraction.
     + `digiBoard`
     + `digiGroup`
     + `digiChannel`
 */
class DigitizerReco: public PluginBase
{
public:
    //---ctors---
    DigitizerReco(){};

    //---dtor---
    ~DigitizerReco(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    map<string, int>            nSamples_;
    vector<string>              channelsNames_;
    map<string, WFClass*>       WFs_;
    DigitizerCalibration        digitizerCalib_; 
};

DEFINE_PLUGIN(DigitizerReco);

#endif
