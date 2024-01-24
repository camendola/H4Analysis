#include "DigitizerReco.h"

//----------Utils-------------------------------------------------------------------------
bool DigitizerReco::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    //---inputs---
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    //---read dV and dT digitizer calibration from file and register it to the shared data
    if(opts.OptExist(instanceName_+".calibration"))
    {
        TFile* calibFile = TFile::Open(opts.GetOpt<string>(instanceName_+".calibration", 0).c_str(), "READ");
        if(calibFile)
        {
            auto calibPtr = calibFile->Get(opts.GetOpt<string>(instanceName_+".calibration", 1).c_str());
            if(calibPtr)
            {
                Log("using calibration "+opts.GetOpt<string>(instanceName_+".calibration", 1)
                    +" from "+opts.GetOpt<string>(instanceName_+".calibration", 0));
                
                digitizerCalib_ = *((DigitizerCalibration*)calibPtr->Clone("dVdTCalibration"));
                RegisterSharedData(&digitizerCalib_, "dVdTCalibration", false);
            }
            else
            {
                Log("calibration "+opts.GetOpt<string>(instanceName_+".calibration", 1)
                    +" not found in "+opts.GetOpt<string>(instanceName_+".calibration", 0)
                    +". No calibration will be applied", WARN);
            }
        }
        else
        {
            Log("impossible to open calibration file "
                +opts.GetOpt<string>(instanceName_+".calibration", 0)
                +". No calibration will be applied", WARN);
        }
    }
    
    //---create channel container (WFClass)
    for(auto& channel : channelsNames_)
    {
        nSamples_[channel] = opts.GetOpt<int>(channel+".nSamples");
        auto tUnit = opts.GetOpt<float>(channel+".tUnit");
        auto polarity = opts.OptExist(channel+".polarity") ? opts.GetOpt<int>(channel+".polarity") : 1;
        if(opts.OptExist(channel+".type"))
        {
            std::cout << " >> channel.type" << std::endl;
            if(opts.GetOpt<string>(channel+".type") == "NINO")
                WFs_[channel] = new WFClassNINO(polarity, tUnit);
            else if(opts.GetOpt<string>(channel+".type") == "Clock")
                WFs_[channel] = new WFClassClock(polarity,tUnit);
            else if(opts.GetOpt<string>(channel+".type") == "LiTeDTU")
                WFs_[channel] = new WFClassLiTeDTU(polarity, tUnit);
        }
        else
            WFs_[channel] = new WFClass(polarity, tUnit);
        
        
        //---set channel calibration if available
        unsigned int digiBd = opts.GetOpt<unsigned int>(channel+".digiBoard");
        unsigned int digiGr = opts.GetOpt<unsigned int>(channel+".digiGroup");
        unsigned int digiCh = opts.GetOpt<unsigned int>(channel+".digiChannel");
        auto ch_key = make_tuple(digiBd, digiGr, digiCh);
        if(digitizerCalib_.find(ch_key) != digitizerCalib_.end())
            WFs_[channel]->SetCalibration(&digitizerCalib_[ch_key]);

        //---register WF in the shared data
        RegisterSharedData(WFs_[channel], channel, false);
    }
    
    return true;
}

bool DigitizerReco::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{        
    //---read the digitizer
    //---set time reference from digitized trigger
    int trigRef=0;
    // FIXME
    // for(int iSample=nSamples_[channel]*8; iSample<nSamples_[channel]*9; ++iSample)
    // {
    //     if(event.digiSampleValue[iSample] < 1000)
    //     {
    //         trigRef = iSample-nSamples_[channel]*8;
    //         break;
    //     }
    // }
    
    //---user channels
    bool evtStatus = true;
    //--- check the time alignement
//    for (int i = 0; i<6; i++){
//        long long int Tdiff = event.evtTime[2*i+1]-event.evtTime[0];
//        //std::cout<< " check desync wrt " << 2*i+1 << " event " << event.evtNumber << " Tdiff " << Tdiff/1000. << std::endl;
//        if(fabs(Tdiff/1000.) > 0.5){
//            evtStatus = false;
//            std::cout<< " found desync wrt " << 2*i+1 << " event " << event.evtNumber << " Tdiff " << Tdiff << std::endl;
//        }
//    }
//
    
    for(auto& channel : channelsNames_)
    {
        //---reset and read new WFs_
        WFs_[channel]->Reset();
	
        auto digiBd = opts.GetOpt<unsigned int>(channel+".digiBoard");
        auto digiGr = opts.GetOpt<unsigned int>(channel+".digiGroup");
        auto digiCh = opts.GetOpt<unsigned int>(channel+".digiChannel");
        auto offset = event.digiMap.at(make_tuple(digiBd, digiGr, digiCh));
        auto max_sample = offset+std::min(nSamples_[channel], event.digiNSamplesMap[make_tuple(digiBd, digiGr, digiCh)]); 
        auto iSample = offset;
        while(iSample < max_sample && event.digiBoard[iSample] != -1)
        {
            //Set the start index cell
            if (iSample==offset)
                WFs_[channel]->SetStartIndexCell(event.digiStartIndexCell[iSample]);

            //---H4DAQ bug: sometimes ADC value is out of bound.
            //---skip everything if one channel is bad
            if(event.digiSampleValue[iSample] > 1e6)
            {
                evtStatus = false;
                WFs_[channel]->AddSample(4095);
            }
            else if (channel == "CLK")
            {
                WFs_[channel]->AddSample(event.digiSampleValue[iSample]);
            }
            else
            {
                WFs_[channel]->AddSample(event.digiSampleValue[iSample], event.digiSampleGain[iSample]);
            }
            iSample++;
        }

        if(opts.OptExist(channel+".useTrigRef") && opts.GetOpt<bool>(channel+".useTrigRef"))
            WFs_[channel]->SetTrigRef(trigRef);
        if(WFs_[channel]->GetCalibration())
            WFs_[channel]->ApplyCalibration();
    }

    if(!evtStatus)
        Log("bad amplitude detected", WARN);

    return evtStatus;
}
