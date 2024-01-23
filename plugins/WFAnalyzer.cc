#include "WFAnalyzer.h"

//----------Utils-------------------------------------------------------------------------
bool WFAnalyzer::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    trg_ = "";

    //---inputs---
    for(auto& src : {"srcInstanceName", "trgInstanceName"})
    {
        if(!opts.OptExist(instanceName_+"."+src))
        {
            Log("no "+string(src)+" plugin specified", ERR);
            return false;
        }
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    if(opts.OptExist(instanceName_+".trgInstanceName"))
      trgInstance_ = opts.GetOpt<string>(instanceName_+".trgInstanceName");
    else
      trgInstance_="";
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    timeRecoTypes_ = opts.GetOpt<vector<string> >(instanceName_+".timeRecoTypes");

    //---load WFs from source instance shared data
    for(auto& channel : channelsNames_)
    {
        std::cout << " >> channel " << channel << std::endl;
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
        {
	    WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        }
        else
            Log("channels samples not found check DigiReco step", WARN); 
    }
    
    //---channels setup
    string templateTag="";
    if(opts.OptExist(instanceName_+".templateTags"))
        for(auto& tag : opts.GetOpt<vector<string> >(instanceName_+".templateTags"))
            for(auto& run : opts.GetOpt<vector<string> >(tag+".runList"))
                if(run == opts.GetOpt<string>("h4reco.run"))
                    templateTag = tag;

    int nSamples = 0;
    for(auto& channel : channelsNames_)
    {
        nSamples += opts.GetOpt<int>(channel+".nSamples");
        if(opts.OptExist(channel+".templateFit.file"))
        {            
            TFile* templateFile = TFile::Open(opts.GetOpt<string>(channel+".templateFit.file", 0).c_str(), ".READ");
            if(templateFile)
            {
                TH1* wfTemplate=(TH1*)templateFile->Get((opts.GetOpt<string>(channel+".templateFit.file", 1)+templateTag).c_str());
                if(wfTemplate)
                {
                    templates_["PHYS"][channel] = (TH1F*) wfTemplate->Clone();
						  //std::cout << templates_["PHYS"][channel]->GetEntries() << std::endl;
                    templates_["PHYS"][channel]->SetDirectory(0);
                    WFs_[channel]->SetTemplate(templates_["PHYS"][channel]);
                    //std::cout << " ... SetTemplate for channel " << channel << std::endl;
                }
                else
                {
                    Log("template "+opts.GetOpt<string>(channel+".templateFit.file", 1)
                        +" not found in "+opts.GetOpt<string>(channel+".templateFit.file", 0), ERR);
                    return false;
                }                
            }
            else
            {
                Log("template file "+opts.GetOpt<string>(channel+".templateFit.file", 0)+" not found", ERR);
                return false;
            }
            templateFile->Close();
        }
        // added from Simone S.
        if(opts.OptExist(channel+".templateFit.laserFile"))
        {
            TFile* templateFile = TFile::Open(opts.GetOpt<string>(channel+".templateFit.laserFile", 0).c_str(), ".READ");
            if(templateFile)
            {
                TH1* wfTemplate=(TH1*)templateFile->Get((opts.GetOpt<string>(channel+".templateFit.laserFile", 1)+templateTag).c_str());
                if(wfTemplate)
                {
                    templates_["LASER"][channel] = (TH1F*) wfTemplate->Clone();
                    templates_["LASER"][channel]->SetDirectory(0);
                }
                else
                {
                    Log("template "+opts.GetOpt<string>(channel+".templateFit.laserFile", 1)
                        +" not found in "+opts.GetOpt<string>(channel+".templateFit.laserFile", 0), ERR);
                    return false;
                }                
            }
            else
            {
                Log("template file "+opts.GetOpt<string>(channel+".templateFit.laserFile", 0)+" not found", ERR);
                return false;
            }
            templateFile->Close();
        }
        //
        if(opts.OptExist(channel+".templateFit.spikeFile"))
        {
            TFile* templateFile = TFile::Open(opts.GetOpt<string>(channel+".templateFit.spikeFile", 0).c_str(), ".READ");
            TH1* wfTemplate=(TH1*)templateFile->Get((opts.GetOpt<string>(channel+".templateFit.spikeFile", 1)+templateTag).c_str());
            spikeTemplates_[channel] = (TH1F*) wfTemplate->Clone();
            spikeTemplates_[channel] -> SetDirectory(0);
            templateFile->Close();
        }
        
	WFClassClock* isclock = dynamic_cast<WFClassClock*>(WFs_[channel]);
	if (isclock != NULL && opts.OptExist(channel+".clkPeriod")) 
            WFs_[channel]->SetPeriod(opts.GetOpt<float>(channel+".clkPeriod"));
        
        //---keep track of all the possible time reco method requested
        for(auto type_name : timeRecoTypes_)
        {
            if(opts.OptExist(channel+"."+type_name))
                timeOpts_[channel+"."+type_name] = opts.GetOpt<vector<float> >(channel+"."+type_name);
        }
    }
    
    //---outputs---
    string digiTreeName = opts.OptExist(instanceName_+".digiTreeName") ?
        opts.GetOpt<string>(instanceName_+".digiTreeName") : "digi";
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    RegisterSharedData(new TTree(digiTreeName.c_str(), "digi_tree"), "digi_tree", storeTree);
    digiTree_ = DigiTree(index, (TTree*)data_.back().obj);
    digiTree_.Init(channelsNames_, timeRecoTypes_);
    //---store channel time and amplitude calibrations---
    int outCh=0;
    for(auto& channel : channelsNames_)
    {
        digiTree_.ampl_calib[outCh] = opts.OptExist(channel+".amplitudeCalib") ? 
            opts.GetOpt<float>(channel+".amplitudeCalib") : 1.;
        digiTree_.time_calib[outCh] = opts.OptExist(channel+".timeCalib") ? 
            opts.GetOpt<float>(channel+".timeCalib") : 1.;
        ++outCh;
    }

    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
    {
        string wfTreeName = opts.OptExist(instanceName_+".wfTreeName") ?
            opts.GetOpt<string>(instanceName_+".wfTreeName") : "wf";
        RegisterSharedData(new TTree(wfTreeName.c_str(), "wf_tree"), "wf_tree", true);
        outWFTree_ = WFTree(nSamples, index, (TTree*)data_.back().obj);
        outWFTree_.Init();
    }

    eventAnalyzer_ = RecoEventAnalyzer(plugins, index);
    
    return true;
}

bool WFAnalyzer::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---reset output event 
    int outCh=0;
    bool fillWFtree=false;
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
        fillWFtree = *digiTree_.index % opts.GetOpt<int>(instanceName_+".WFtreePrescale") == 0;

    //---Check if trigger bit has changed from previous event        
    auto ctrg = "PHYS";
    
    if (trgInstance_!="")
      ctrg=((TObjString*)(plugins[trgInstance_]->GetSharedData(trgInstance_+"_trg_bit", "", false))[0].obj)->GetString().Data();

    if(ctrg != trg_ && 
       templates_.find(ctrg) != templates_.end())
      {
	for(auto& channel : channelsNames_)
	  {
	    if(templates_[ctrg].find(channel) != templates_[ctrg].end())
	      WFs_[channel]->SetTemplate(templates_[ctrg][channel]);       
	  }                                
    trg_ = ctrg;
      }

    //---compute reco variables
    for(auto& channel : channelsNames_)
    {
        //---skip dead channels or channels with too few samples
        if(WFs_.find(channel) == WFs_.end() || WFs_[channel]->GetNSample() < opts.GetOpt<int>(channel+".signalWin", 0) + 10)
        {
            digiTree_.FillVoidChannel(outCh);
            ++outCh;
            continue;
        }

        //---subtract a specified channel if requested
        if(opts.OptExist(channel+".subtractChannel") && WFs_.find(opts.GetOpt<string>(channel+".subtractChannel")) != WFs_.end())
            *WFs_[channel] -= *WFs_[opts.GetOpt<string>(channel+".subtractChannel")];        
	
	if(opts.OptExist(channel+".baselineWin"))
        { 
	    WFs_[channel]->SetBaselineWindow(opts.GetOpt<int>(channel+".baselineWin", 0), 
					     opts.GetOpt<int>(channel+".baselineWin", 1));
        }
	if(opts.OptExist(channel+".baselineInt"))
        {
	    WFs_[channel]->SetBaselineIntegralWindow(opts.GetOpt<int>(channel+".baselineInt", 0),
						     opts.GetOpt<int>(channel+".baselineInt", 1));
        }
	
	WFs_[channel]->SetSignalWindow(opts.GetOpt<int>(channel+".signalWin", 0), 
				       opts.GetOpt<int>(channel+".signalWin", 1));
	
        WFs_[channel]->SetSignalIntegralWindow(opts.GetOpt<int>(channel+".signalInt", 0),
                                               opts.GetOpt<int>(channel+".signalInt", 1));
	WFBaseline baselineInfo; 
	if (opts.OptExist(channel+".baseline"))
	  baselineInfo = WFs_[channel]->SubtractBaseline(opts.GetOpt<float>(channel+".baseline"));
	else if (opts.OptExist(channel+".baselineFit"))
	  baselineInfo = WFs_[channel]->SubtractBaselineFit();
	else
	  baselineInfo = WFs_[channel]->SubtractBaseline();
	
	//FIXME MAREMMA MAIALA
        WFFitResults interpolAmpMax;
        if(opts.OptExist(channel+".signalWin", 4))
        {
            string max_function = opts.OptExist(channel+".signalWin", 4) ? opts.GetOpt<string>(channel+".signalWin", 4) : "pol2";
            int nParams = opts.OptExist(channel+".signalWin", 5) ? opts.GetOpt<int>(channel+".signalWin", 5) : 0;
            std::vector<float> max_params;
            for(int ipar = 0; ipar < nParams; ++ipar)
                max_params.push_back( opts.GetOpt<float>(channel+".signalWin", 6+ipar) );
            interpolAmpMax = WFs_[channel]->GetInterpolatedAmpMax(-1,-1,
                                                                  opts.GetOpt<int>(channel+".signalWin", 2),
                                                                  opts.GetOpt<int>(channel+".signalWin", 3), max_function, max_params);
        }
        else
        {
            string max_function = opts.OptExist(channel+".signalWin", 3) ? opts.GetOpt<string>(channel+".signalWin", 3) : "pol2";
            interpolAmpMax = WFs_[channel]->GetInterpolatedAmpMax(-1,-1,
                                                                  opts.GetOpt<int>(channel+".signalWin", 2)/2,
                                                                  opts.GetOpt<int>(channel+".signalWin", 2)/2,
                                                                  max_function);
				//if(channel == "C2") std::cout << "  MAX amp " << interpolAmpMax.ampl << " at " << interpolAmpMax.time << std::endl;
        }
        digiTree_.pedestal[outCh] = baselineInfo.baseline;
	WFClassLiTeDTU* islitedtu = dynamic_cast<WFClassLiTeDTU*>(WFs_[channel]);
	if (islitedtu != NULL) digiTree_.gain[outCh] = WFs_[channel]->GetGain();
        digiTree_.b_charge[outCh] = WFs_[channel]->GetIntegral(opts.GetOpt<int>(channel+".baselineInt", 0), 
                                                               opts.GetOpt<int>(channel+".baselineInt", 1));        
        digiTree_.b_slope[outCh] = baselineInfo.slope;
        digiTree_.b_rms[outCh] = baselineInfo.rms;
        digiTree_.maximum[outCh] = WFs_[channel]->GetAmpMax();
        digiTree_.time_maximum[outCh] = WFs_[channel]->GetTimeCF(1).time;
        digiTree_.amp_max[outCh] = interpolAmpMax.ampl;
        digiTree_.time_max[outCh] = interpolAmpMax.time;
        digiTree_.chi2_max[outCh] = interpolAmpMax.chi2;
        digiTree_.charge_tot[outCh] = WFs_[channel]->GetModIntegral(opts.GetOpt<int>(channel+".baselineInt", 1), 
                                                                    WFs_[channel]->GetNSample());
        digiTree_.charge_sig[outCh] = WFs_[channel]->GetSignalIntegral(opts.GetOpt<int>(channel+".signalInt", 0), 
                                                                       opts.GetOpt<int>(channel+".signalInt", 1));
        //---compute time with all the requested time reconstruction method
        for(unsigned int iT=0; iT<timeRecoTypes_.size(); ++iT)
        {
            //---compute time with selected method or store default value (-99)
            if(timeOpts_.find(channel+"."+timeRecoTypes_[iT]) != timeOpts_.end() && WFs_[channel]->GetAmpMax() > 100.)
            {
                WFFitResults timeInfo = WFs_[channel]->GetTime(timeRecoTypes_[iT], timeOpts_[channel+"."+timeRecoTypes_[iT]]);
                digiTree_.time[outCh+iT*channelsNames_.size()] = timeInfo.time;
                digiTree_.time_error[outCh+iT*channelsNames_.size()] = timeInfo.error;                
                digiTree_.time_chi2[outCh+iT*channelsNames_.size()] = timeInfo.chi2;
                digiTree_.time_slope[outCh+iT*channelsNames_.size()] = timeInfo.slope;
            }
            else
            {
                digiTree_.time[outCh+iT*channelsNames_.size()] = -99;
                digiTree_.time_error[outCh+iT*channelsNames_.size()] = -99;
                digiTree_.time_chi2[outCh+iT*channelsNames_.size()] = -99;
                digiTree_.time_slope[outCh+iT*channelsNames_.size()] = -99;
            }
        }
        digiTree_.period[outCh] = WFs_[channel]->GetPeriod();
        //---template fit (only specified channels)
        if(opts.OptExist(channel+".templateFit.file"))
        {
				//std::cout << " Input to TemplateFit() : " << (opts.OptExist(channel+".templateFit.amplitudeThreshold") ? opts.GetOpt<float>(channel+".templateFit.amplitudeThreshold") : 0 )<< "\t" << opts.GetOpt<float>(channel+".templateFit.fitWin", 0) << "\t" << opts.GetOpt<int>(channel+".templateFit.fitWin", 1) << "\t" << opts.GetOpt<int>(channel+".templateFit.fitWin", 2) << std::endl;
            auto fitResults = WFs_[channel]->TemplateFit(
                opts.OptExist(channel+".templateFit.amplitudeThreshold") ? opts.GetOpt<float>(channel+".templateFit.amplitudeThreshold") : 0,
                opts.GetOpt<float>(channel+".templateFit.fitWin", 0),
                opts.GetOpt<int>(channel+".templateFit.fitWin", 1),
                opts.GetOpt<int>(channel+".templateFit.fitWin", 2));
            digiTree_.fit_ampl[outCh] = fitResults.ampl;
            digiTree_.fit_time[outCh] = fitResults.time;
            digiTree_.fit_terr[outCh] = fitResults.error;            
            digiTree_.fit_chi2[outCh] = fitResults.chi2;
            digiTree_.fit_period[outCh] = WFs_[channel]->GetTemplateFitPeriod();

            //if(channel == "C2" || channel == "MCP1" || channel == "MCP2") //std::cout << "--> Channel " << channel <<  " EV " << event.evtNumber << " Fit_ampl  ="<< fitResults.ampl << "\t Chi2 " << fitResults.chi2 << std::endl;

            WFFitResultsScintPlusSpike fitResultsScintPlusSpike{-1, -1000, -1, -1000, -1};
            if(opts.OptExist(channel+".templateFit.spikeFile"))
            {
                WFs_[channel]->SetTemplateScint(templates_["PHYS"][channel]);
                WFs_[channel]->SetTemplateSpike(spikeTemplates_[channel]);
                fitResultsScintPlusSpike = WFs_[channel]->TemplateFitScintPlusSpike(
                    opts.OptExist(channel+".templateFit.amplitudeThreshold") ? opts.GetOpt<float>(channel+".templateFit.amplitudeThreshold") : 0,
                    opts.GetOpt<float>(channel+".templateFit.fitWin", 0),
                    opts.GetOpt<int>(channel+".templateFit.fitWin", 1),
                    opts.GetOpt<int>(channel+".templateFit.fitWin", 2));
                digiTree_.fit_ampl_scint[outCh] = fitResultsScintPlusSpike.ampl_scint;
                digiTree_.fit_time_scint[outCh] = fitResultsScintPlusSpike.time_scint;
                digiTree_.fit_ampl_spike[outCh] = fitResultsScintPlusSpike.ampl_spike;
                digiTree_.fit_time_spike[outCh] = fitResultsScintPlusSpike.time_spike;
                digiTree_.fit_chi2_scint_plus_spike[outCh] = fitResultsScintPlusSpike.chi2;
                digiTree_.fit_converged_scint_plus_spike[outCh] = fitResultsScintPlusSpike.converged;
            }
        }            
	
        //---WFs---
        if(fillWFtree)
        {
            auto analizedWF = WFs_[channel]->GetSamples();
            auto sampleTimes = WFs_[channel]->GetTimes();
            for(unsigned int jSample=0; jSample<analizedWF->size(); ++jSample)
            {
                outWFTree_.WF_ch.push_back(outCh);
                outWFTree_.WF_time.push_back(sampleTimes->at(jSample));
                outWFTree_.WF_val.push_back(analizedWF->at(jSample));
            }
        }
        //---increase output tree channel counter
        ++outCh;
    }

    //---fill the output trees 
    //---reco var
    digiTree_.Fill();
    //---WFs
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
        outWFTree_.Fill();
    
    return true;
}
