#ifndef __H4_RECO__
#define __H4_RECO__

#include <netdb.h>
#include <unistd.h>
#include <csignal>
#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <experimental/filesystem>

#include "TFile.h"
#include "TChain.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/PluginLoader.h"
#include "interface/PluginBase.h"
#include "interface/DataLoader.h"
#include "interface/RecoTree.h"

namespace fs = std::experimental::filesystem;

//**********UTILS*************************************************************************
std::string getMachineDomain()
{
    char hn[254];
    char *dn;
    struct hostent *hp;

    gethostname(hn, 254);
    hp = gethostbyname(hn);
    dn = strchr(hp->h_name, '.');
    if ( dn != NULL ) 
    {
        return std::string(++dn);
    }
    else 
    {
        return "";
    }
}

//----------Simple function to track memory and CPU usage---------------------------------
void TrackProcess(float* cpu, float* mem, float* vsz, float* rss)
{
    string dummy1, dummy2, dummy3, dummy4, dummy5;
    string time;
    //---get cpu/mem info
    int pid = getpid();
    string ps_command = "ps up "+to_string(pid)+" >.H4Reco_proc.tmp";
    system(ps_command.c_str());
    ifstream proc_tmp(".H4Reco_proc.tmp", ios::in);
    getline(proc_tmp, dummy1);
    proc_tmp >> dummy1 >> dummy2 >> cpu[0] >> mem[0] >> vsz[0] >> rss[0]
             >> dummy3 >> dummy4 >> dummy5 >> time;
    vsz[0] = vsz[0]/1000;
    rss[0] = rss[0]/1000;
    proc_tmp.close();
    if(cpu[0]>cpu[1])
        cpu[1] = cpu[0];
    if(mem[0]>mem[1])
        mem[1] = mem[0];
    if(vsz[0]>vsz[1])
        vsz[1] = vsz[0];
    if(rss[0]>rss[1])
        rss[1] = rss[0];
    //---print statistics
    cout << "\033[0m""-----Machine stats---current/max-----" << endl
         << "CPU(%): " << cpu[0] << "/" << cpu[1] << endl
         << "MEM(%): " << mem[0] << "/" << mem[1] << endl
         << "VSZ(M): " << vsz[0] << "/" << vsz[1] << endl
         << "RSS(M): " << rss[0] << "/" << rss[1] << endl
         << "time lasted: " << time << endl;
}

//----------Exception handler-------------------------------------------------------------
//---This fuction is meant to provide debug information by printing
//---the last called Plugin::Method while re-throwing the original exception
void HandleException(std::exception_ptr eptr, PluginBase* plugin)
{
    try
    {
        if (eptr)
        {
            std::rethrow_exception(eptr);
        }
    }
    catch(const std::exception& e)
    {
        auto plugin_type = std::regex_replace(typeid(plugin).name(), std::regex(".*[0-9]+"), "");
        std::cout << "\033[1;31m" << ">>>>> H4Reco ERROR! <<<<<" << "\033[0m" << std::endl
                  << "Error in: " << "\033[1;33m" << plugin_type << "::" << plugin->GetCurrentMethod() << "\033[0m" << std::endl
                  << "Plugin type: " << "\033[1;33m" << plugin->GetPluginType() << "\033[0m" << std::endl
                  << "Instance name: " << "\033[1;33m" << plugin->GetInstanceName() << "\033[0m" << std::endl
                  << "Caught exception: " << e.what() << std::endl;
        exit(-1);
    }
}

//**********MAIN**************************************************************************
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << argv[0] << " cfg file " << "[run] " << "[first spill] " << "[number of spills] " << endl;
        return -1;
    }

    //---exception handler---
    std::exception_ptr eptr;
    
    //---memory consumption tracking---
    float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};

    //---load options---    
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);

    //-----input setup-----    
    int spill=-1;
    if(argc > 2)
    {
        vector<string> run(1, argv[2]);
        opts.SetOpt("h4reco.run", run);
    }
    if(argc > 3) {
        spill = atoi(argv[3]);
        vector<string> files(1, "1");
        opts.SetOpt("h4reco.maxFiles", files);
    }
    if(argc > 4)
    {
        vector<string> files(1, argv[4]);
        opts.SetOpt("h4reco.maxFiles", files);
    }
    auto out_file_name = opts.GetOpt<string>("h4reco.outNameSuffix");
    bool storeTree = opts.OptExist("h4reco.storeTree") ?
        opts.GetOpt<bool>("h4reco.storeTree") : true;

    string run = opts.GetOpt<string>("h4reco.run");
    int totLoops= opts.OptExist("h4reco.totLoops") ? opts.GetOpt<int>("h4reco.totLoops") : 1;


    //-----Load raw data-----
    vector<string> spillOpt(1, to_string(spill));
    opts.SetOpt("h4reco.firstSpill", spillOpt);
    DataLoader dataLoader(opts);

    //-----output setup-----
    uint64 index=0;
    if(spill == -1)
        out_file_name += run+".root";
    else if(opts.GetOpt<int>("h4reco.maxFiles") > 1)
        out_file_name += run+"_"+spillOpt.back()+".root";
    else
        out_file_name += "/"+run+"/"+to_string(spill)+".root";
    fs::create_directories(fs::absolute(fs::path(out_file_name.substr(0, out_file_name.find_last_of("/")+1))));
    auto* outROOT = new TFile(out_file_name.c_str(), "RECREATE");
    outROOT->cd();
    TDirectory* outDIR=outROOT;

    RecoTree mainTree(&index);

    //---Get plugin sequence---
    PluginLoader<PluginBase>* loader;
    vector<PluginLoader<PluginBase>* > pluginLoaders;    
    map<string, PluginBase*> pluginMap;
    vector<PluginBase*> pluginSequence;
    vector<string> pluginList = opts.GetOpt<vector<string> >("h4reco.pluginList");    
    //---plugin creation
    pluginLoaders.reserve(pluginList.size());
    for(auto& plugin : pluginList)
    {
        cout << ">>> Loading plugin <" << plugin << ">" << endl;
        //---create loader 
        loader = new PluginLoader<PluginBase>(opts.GetOpt<string>(plugin+".pluginType"));
        pluginLoaders.push_back(loader);
        pluginLoaders.back()->Create();
        //---get instance and put it in the plugin sequence   
        PluginBase* newPlugin = pluginLoaders.back()->CreateInstance(plugin);
        if(newPlugin)
        {
            pluginSequence.push_back(newPlugin);
            pluginMap[plugin] = pluginSequence.back();
        }
        else
        {
            cout << ">>> ERROR: plugin type " << opts.GetOpt<string>(plugin+".pluginType") << " is not defined." << endl;
            return 0;
        }
    }

    //---begin
    for(auto& plugin : pluginSequence)
    {
        //---call Begin() methods and check the return status
        try
        {
            //---fake call to base class for debug porpouses
            plugin->PluginBase::Begin(pluginMap, opts, &index);
            //---real call
            bool r_status = plugin->Begin(pluginMap, opts, &index);
            if(!r_status)
            {
                cout << ">>> ERROR: plugin returned bad flag from Begin() call: " << plugin->GetInstanceName() << endl;
                exit(-1);
            }
        }
        catch(...)
        {
            eptr = std::current_exception();
        }
        HandleException(eptr, plugin);
    }
            

    int iLoop=0;
    int maxEvents = opts.OptExist("h4reco.maxEvents") ? opts.GetOpt<int>("h4reco.maxEvents") : -1;
    bool isSim = opts.OptExist("h4reco.generateEvents") ? opts.GetOpt<bool>("h4reco.generateEvents") : false;

    if(isSim)
        cout << "\033[1;36m" << ">>> Processing H4DAQ simulation <<<" << "\033[0m" << endl;
    else
        cout << "\033[1;36m" << ">>> Processing H4DAQ run #" << run << " <<<" << "\033[0m" << endl;

    bool endLoop=true;
    while(iLoop<totLoops && endLoop)
    {
	int nEvents = 0;

	if (totLoops>1)
        {
	    cout << "\033[1;36m" << ">>> Starting iteration #" << iLoop << " <<<" << "\033[0m" << endl;
	    outDIR=outROOT->mkdir(("Iter_"+to_string(iLoop)).c_str());
        }

	//---begin loop
	mainTree.tree_->Reset();
	for(auto& plugin : pluginSequence)
        {
	    //---call BeginLoop() methods and check the return status
            try
            {
                //---fake call to base class for debug porpouses
                plugin->PluginBase::BeginLoop(iLoop, pluginMap, opts);
                //---real call
                bool r_status = plugin->BeginLoop(iLoop, pluginMap, opts);
                if(!r_status)
                {
                    cout << ">>> ERROR: plugin returned bad flag from BeginLoop() call: " << plugin->GetInstanceName() << endl;
                    exit(-1);
                }
            }
            catch(...)
            {
                eptr = std::current_exception();
            }
            HandleException(eptr, plugin);
	    
	    //---Get plugin shared data
	    for(auto& shared : plugin->GetSharedData("", "TTree", true))
            {
	    	TTree* tree = (TTree*)shared.obj;
		tree->Reset();
	    	tree->SetDirectory(outDIR);
            }
        }


	//---events loop
	while((dataLoader.NextEvent() && (nEvents < maxEvents || maxEvents == -1)) || (isSim && (nEvents < maxEvents)))
        {
       bool time_aligned = true;
       vector<int> VFE_ch = {1,3,5,7,9};
       for(auto& ch : VFE_ch){   
          //std::cout << "VFE channel " << ch << std::endl;
          if((long long int)(dataLoader.GetTree().evtTime[0] - dataLoader.GetTree().evtTime[ch]) > 1500) time_aligned = false;
          if((long long int)(dataLoader.GetTree().evtTime[0] - dataLoader.GetTree().evtTime[ch]) <  200) time_aligned = false;
       }
       if((long long int)(dataLoader.GetTree().evtTime[11] - dataLoader.GetTree().evtTime[0]) > 200.) time_aligned = false;
       //if(!time_aligned) continue;
	    if(dataLoader.FirstEventInSpill())
            {
                cout << "\033[1;36m" << ">>> Processed spills: " << dataLoader.GetNFilesProcessed() << "/" << dataLoader.GetNFiles() << endl;
                cout << ">>> Processed events: " << nEvents << "\033[0m" << endl;
		TrackProcess(cpu, mem, vsz, rss);
            }
        
	    //---set index value run*1e10+spill*1e4+event
	    index = dataLoader.GetTree().runNumber*1e9 + dataLoader.GetTree().spillNumber*1e5 + dataLoader.GetTree().evtNumber;
	    
	    //---call ProcessEvent for each plugin and check the return status
	    bool status=true;
	    for(auto& plugin : pluginSequence)
            {
                try
                {
                    //---fake call to base class for debug porpouses
                    plugin->PluginBase::ProcessEvent(dataLoader.GetTree(), pluginMap, opts);
                    //---real call + check for filters
                    if(status)
                        status &= plugin->ProcessEvent(dataLoader.GetTree(), pluginMap, opts);
                }
                catch(...)
                {
                    eptr = std::current_exception();
                }
                HandleException(eptr, plugin);
            }

	    //    cbasile: check if the time-stamps synchronization hold
            //int nDigitizers = (int)(dataLoader.GetTree().nEvtTimes/2); 
            //std::cout << " --> nDigitizers " << nDigitizers << std::endl; 
            //for(int iD=0; iD<nDigitizers; ++iD){
            //    if(fabs((long long int)(dataLoader.GetTree().evtTime[2*iD+1] - dataLoader.GetTree().evtTime[0])) > 500.) status = false;   
            //}
        
	    //---Fill the main tree with info variables and increase event counter
            mainTree.time_stamps.clear();
            for(int iT=0; iT<dataLoader.GetTree().nEvtTimes; ++iT)
                mainTree.time_stamps.push_back(dataLoader.GetTree().evtTime[iT]);
            mainTree.evt_flag = status;
            mainTree.run = dataLoader.GetTree().runNumber;
            mainTree.spill = dataLoader.GetTree().spillNumber;
            mainTree.event = dataLoader.GetTree().evtNumber;
            mainTree.Fill();
            ++nEvents;
        }
    cout << ">>> TOTAL processed events: " << nEvents << "\033[0m" << endl;

	//---end
	for(auto& plugin : pluginSequence)
        {
	    //---call endjob for each plugin
            try
            {
                //---fake call to base class for debug porpouses
                plugin->PluginBase::EndLoop(iLoop, pluginMap, opts);
                //---real call
                endLoop &= plugin->EndLoop(iLoop, pluginMap, opts);
            }
            catch(...)
            {
                eptr = std::current_exception();
            }
            HandleException(eptr, plugin);

	    //---get permanent data from each plugin and store them in the out file
	    for(auto& shared : plugin->GetSharedData())
            {
	    	if(shared.obj->IsA()->GetName() == string("TTree") && storeTree)
                {
	    	    TTree* currentTree = (TTree*)shared.obj;
	    	    outDIR->cd();
                    currentTree->BuildIndex("index");
	    	    currentTree->Write(currentTree->GetName(), TObject::kOverwrite);
	    	    mainTree.AddFriend(currentTree->GetName());
                }
	    	else
                {
	    	    outDIR->cd();
	    	    shared.obj->Write(shared.tag.c_str(), TObject::kOverwrite);
                }
            }
        }

	outDIR->cd();
	if(storeTree)
            mainTree.Write();
	++iLoop;
    }

    //---end
    for(auto& plugin : pluginSequence)
    {
        //---call endjob for each plugin        
        try
        {
            //---fake call to base class for debug porpouses
            plugin->PluginBase::End(pluginMap, opts);
            //---real call
            bool r_status = plugin->End(pluginMap, opts);

            if(!r_status)
            {
                cout << ">>> ERROR: plugin returned bad flag from End() call: " << plugin->GetInstanceName() << endl;
                exit(-1);
            }
        }
        catch(...)
        {
            eptr = std::current_exception();
        }
        HandleException(eptr, plugin);
    }
    
    //---close
    outROOT->cd();
    opts.Write("cfg");
    outROOT->Close();
    for(auto& loader : pluginLoaders)
        loader->Destroy();

    //---info
    TrackProcess(cpu, mem, vsz, rss);

    exit(0);
}    

#endif
