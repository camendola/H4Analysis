#include "TriggerTypeFilter.h"

//----------Utils-------------------------------------------------------------------------
bool TriggerTypeFilter::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    triggerBoard_ = opts.GetOpt<unsigned int>(instanceName_+".triggerBoard");

    filterEvents_ = opts.GetOpt<bool>(instanceName_+".filterEvents");
    if(filterEvents_)
    {
        filterName_ = opts.GetOpt<string>(instanceName_+".filterName");
        Log("Only events on type "+filterName_+" will be processed"); 
    }
    
    auto trg_masks = opts.GetOpt<vector<float> >(instanceName_+".triggerMasks");
    auto trg_names = opts.GetOpt<vector<string> >(instanceName_+".triggerNames");

    assert(trg_masks.size() == trg_names.size());
    std::transform(trg_masks.begin(), trg_masks.end(), trg_names.begin(), 
                   std::inserter(maskToName_, maskToName_.end()), [](float mask, string name)
                   {
                       return std::make_pair(int(mask), name);
                   });

    //---outputs---
    string trgTreeName = opts.OptExist(instanceName_+".treeName") ?
        opts.GetOpt<string>(instanceName_+".treeName") : "trg";
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    //---register tree
    RegisterSharedData(new TTree(trgTreeName.c_str(), "trg_tree"), "trg_tree", storeTree);
    trgTree_ = TrgTree(index, (TTree*)data_.back().obj);
    trgTree_.Init(maskToName_);
    
    //---register trigger bit
    trg_ = "";
    RegisterSharedData(&trg_, "trg_bit", false);

    return true;
}

bool TriggerTypeFilter::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{        
    //---check trigger type
    
    auto pos = std::find(event.triggerWordsBoard, event.triggerWordsBoard+event.nTriggerWords,
                         triggerBoard_);
    if(pos == event.triggerWordsBoard+event.nTriggerWords)
        Log("Trigger board not found", WARN);
    else
    {
        for(auto const& [mask, name] : maskToName_)
        {
            //---A bit set in trigger word means that the veto was set for that particular trigger
            //   hence the active trigger is the only one for which the bit is not set
            if((event.triggerWords[std::distance(event.triggerWordsBoard, pos)] & mask) != mask)
            {
                trg_ = name.c_str();
                trgTree_.trg = mask;
            }
        }
    }        
    
    trgTree_.Fill();

    return filterEvents_ ? trg_.GetString().Data() == filterName_ : true;
}
