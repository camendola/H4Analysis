#ifndef __WFCLASS_LITEDTU_H__
#define __WFCLASS_LITEDTU_H__

#include "WFClass.h"
#include "TGraph.h"

using namespace std;

class WFClassLiTeDTU : public WFClass
{
public:
    //---ctors---
    WFClassLiTeDTU() {};
    WFClassLiTeDTU(int polarity, float tUnit, DigiChannelCalibration* calibration=NULL);
    //---setters---
    void                           SetTemplateScint(TH1* templateWF=NULL) override; 
    void                           SetTemplateSpike(TH1* templateWF=NULL) override;
    
    //---getters---
    inline int                     GetGain() { return *std::max_element(gain_.begin(), gain_.end()); };

    //---operators---
    WFClassLiTeDTU&                operator=(const WFClassLiTeDTU& origin);

    //---utils---
    void                           Reset() override;    
    double                         TemplatesChi2(const double* par=NULL);
    void                           AddSample(float sample, float gain) override; 
    WFFitResults                   TemplateFit(float ampl_threshold=0, float offset=0., int lW=0, int hW=0) override;
    WFFitResultsScintPlusSpike     TemplateFitScintPlusSpike(float amp_threshold=0., float offset=0., int lW=0, int hW=0) override; 
    
 protected: 
    float          tmplFitTimeScint_;
    float          tmplFitAmpScint_;
    float          tmplFitTimeSpike_;
    float          tmplFitAmpSpike_;
    bool           tmplFitConverged_;
    float          tmplTimeMaxScint_;
    float          tmplTimeMaxSpike_;
    
    ROOT::Math::Interpolator* interpolatorScint_;
    ROOT::Math::Interpolator* interpolatorSpike_;
  
};

#endif
