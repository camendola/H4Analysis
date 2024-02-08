#ifndef __WFCLASS_CLOCK_H__
#define __WFCLASS_CLOCK_H__

#include "WFClass.h"
#include "TGraph.h"

using namespace std;

class WFClassClock : public WFClass
{
public:
    //---ctors---
    WFClassClock() {};
    WFClassClock(int polarity, float tUnit);

    WFBaseline             SubtractBaselineFit(int min=-1, int max=-1) override; 
    void                   SetTemplate(TH1* templateWF=NULL) override;

    //---getters---
    void         AddSample(float sample) override;
    WFFitResults GetTime(string method, vector<float>& params) override;
    WFFitResults GetTimeLE(float thr = 0, int nmFitSamples=2, int npFitSamples=2, int min=-1, int max=-1) override;
    WFFitResults GetTimeCLK(float wleft=-1.3, float wright=1.3, int min=130, int max=900);    
    WFFitResults TemplateFit(float ampl_threshold=0, float offset=0., int lW=0, int hW=0) override;
    float        GetPeriod() override { return clkPeriod_; };
    float        GetTemplateFitPeriod() override { return tmplFitPeriod_; };
    inline void   SetPeriod(float period) override { extPeriod_ = period; };

protected:
    float extPeriod_; 
    double       TemplateChi2(const double* par=NULL) override;

private:
    float clkPeriod_;
    float tmplFitPeriod_;
};

#endif
