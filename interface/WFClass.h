#ifndef __WFCLASS_H__
#define __WFCLASS_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TError.h"

#include "interface/H4Tree.h"

using namespace std;

#define MAX_INTERPOLATOR_POINTS 10000

struct WFBaseline
{
    float baseline;
    float rms;
    float slope;
    float k;
    float chi2;
};

struct WFFitResults
{
    double ampl;
    double time;
    double error;
    double chi2;
    double slope;
};

struct WFFitResultsScintPlusSpike
{
    double ampl_scint;
    double time_scint;
    double ampl_spike;
    double time_spike;
    double chi2;
    bool converged;
};

//---Single channel calibration
//   The index correspond to startCellIndex
struct calibrationParameters
{
    double deltaV;
    double slopeV;
    double quadraticV;
    double deltaT;
};
typedef std::vector<calibrationParameters> DigiChannelCalibration;

//---Digitizer calibration container:
//   just a ROOT wrapper for an unordered map.
class DigitizerCalibration :
    public TNamed,
    public std::unordered_map<const bgc_key_t, DigiChannelCalibration, key_hash, key_equal>
{
public:

    DigitizerCalibration() : TNamed(), std::unordered_map<const bgc_key_t, DigiChannelCalibration, key_hash, key_equal>()
        {
            clear();
        }

    virtual ~DigitizerCalibration()
        {
        }

    ClassDef(DigitizerCalibration, 1)
};

class WFClass : public TObject
{
public:
    //---ctors---
    WFClass():
        uncalibSamples_(), samples_(uncalibSamples_)
        {};
    WFClass(int polarity, float tUnit, DigiChannelCalibration* calibration=NULL);
    //---dtor---
    ~WFClass() {};

    //crop method
    void CropWF(WFClass& WF, int firstSample=-1, int lastSample=-1);

    //---getters---
    inline const vector<double>*   GetSamples() {return &samples_;};
    inline const vector<double>*   GetTimes() {return &times_;};
    inline const vector<double>*   GetOriginalSamples() {return &uncalibSamples_;};    
    inline int                     GetStartIndexCell() {return startIndexCell_;}
    inline int                     GetBWinMin() {return bWinMin_;}
    inline int                     GetBWinMax() {return bWinMax_;}
    inline int                     GetBIntWinMin() {return bIntWinMin_;}
    inline int                     GetBIntWinMax() {return bIntWinMax_;}
    inline float                   GetBaseline() {return baseline_;}
    inline float                   GetBaselineRMS() {return bRMS_;}
    inline int                     GetNSample() {return samples_.size();};
    inline int                     GetMaxSample() {return maxSample_;};
    inline float                   GetTUnit() {return tUnit_;};
    inline int                     GetSWinMin() {return sWinMin_;}
    inline int                     GetSWinMax() {return sWinMax_;}
    inline int                     GetSIntWinMin() {return sIntWinMin_;}
    inline int                     GetSIntWinMax() {return sIntWinMax_;}
    inline float                   GetFitAmpMax() {return fitAmpMax_;};
    inline float                   GetFitTimeMax() {return fitTimeMax_*tUnit_;};
    inline float                   GetLEThr() {return leThr_;};
    inline float                   GetLETime() {return leTime_;};
    inline float                   GetLESlope() {return leSlope_;};
    inline float                   GetTEThr() {return teThr_;};
    inline float                   GetTETime() {return teTime_;};
    inline float                   GetTESlope() {return teSlope_;};
    inline float                   GetCFFrac() {return cfFrac_;};
    inline float                   GetCFTime() {return cfTime_;};
    inline float                   GetCFSlope() {return cfSlope_;};
    inline TF1*                    GetAmpFunc() { return f_max_; };
    inline TF1*                    GetFitFunc() { return f_fit_; };
    inline DigiChannelCalibration* GetCalibration() { return calibration_; };
    inline int                     GetGain() { return *std::max_element(gain_.begin(), gain_.end()); };
    float                          GetAmpMax(int min=-1, int max=-1);
    WFFitResults                   GetInterpolatedSample(int sample, int samplesLeft=-1, int samplesRight=-1);
    WFFitResults                   GetInterpolatedAmpMax(int min=-1, int max=-1, int nmFitSamples=7, int npFitSamples=7, string function="pol2", vector<float> params=vector<float>{});
    virtual WFFitResults           GetTime(string method, vector<float>& params); 
    virtual WFFitResults           GetTimeCF(float frac, int nFitSamples=5, int min=-1, int max=-1);
    virtual WFFitResults           GetTimeLE(float thr, int nmFitSamples=1, int npFitSamples=3, int min=-1, int max=-1);
    virtual WFFitResults           GetTimeTE(float thr, int nmFitSamples=1, int npFitSamples=3, int min=-1, int max=-1);
    float                          GetIntegral(int min=-1, int max=-1);
    float                          GetModIntegral(int min=-1, int max=-1);
    virtual float                  GetSignalIntegral(int riseWin, int fallWin);
    virtual float                  GetPeriod() { return -1; };
    virtual float                  GetTemplateFitPeriod() { return -1; };
    
    //---setters---
    inline void                    SetTrigRef(float trigRef){trigRef_ = trigRef;};
    inline void                    SetCalibration(DigiChannelCalibration* calib){calibration_=calib;};
    inline void                    SetStartIndexCell(int cell){startIndexCell_=cell;};
    virtual void                   SetTemplateScint(TH1* templateWF=NULL) { return; };
    virtual void                   SetTemplateSpike(TH1* templateWF=NULL) { return; };
    void                           SetSignalWindow(int min, int max);
    void                           SetSignalIntegralWindow(int min, int max);
    void                           SetBaselineWindow(int min, int max);
    void                           SetBaselineIntegralWindow(int min, int max);
    virtual void                   SetTemplate(TH1* templateWF=NULL);
    virtual inline void            SetPeriod(float period) { return; };

    //---utils---
    virtual void                           Reset();
    bool                                   ApplyCalibration();
    virtual void                           AddSample(float sample);
    virtual void                           AddSample(float sample, float gain) {AddSample(sample);};
    WFBaseline                             SubtractBaseline(int min=-1, int max=-1);
    WFBaseline                             SubtractBaseline(float baseline);
    virtual WFBaseline                     SubtractBaselineFit(int min=-1, int max=-1) {return WFBaseline{0, 0, 0, 0, 999};};
				           
    virtual WFFitResults                   TemplateFit(float amp_threshold=0., float offset=0., int lW=0, int hW=0);
				           
    double                                 AnalyticFit(TF1* f, int lW, int hW);
    void                                   EmulatedWF(WFClass& wf, float rms, float amplitude, float time);
    void                                   FFT(WFClass& wf, float tau, int cut);
    void                                   Print();
    virtual WFFitResultsScintPlusSpike     TemplateFitScintPlusSpike(float amp_threshold=0., float offset=0., int lW=0, int hW=0) { return {0}; };
    //---operators---
    WFClass&                       operator=(const WFClass& origin);
    WFClass                        operator-(const WFClass& sub);
    WFClass                        operator+(const WFClass& add);
    WFClass&                       operator-=(const WFClass& sub);
    WFClass&                       operator+=(const WFClass& add);
protected:
    //---utils---
    float                          BaselineRMS();
    float                          LinearInterpolation(float& A, float& B, const int& min, const int& max, const int& skipSample=-1);
    virtual  double                TemplateChi2(const double* par=NULL);
    double                         AnalyticChi2(const double* par=NULL);
    
protected:
    vector<double>  uncalibSamples_;
    vector<double>  calibSamples_;
    vector<double>  samples_;
    vector<double>  gain_;
    vector<double>  times_;

    DigiChannelCalibration* calibration_;

    int            startIndexCell_;
    float          tUnit_;
    int            polarity_;
    float          trigRef_;
    int            sWinMin_;
    int            sWinMax_;
    int            sIntWinMin_;
    int            sIntWinMax_;
    int            bWinMin_;
    int            bWinMax_;
    int            bIntWinMin_;
    int            bIntWinMax_;
    int            maxSample_;
    float          fitAmpMax_;
    float          fitTimeMax_;
    float          fitChi2Max_;
    float          baseline_;
    float          bRMS_;
    int            cfSample_;
    float          cfFrac_;
    float          cfTime_;
    float          cfSlope_;
    int            leSample_;
    float          leThr_;
    float          leTime_;
    float          leSlope_;
    int            teSample_;
    float          teThr_;
    float          teTime_;
    float          teSlope_;
    float          chi2cf_;
    float          chi2le_;
    float          chi2te_;
    int            fWinMin_;
    int            fWinMax_;
    float          tmplFitTime_;
    float          tmplFitTimeErr_;
    float          tmplFitAmp_;
    float          tmplFitAmpShift_;
    TF1*           f_max_;
    TF1*           f_fit_;
    float          interpolatorMin_;
    float          interpolatorMax_;         
    ROOT::Math::Interpolator* interpolator_;
};

#endif
