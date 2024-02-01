#include "interface/WFClassLiTeDTU.h"

//**********Constructors******************************************************************
//
WFClassLiTeDTU::WFClassLiTeDTU(int polarity, float tUnit, DigiChannelCalibration* calibration):
    WFClass(polarity, tUnit), tmplFitTimeScint_(-1), tmplFitAmpScint_(-1), tmplFitTimeSpike_(-1), tmplFitAmpSpike_(-1),
    tmplFitConverged_(false), tmplTimeMaxScint_(0), tmplTimeMaxSpike_(0), interpolatorScint_(NULL), interpolatorSpike_(NULL)
{}

//**********Getters***********************************************************************

//---------Add waveform sample to the list of uncalibrated samples------------------------
//---sample is inserted at the end of uncalibSamples_
//---the times vector is filled with the uncalibrated sample time computed from the time unit
//---a gain for each sample can be added. This is stored in a separate vector as well multiplied to the sample value.
void WFClassLiTeDTU::AddSample(float sample, float gain) 
{
    uncalibSamples_.push_back(polarity_*sample*gain);
    gain_.push_back(gain);
    times_.push_back( (samples_.size()-1.)*tUnit_ );
    samples_ = uncalibSamples_;
};


//----------template fit to the WF--------------------------------------------------------
WFFitResults WFClassLiTeDTU::TemplateFit(float amp_threshold, float offset, int lW, int hW)
{
    double tmplFitChi2=0;
    if(tmplFitAmp_ == -1)
    {
        amp_threshold *= gain_[maxSample_];
        if(samples_[maxSample_]>amp_threshold)
        {
            //---set template fit window around maximum, [min, max)
            BaselineRMS();
            GetAmpMax();    
            fWinMin_ = maxSample_ + int(offset/tUnit_) - lW;
            fWinMax_ = maxSample_ + int(offset/tUnit_) + hW;
            //---setup minimization
            auto t0 = GetInterpolatedAmpMax().time;
            ROOT::Math::Functor chi2(this, &WFClassLiTeDTU::TemplateChi2, 2);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(100000);
            minimizer->SetMaxIterations(1000);
            minimizer->SetTolerance(1e-4);
            minimizer->SetPrintLevel(0);
            minimizer->SetFunction(chi2);
            minimizer->SetLimitedVariable(0, "amplitude", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
            minimizer->SetLimitedVariable(1, "deltaT", t0, 1e-3, times_[fWinMin_], times_[fWinMax_]);        
            //---fit
            minimizer->Minimize();
            tmplFitAmp_ = minimizer->X()[0];
            tmplFitTime_ = minimizer->X()[1];
            tmplFitTimeErr_ = minimizer->Errors()[1];

            delete minimizer;
        }
        else
            return GetInterpolatedAmpMax();
    }    

    return WFFitResults{tmplFitAmp_, tmplFitTime_, tmplFitTimeErr_, TemplateChi2()/(fWinMax_-fWinMin_+1-2), 0};
}

WFFitResultsScintPlusSpike WFClassLiTeDTU::TemplateFitScintPlusSpike(float amp_threshold, float offset, int lW, int hW)
{
    if(tmplFitAmpScint_ == -1 && tmplFitAmpSpike_ == -1)
    {
        amp_threshold *= gain_[maxSample_];
        if (samples_[maxSample_] > amp_threshold) {
            //---set template fit window around maximum, [min, max)
            BaselineRMS();
            GetAmpMax();
            fWinMin_ = maxSample_ + int(offset/tUnit_) - lW;
            fWinMax_ = maxSample_ + int(offset/tUnit_) + hW;
            float deltaTPeakShift = -0.5 * (tmplTimeMaxScint_ + tmplTimeMaxSpike_);
            //---setup minimization
            ROOT::Math::Functor chi2(this, &WFClassLiTeDTU::TemplatesChi2, 4);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(100000);
            minimizer->SetMaxIterations(1000);
            minimizer->SetTolerance(1e-3);
            minimizer->SetPrintLevel(0);
            minimizer->SetFunction(chi2);
            minimizer->SetLimitedVariable(0, "amplitude_scint", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
            minimizer->SetLimitedVariable(1, "deltaT_scint", maxSample_*tUnit_+deltaTPeakShift, 1e-2, fWinMin_*tUnit_+deltaTPeakShift, fWinMax_*tUnit_+deltaTPeakShift);
            minimizer->SetLimitedVariable(2, "amplitude_spike", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
            minimizer->SetLimitedVariable(3, "deltaT_scint_minus_deltaT_spike", tmplTimeMaxScint_ - tmplTimeMaxSpike_, 1e-2, 0., fWinMax_*tUnit_);
            //---fit
            tmplFitConverged_ = minimizer->Minimize();

            //---try a second fit from a different starting point if the previous attempt failed
            if (not tmplFitConverged_) {
                minimizer->SetVariableValue(0, GetAmpMax()/2.);
                minimizer->SetVariableValue(1, maxSample_*tUnit_+deltaTPeakShift);
                minimizer->SetVariableValue(2, GetAmpMax()/2.);
                minimizer->SetVariableValue(3, tmplTimeMaxScint_ - tmplTimeMaxSpike_);
                tmplFitConverged_ = minimizer->Minimize();
            }

            const auto minparams = minimizer->X();
            tmplFitAmpScint_ = minparams[0];
            tmplFitTimeScint_ = minparams[1];
            tmplFitAmpSpike_ = minparams[2];
            tmplFitTimeSpike_ = minparams[1] - minparams[3];

            delete minimizer;
        } else { // fall back to interpolated values without spike amplitude
            const auto interpolFitResult = GetInterpolatedAmpMax();
            tmplFitAmpScint_ = interpolFitResult.ampl;
            tmplFitTimeScint_ = interpolFitResult.time;
            tmplFitAmpSpike_ = 0.;
            tmplFitTimeSpike_ = tmplFitTimeScint_;
            tmplFitConverged_ = false;
        }
    }

    return WFFitResultsScintPlusSpike{tmplFitAmpScint_, tmplFitTimeScint_, tmplFitAmpSpike_, tmplFitTimeSpike_, TemplatesChi2()/(fWinMax_-fWinMin_-4), tmplFitConverged_};
}

//----------Set the scintillation plus spike fit templates-------------------------------
void WFClassLiTeDTU::SetTemplateScint(TH1* templateWF)
{
    //---check input
    if(!templateWF)
    {
        cout << ">>>ERROR: scintillation template passed as input does not exist" << endl;
        return;
    }

    //---reset template fit variables
    if(interpolatorScint_)
        return;

    interpolatorScint_ = new ROOT::Math::Interpolator(0, ROOT::Math::Interpolation::kCSPLINE);
    tmplFitTimeScint_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tmplTimeMaxScint_ = tmplFitTimeScint_;
    tmplFitAmpScint_ = -1;

    //---fill interpolator data
    vector<double> x, y;
    for(int iBin=1; iBin<=templateWF->GetNbinsX(); ++iBin)
    {
        x.push_back(templateWF->GetBinCenter(iBin)-tmplFitTimeScint_);
        y.push_back(templateWF->GetBinContent(iBin));
    }
    interpolatorScint_->SetData(x, y);

    return;
}

void WFClassLiTeDTU::SetTemplateSpike(TH1* templateWF)
{
    //---check input
    if(!templateWF)
    {
        cout << ">>>ERROR: spike template passed as input does not exist" << endl;
        return;
    }

    //---reset template fit variables
    if(interpolatorSpike_)
        return;

    interpolatorSpike_ = new ROOT::Math::Interpolator(0, ROOT::Math::Interpolation::kCSPLINE);
    tmplFitTimeSpike_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tmplTimeMaxSpike_ = tmplFitTimeSpike_;
    //tmplFitTimeSpike_ = tmplFitTimeScint_;
    tmplFitAmpSpike_ = -1;

    ////Pulse shape for testing
    //TF1* spikeShape = new TF1("spikeShape", "1/0.999983*(exp(-0.5*((x)/5.1)^2) - 0.072*exp(-0.5*((x-20)/4.9)^2))", 0, 1000);
    //TF1* spikeShape = new TF1("spikeShape", "0", 0, 1000);

    //---fill interpolator data
    vector<double> x, y;
    for(int iBin=1; iBin<=templateWF->GetNbinsX(); ++iBin)
    {
        x.push_back(templateWF->GetBinCenter(iBin)-tmplFitTimeSpike_);
        y.push_back(templateWF->GetBinContent(iBin));
        //y.push_back(spikeShape->Eval(templateWF->GetBinCenter(iBin)-tmplFitTimeSpike_));
    }
    interpolatorSpike_->SetData(x, y);

    return;
}


//----------Reset: get new set of sample, keep interpolator-------------------------------
void WFClassLiTeDTU::Reset()
{
    startIndexCell_=-1;
    trigRef_=-1;
    sWinMin_=-1;
    sWinMax_=-1;
    sIntWinMin_=-1;
    sIntWinMax_=-1;
    bWinMin_=-1;
    bWinMax_=-1;
    bIntWinMin_=-1;
    bIntWinMax_=-1;
    maxSample_=-1;
    fitAmpMax_=-1;
    fitTimeMax_=-1;
    fitChi2Max_=-1;
    baseline_=-1;
    bRMS_=-1;
    cfSample_=-1;
    cfFrac_=-1;
    cfTime_=-1;
    cfSlope_=-1;
    leSample_=-1;
    leThr_=-1;
    leTime_=-1;
    leSlope_=-1;
    teSample_=-1;
    teThr_=-1;
    teTime_=-1;
    teSlope_=-1;
    chi2cf_=-1;
    chi2le_=-1;
    chi2te_=-1;
    fWinMin_=-1;
    fWinMax_=-1;
    tmplFitTime_=-1;
    tmplFitTimeErr_=-1;
    tmplFitAmp_=-1;
    tmplFitAmpShift_=0;
    tmplFitTimeScint_ = -1;
    tmplFitAmpScint_ = -1;
    tmplFitTimeSpike_ = -1;
    tmplFitAmpSpike_ = -1;
    tmplFitConverged_ = false;
    //interpolatorMin_=-1;
    //interpolatorMax_=-1;
    uncalibSamples_.clear();
    calibSamples_.clear();
    gain_.clear();
    times_.clear();

    samples_ = uncalibSamples_;
} 


//----------chi2 for scintillation plus spike template fit---------------------------------------------------------
double WFClassLiTeDTU::TemplatesChi2(const double* par)
{
    double chi2 = 0;
    double delta = 0;
    for(int iSample=fWinMin_; iSample<=fWinMax_; ++iSample)
    {
        if(iSample < 0 || iSample >= int(samples_.size()))
        {
            cout << ">>>WARNING: template fit out of samples range (chi2 set to -1)" << endl;
            chi2 += 9999;
        }
        else
        {
            //---fit: par[0]*ref_shape_scint(t-par[1]) + par[2]*ref_shape_spike(t-par[3])
            //---with par[0]=amplitude scintillation, par[1]=DeltaT scintillation
            //---and  par[2]=amplitude spike, par[3]=DeltaT spike
            //---if not fitting return chi2 value of best fit
            const auto iSampleTime = iSample*tUnit_;
            if(par)
                //delta = (samples_[iSample] - par[0]*interpolatorScint_->Eval(iSampleTime-par[1]) - par[2]*interpolatorSpike_->Eval(iSampleTime-par[3]))/bRMS_;
                delta = (samples_[iSample] - par[0]*interpolatorScint_->Eval(iSampleTime-par[1]) - par[2]*interpolatorSpike_->Eval(iSampleTime-(par[1]-par[3])))/bRMS_;
            else
                delta = (samples_[iSample] - tmplFitAmpScint_*interpolatorScint_->Eval(iSampleTime-tmplFitTimeScint_) - tmplFitAmpSpike_*interpolatorSpike_->Eval(iSampleTime-tmplFitTimeSpike_))/bRMS_;
            chi2 += delta*delta;
        }
    }

    return chi2;
}
//**********operators*********************************************************************
//----------assignment--------------------------------------------------------------------
WFClassLiTeDTU& WFClassLiTeDTU::operator=(const WFClassLiTeDTU& origin)
{
    uncalibSamples_ = origin.uncalibSamples_;
    calibSamples_ = origin.calibSamples_;
    if(&origin.samples_ == &origin.uncalibSamples_)
        samples_ = uncalibSamples_;
    else
        samples_ = calibSamples_;
    tUnit_ = origin.tUnit_;
    polarity_ = origin.polarity_;
    sWinMin_ = origin.sWinMin_;
    sWinMax_ = origin.sWinMax_;
    bWinMin_ = origin.bWinMin_;
    bWinMax_ = origin.bWinMax_;
    maxSample_ = origin.maxSample_;
    fitAmpMax_ = origin.fitAmpMax_;
    baseline_ = origin.baseline_;
    bRMS_ = origin.bRMS_;
    cfSample_ = origin.cfSample_;
    cfFrac_ = origin.cfFrac_;
    cfTime_ = origin.cfTime_;
    leSample_ = origin.leSample_;
    leThr_ = origin.leThr_;
    leTime_ = origin.leTime_;
    chi2cf_ = origin.chi2cf_;
    chi2le_ = origin.chi2le_;
    fWinMin_ = origin.fWinMin_;
    fWinMax_ = origin.fWinMax_;
    tmplFitTime_ = origin.tmplFitTime_;
    tmplFitTimeErr_ = origin.tmplFitTimeErr_;
    tmplFitAmp_ = origin.tmplFitAmp_;
    tmplFitAmpShift_ = origin.tmplFitAmpShift_;    
    interpolator_ = NULL;
    interpolatorScint_ = NULL;
    interpolatorSpike_ = NULL;

    return *this;
}
