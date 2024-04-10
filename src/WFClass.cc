#include "interface/WFClass.h"

#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TGraphErrors.h"
#include "TROOT.h"

//**********Constructors******************************************************************
WFClass::WFClass(int polarity, float tUnit, DigiChannelCalibration* calibration):
    uncalibSamples_(), calibSamples_(), samples_(uncalibSamples_),
    tUnit_(tUnit), polarity_(polarity), trigRef_(0), sWinMin_(-1), sWinMax_(-1), 
    bWinMin_(-1), bWinMax_(-1),  maxSample_(-1), fitAmpMax_(-1), fitTimeMax_(-1),
    fitChi2Max_(-1), baseline_(-1), bRMS_(-1), cfSample_(-1), cfFrac_(-1), cfTime_(-1),
    leSample_(-1), leTime_(-1), chi2cf_(-1), chi2le_(-1),
    fWinMin_(-1), fWinMax_(-1), tmplFitTime_(-1), tmplFitTimeErr_(-1), tmplFitAmp_(-1), tmplFitAmpShift_(0),
    f_max_(NULL), f_fit_(NULL), interpolator_(NULL)
{
    calibration_ = calibration;
    gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");
}
//**********Getters***********************************************************************

//----------Get the max/min amplitude wrt polarity----------------------------------------
float WFClass::GetAmpMax(int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return -1;
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---return the max if already computed
    else if(maxSample_ != -1)
        return samples_.at(maxSample_);
    // FIXME
    else if(samples_.size() == 0)
        return -1;
    
    //---find the max
    maxSample_=sWinMin_;
    for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        if(samples_.at(iSample) > samples_.at(maxSample_)) 
            maxSample_ = iSample;
    }    
    return samples_.at(maxSample_);
}

WFFitResults WFClass::GetInterpolatedSample(int sample, int samplesLeft, int samplesRight)
{
    if (sample>0 && ( samplesLeft>0 || samplesRight>0 ) )
    {
        float A,B, chi2;
        chi2 = LinearInterpolation(A, B, sample-samplesLeft, sample+samplesRight, sample); //get interpolated sample from left & right samples
        return WFFitResults{ A, 0, bRMS_, chi2, B};
    }
    else
        return {-1, -1000, -1, 0};
}

//----------Get the interpolated max/min amplitude wrt polarity---------------------------
WFFitResults WFClass::GetInterpolatedAmpMax(int min, int max, int nmFitSamples, int npFitSamples, string function, vector<float> params)
{
    //---check if already computed
    if(min==-1 && max==-1 && fitAmpMax_!=-1)
        return WFFitResults{fitAmpMax_, fitTimeMax_, -1, fitChi2Max_, 0};
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return {-1, -1000, -1, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---get maximum sample
    else if(maxSample_ == -1) 
        GetAmpMax(min, max); 

    //---fit the max
    double times[nmFitSamples+npFitSamples+1];
    double samples[nmFitSamples+npFitSamples+1];
    double xerr[nmFitSamples+npFitSamples+1];
    double yerr[nmFitSamples+npFitSamples+1];

    if(f_max_==NULL)
        f_max_ = new TF1("f_max", function.c_str(), 0., tUnit_*samples_.size(), TF1::EAddToList::kNo);
    else
        f_max_->SetRange(times_[maxSample_-nmFitSamples]-0.5*tUnit_, times_[maxSample_+npFitSamples]+0.5*tUnit_);
    for(unsigned int ipar = 0; ipar < params.size(); ++ipar)
        f_max_->SetParameter(ipar, params.at(ipar));
    
    int bin=0;
    double brms=BaselineRMS();
    for(int iSample=std::max(0, maxSample_-nmFitSamples); iSample<=std::min(GetNSample()-1, maxSample_+npFitSamples); ++iSample)
    {
        if (iSample < times_.size() && iSample < samples_.size()) {
            times[bin] = times_[iSample];
            samples[bin] = samples_[iSample];
            xerr[bin] = 0.;
            yerr[bin] = brms;
            ++bin;
        }
    }
    
    TGraphErrors h_max(nmFitSamples+npFitSamples+1, times, samples, xerr, yerr);
    
    if(h_max.GetMaximum() != 0)
    {
        auto fit_result = h_max.Fit(f_max_, "QRSO");
        fitTimeMax_ = f_max_->GetMaximumX();
        fitAmpMax_ = f_max_->GetMaximum();
        fitChi2Max_ = f_max_->GetChisquare() / f_max_->GetNDF();
    }
    else
    {
        fitTimeMax_ = -1;
        fitAmpMax_ = 1000;
        fitChi2Max_ = -1;
    }

    return WFFitResults{fitAmpMax_, fitTimeMax_, bRMS_, fitChi2Max_, 0};
}

//----------Get time with the specified method--------------------------------------------
WFFitResults WFClass::GetTime(string method, vector<float>& params)
{
    //---CFD
    if(method.find("CFD") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR, WFClass: to few arguments passed for CFD time computation" << endl;
        else if(params.size()<2)
            return GetTimeCF(params[0]);
        else if(params.size()<3)
            return GetTimeCF(params[0], params[1]);
        else
            return GetTimeCF(params[0], params[1], params[2], params[3]);
    }
    
    //---LED
    else if(method.find("LED") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR, WFClass: to few arguments passed for LED time computation" << endl;
        else if(params.size()<2)
            return GetTimeLE(params[0]);
        else if(params.size()<4)
            return GetTimeLE(params[0], params[1], params[2]);
        else
            return GetTimeLE(params[0], params[1], params[2], params[3], params[4]);
    }
    
    //---TED
    else if(method.find("TED") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR, WFClass: to few arguments passed for TED time computation" << endl;
        else if(params.size()<2)
            return GetTimeTE(params[0]);
        else if(params.size()<4)
            return GetTimeTE(params[0], params[1], params[2]);
        else
            return GetTimeTE(params[0], params[1], params[2], params[3], params[4]);
    }
    
    cout << ">>>ERROR, WFClass: time reconstruction method <" << method << "> not supported" << endl;
    return WFFitResults{-1, -1000, -1, -1, 0};
}

//----------Get CF time for a given fraction and in a given range-------------------------
WFFitResults WFClass::GetTimeCF(float frac, int nFitSamples, int min, int max)
{
    float A=0, B=0;
    if(frac != cfFrac_ || cfSample_ != -1)
    {
        //---setups---
        int tStart=min;
        if(tStart == -1)
            tStart=sWinMin_ == -1 ? 0 : sWinMin_;
        cfSample_ = tStart;
        cfFrac_ = frac;
        if(fitAmpMax_ == -1)
            GetInterpolatedAmpMax(min, max);
        if(frac == 1) 
            return WFFitResults{fitAmpMax_, times_[maxSample_], 1, 0};
        
        //---find first sample above Amax*frac
        for(int iSample=maxSample_; iSample>tStart; --iSample)
        {
            if(samples_.at(iSample) < fitAmpMax_*frac) 
            {
                cfSample_ = iSample;
                break;
            }
        }
        //---interpolate -- A+Bx = frac * amp
        chi2cf_ = LinearInterpolation(A, B, cfSample_-(nFitSamples-1)/2, cfSample_+(nFitSamples-1)/2);
        cfTime_ = (fitAmpMax_ * frac - A) / B;
        cfSlope_ = B;
    }

    return WFFitResults{fitAmpMax_*frac, cfTime_, -1, chi2cf_, cfSlope_};
}

//----------Get leading edge time at a given threshold and in a given range---------------
WFFitResults WFClass::GetTimeLE(float thr, int nmFitSamples, int npFitSamples, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{leThr_, -1000, -1, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---compute LED time value 
    float A=0, B=0;
    if(thr != leThr_ || leSample_ == -1)
    {
        //---find first sample above thr
        leThr_ = thr;
        leSample_ = -1;
        for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
        {
            if(samples_.at(iSample) > leThr_) 
            {
                leSample_ = iSample;
                break;
            }
        }
        if(leSample_>0)
        {
            //---interpolate -- A+Bx = amp
            chi2le_ = LinearInterpolation(A, B, leSample_-nmFitSamples, leSample_+npFitSamples);
            leTime_ = (leThr_ - A) / B;
            leSlope_ = B;
        }
        else
        {
            leTime_ = -1000;
            chi2le_ = -1;
            leSlope_ = 0;
        }
    }

    return WFFitResults{leThr_, leTime_, -1, chi2le_, leSlope_};
}

//----------Get trailing edge time at a given threshold and in a given range---------------
WFFitResults WFClass::GetTimeTE(float thr, int nmFitSamples, int npFitSamples, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{teThr_, -1000, -1, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---compute TED time value 
    float A=0, B=0;
    if(thr != teThr_ || teSample_ == -1)
    {
        //---find first sample above thr
        teThr_ = thr;
        teSample_ = -1;
        for(int iSample=maxSample_; iSample<samples_.size(); ++iSample)
        {
            if(samples_.at(iSample) < teThr_) 
            {
                teSample_ = iSample;
                break;
            }
        }
        if(teSample_>0)
        {
            //---interpolate -- A+Bx = amp
            chi2te_ = LinearInterpolation(A, B, teSample_-nmFitSamples, teSample_+npFitSamples);
            teTime_ = (teThr_ - A) / B;
            teSlope_ = B;
        }
        else
        {
            leTime_ = -1000;
            chi2le_ = -1;
            teSlope_ = 0;
        }
    }

    return WFFitResults{teThr_, teTime_, -1, chi2te_, teSlope_};
}

//----------Get the waveform integral in the given range----------------------------------
float WFClass::GetIntegral(int min, int max)
{
    //---compute integral
    float integral=0;
    for(int iSample=min; iSample<max; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        integral += samples_.at(iSample);
    }

    return integral;
}

//----------Get the signal integral around the the max-------------------------------------
float WFClass::GetSignalIntegral(int riseWin, int fallWin)
{
    //---compute position of the max
    if(maxSample_ == -1)
        GetAmpMax();

    //---compute integral
    float integral=0;
    for(int iSample=maxSample_-riseWin; iSample<maxSample_+fallWin; ++iSample)
    {
        //---if signal window goes out of bound return a bad value
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        integral += samples_.at(iSample);
    }

    return integral;
}


//----------Get the integral of Abs(WF) over the given range------------------------------
float WFClass::GetModIntegral(int min, int max)
{   
    float integral=0;
    for(int iSample=min; iSample<max; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        if(samples_.at(iSample) < 0)
            integral -= samples_.at(iSample);
        else
            integral += samples_.at(iSample);
    }
    return integral;
}

//**********Setters***********************************************************************

//----------Set the signal windows---------------------------------------------------------
void WFClass::SetSignalWindow(int min, int max)
{
    sWinMin_ = std::max(int(min + trigRef_), 0);
    sWinMax_ = std::min(int(max + trigRef_), int(samples_.size()));
    if (sWinMax_ < sWinMin_) {
        sWinMin_ = std::max(sWinMax_ - 1, 0);
    }
}

void WFClass::SetSignalIntegralWindow(int min, int max)
{
    sIntWinMin_ = min + trigRef_;
    sIntWinMax_ = max + trigRef_;
}

//----------Set the baseline windows-------------------------------------------------------
void WFClass::SetBaselineWindow(int min, int max)
{
    bWinMin_ = std::max(min, 0);
    bWinMax_ = std::min(max, int(samples_.size()));
    if (bWinMax_ < bWinMin_) {
        bWinMin_ = std::max(bWinMax_ - 1, 0);
    }
}

void WFClass::SetBaselineIntegralWindow(int min, int max)
{
    bIntWinMin_ = std::max(min, 0);
    bIntWinMax_ = std::min(max, int(samples_.size()));
    if (bIntWinMax_ < bIntWinMin_) {
        bIntWinMin_ = std::max(bIntWinMax_ - 1, 0);
    }
}

//----------Set the fit template----------------------------------------------------------
void WFClass::SetTemplate(TH1* templateWF)
{
    //---check input
    if(!templateWF)
    {
        cout << ">>>ERROR: template passed as input does not exist" << endl;
        return;
    }

    //---reset template fit variables
    if(interpolator_)
        delete interpolator_;

    interpolator_ = new ROOT::Math::Interpolator(0, ROOT::Math::Interpolation::kCSPLINE);
    tmplFitTime_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tmplFitAmp_ = -1;

    //---fill interpolator data
    vector<double> x, y;
    for(int iBin=1; iBin<=templateWF->GetNbinsX(); ++iBin)
    {
        x.push_back(templateWF->GetBinCenter(iBin)-tmplFitTime_);
        y.push_back(templateWF->GetBinContent(iBin));
    }
    interpolator_->SetData(x, y);
    interpolatorMin_ = templateWF->GetBinCenter(1)-tmplFitTime_;
    interpolatorMax_ = templateWF->GetBinCenter(templateWF->GetNbinsX())-tmplFitTime_;

    return;
}


//**********Utils*************************************************************************

//----------Reset: get new set of sample, keep interpolator-------------------------------
void WFClass::Reset()
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
    uncalibSamples_.clear();
    calibSamples_.clear();
    gain_.clear();
    times_.clear();

    samples_ = uncalibSamples_;
} 

//----------Apply channel Amp and Time intercalibration-----------------------------------
bool WFClass::ApplyCalibration()
{
    if(!calibration_)
    {
        std::cout << "[WFClass]::ERROR::No calibration available" << std::endl;
        return false;
    }

    //---calibration is applied per cellIndex (not sampleIndex)
    calibSamples_.resize(samples_.size());
    for (int iSample=0; iSample<samples_.size(); ++iSample)
    {
        int cellIndex=(iSample+startIndexCell_)%samples_.size();
        auto& cc = calibration_->at(cellIndex);
        auto us = polarity_*uncalibSamples_.at(iSample);
        times_.at(iSample) = iSample * tUnit_ - cc.deltaT; // time calibration
        calibSamples_.at(iSample) =
            polarity_*(us +                      // uncalib sample
                       cc.deltaV     +           // constant DeltaV
                       cc.slopeV     * us +      // linear term
                       cc.quadraticV * us * us); // quadratic term
    }

    samples_ = calibSamples_;
    
    return true;
}
//---------Add waveform sample to the list of uncalibrated samples------------------------
//---sample is inserted at the end of uncalibSamples_
//---the times vector is filled with the uncalibrated sample time computed from the time unit
void WFClass::AddSample(float sample)
{
    uncalibSamples_.push_back(polarity_*sample); 
    times_.push_back( (samples_.size()-1.)*tUnit_ );
    samples_ = uncalibSamples_;
};



//---------estimate the baseline in a given range and then subtract it from the signal----
WFBaseline WFClass::SubtractBaseline(int min, int max)
{
    if(min!=-1 && max==-1)
    {
        bWinMin_=min;
        bWinMax_=max;
    }
    //---compute baseline
    float baseline_=0;
    int nSamples=0;
    for(int iSample=bWinMin_; iSample<bWinMax_; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        baseline_ += samples_.at(iSample);
	++nSamples; 
    }
    baseline_ = baseline_/((float)nSamples);
    //---subtract baseline
    for(unsigned int iSample=0; iSample<samples_.size(); ++iSample)
        samples_.at(iSample) = (samples_.at(iSample) - baseline_);    
    //---interpolate baseline
    BaselineRMS();
    float A=0, B=0;
    float chi2 = LinearInterpolation(A, B, bWinMin_, bWinMax_);
    
    return WFBaseline{baseline_, bRMS_, A, B, chi2};
}

//---------subtract a custom baseline from the signal----
WFBaseline WFClass::SubtractBaseline(float baseline)
{
    //---subtract baseline
    for(unsigned int iSample=0; iSample<samples_.size(); ++iSample)
        samples_.at(iSample) = (samples_.at(iSample) - baseline);    
    //---interpolate baseline
    float A=0, B=0;
    float chi2 = LinearInterpolation(A, B, bWinMin_, bWinMax_);    
    return WFBaseline{baseline, -1, A, B, chi2};
}


//----------template fit to the WF--------------------------------------------------------
WFFitResults WFClass::TemplateFit(float amp_threshold, float offset, int lW, int hW)
{
    double tmplFitChi2=0;
    if(tmplFitAmp_ == -1)
    {
        if(samples_[maxSample_]>amp_threshold)
        {
            //---set template fit window around maximum, [min, max)
            BaselineRMS();
            GetAmpMax();    
            fWinMin_ = maxSample_ + int(offset/tUnit_) - lW;
            fWinMax_ = maxSample_ + int(offset/tUnit_) + hW;
            //---setup minimization
            auto t0 = GetInterpolatedAmpMax().time;
            ROOT::Math::Functor chi2(this, &WFClass::TemplateChi2, 2);
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

//----------analytic fit to the WF--------------------------------------------------------
double WFClass::AnalyticFit(TF1* f, int lW, int hW)
{
    f_fit_ = f;

    fWinMin_ = lW;
    fWinMax_ = hW;
    //---setup minimization
    ROOT::Math::Functor chi2(this, &WFClass::AnalyticChi2, f_fit_->GetNumberFreeParameters());
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(1000);
    minimizer->SetTolerance(1e-6);
    minimizer->SetPrintLevel(0);
    minimizer->SetFunction(chi2);
    for (int iPar=0;iPar<f_fit_->GetNumberFreeParameters();++iPar)
    {
        double pmin,pmax;
        f_fit_->GetParLimits(iPar,pmin,pmax);
        minimizer->SetLimitedVariable(iPar, Form("p%d",iPar), f_fit_->GetParameter(iPar), (pmax-pmin)/1000., pmin, pmax );
    }
    //---fit
    minimizer->Minimize();
    for (int iPar=0;iPar<f_fit_->GetNumberFreeParameters();++iPar)
    {
        f_fit_->SetParameter(iPar,minimizer->X()[iPar]);
        f_fit_->SetParError(iPar,sqrt(minimizer->CovMatrix(iPar,iPar)));
    }

    return minimizer->MinValue();
}

void WFClass::EmulatedWF(WFClass& wf,float rms, float amplitude, float time)
{
    TRandom3 rnd(0);

    wf.Reset();

    if (tmplFitTime_ == -1)
    {
        std::cout << "ERROR: no TEMPLATE for this WF" << std::endl;
        return;
    }

    for (unsigned int i=0; i<samples_.size();++i)
    {
        float emulatedSample=amplitude*interpolator_->Eval(i*tUnit_-tmplFitTime_-(time-tmplFitTime_));
        emulatedSample+=rnd.Gaus(0,rms);
        wf.AddSample(emulatedSample);
    }
}


void WFClass::FFT(WFClass& wf, float tau, int cut)
{
    if(samples_.size() == 0)
    {
        std::cout << "ERROR: EMPTY WF" << std::endl;
        return;
    }

    wf.Reset();

    int n=samples_.size();
    TVirtualFFT *vfft = TVirtualFFT::FFT(1,&n,"C2CFORWARD");

    Double_t orig_re[n],orig_im[n];
    for(int i=0;i<n;i++) 
    {
        orig_re[i]=samples_[i];
        if(i>1000) orig_re[i]=orig_re[999];// DIGI CAENV1742 NOT USABLE
        orig_im[i]=0;
    }
    vfft->SetPointsComplex(orig_re,orig_im);
    vfft->Transform();
    Double_t re[n],im[n];
    vfft->GetPointsComplex(re,im);

    TVirtualFFT *vinvfft = TVirtualFFT::FFT(1,&n,"C2CBACKWARD M K");
    Double_t cut_re[n],cut_im[n];

    for(int i=0;i<n;i++) 
    {
        if( i> cut-1 && i<n-cut) 
        {
            int delta = TMath::Min(i-cut-1,n-cut-i); 
            double dump=TMath::Exp(-delta/tau);
            cut_im[i]=im[i]*dump;
            cut_re[i]=re[i]*dump;
            continue;
        }
        cut_re[i]= re[i];
        cut_im[i]= im[i];
    }

    vinvfft->SetPointsComplex(cut_re,cut_im);
    vinvfft->Transform();
    Double_t inv_re[n],inv_im[n];
    vinvfft->GetPointsComplex(inv_re,inv_im);

    for(int i=0;i<n ;i++)
        wf.AddSample(inv_re[i]/n);

    delete vinvfft;
    delete vfft;

    return;
}

//----------compute baseline RMS (noise)--------------------------------------------------
float WFClass::BaselineRMS()
{
    if(bRMS_ != -1)
        return bRMS_;

    int nSample=0;
    float sum=0, sum2=0;
    for(int iSample=bWinMin_; iSample<bWinMax_; ++iSample)
    {
        ++nSample;
        sum += samples_.at(iSample);
        sum2 += samples_.at(iSample)*samples_.at(iSample);
    }
    
    bRMS_=sqrt(sum2/nSample - pow(sum/nSample, 2));
    return bRMS_;
}

//----------Linear interpolation util-----------------------------------------------------
float WFClass::LinearInterpolation(float& A, float& B, const int& min, const int& max, const int& sampleToSkip)
{
    TLinearFitter lf(1);
    lf.SetFormula("hyp1");
    int usedSamples=0;
    for(int iSample=min; iSample<=max; ++iSample)
    {
        if(iSample<0 || iSample>=int(samples_.size())) 
            continue;
        if (sampleToSkip>=0 && iSample==sampleToSkip)
            continue;
        double x=times_.at(iSample);
        lf.AddPoint(&x,samples_.at(iSample));
        ++usedSamples;
    }
  
    lf.Eval();
    A=lf.GetParameter(0);
    B=lf.GetParameter(1);
    float chi2=lf.GetChisquare();
    
    return chi2/(usedSamples-2);

}

//----------chi2 for template fit---------------------------------------------------------
double WFClass::TemplateChi2(const double* par)
{
    double chi2 = 0;
    double delta2 = 0;
#ifdef PARALLEL
#pragma omp parallel for reduction(+:chi2)
#endif    
    for(int iSample=fWinMin_; iSample<=fWinMax_; ++iSample)
    {
        if(iSample < 0 || iSample >= int(samples_.size()) ||
           (par && (times_.at(iSample)-par[1] < interpolatorMin_ || times_.at(iSample)-par[1] > interpolatorMax_) ) )
        {
            //cout << ">>>WARNING: template fit out of samples rage (chi2 set to 9999)" << endl;
            chi2 += 9999;
        }
        else
        {
            //---fit: par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
            //---if not fitting return chi2 value of best fit            
            if(par)
            {
                //auto deriv = par[0]*interpolator_->Deriv(times_[iSample]-par[1]);
                auto err2 = bRMS_*bRMS_;// + pow(tUnit_/sqrt(12)*deriv/2, 2);
                delta2 = pow((samples_.at(iSample) - par[0]*interpolator_->Eval(times_[iSample]-par[1])), 2)/err2;
            }
            else
            {
                //auto deriv = tmplFitAmp_*interpolator_->Deriv(times_[iSample]-tmplFitTime_);
                auto err2 = bRMS_*bRMS_;// + pow(tUnit_/sqrt(12)*deriv/2, 2);
                delta2 = pow((samples_.at(iSample) - tmplFitAmp_*interpolator_->Eval(times_[iSample]-tmplFitTime_)), 2)/err2;
            }
            chi2 += delta2;
        }
    }

    return chi2;
}

//----------chi2 for analytic fit---------------------------------------------------------
double WFClass::AnalyticChi2(const double* par)
{

    if (f_fit_ == NULL)
        return -1;
  
    double chi2 = 0;
    double delta = 0;

    if(par)
    {
        for (int iPar=0; iPar<f_fit_->GetNumberFreeParameters();++iPar)
            f_fit_->SetParameter(iPar,*(par+iPar));
    }
    
    for(int iSample=fWinMin_; iSample<=fWinMax_; ++iSample)
    {
        if(iSample < 0 || iSample >= int(samples_.size()))
        {
            cout << ">>>WARNING: analytic fit out of samples rage (chi2 set to -1)" << endl;
            chi2 += 9999;
        }
        else
        {
            //----IMPROVE-ME
            delta = (samples_.at(iSample) - f_fit_->Eval(times_.at(iSample)))/10.;
            chi2 += delta*delta;
        }
    }
      
    return chi2;
}



void WFClass::Print()
{
    std::cout << "+++ DUMP WF +++" << std::endl;
    for (unsigned int i=0; i<samples_.size(); ++i)
        std::cout << "SAMPLE " << i << ": " << samples_[i] << std::endl;
}

//**********operators*********************************************************************
//----------assignment--------------------------------------------------------------------
WFClass& WFClass::operator=(const WFClass& origin)
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
    return *this;
}

//----------subtraction-------------------------------------------------------------------
WFClass WFClass::operator-(const WFClass& sub)
{
    if(tUnit_ != sub.tUnit_)
        return *this;

    WFClass diff(1, tUnit_);
    for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
        diff.AddSample(samples_.at(iSample) - sub.samples_.at(iSample));

    return diff;
}

//----------addition----------------------------------------------------------------------
WFClass WFClass::operator+(const WFClass& sub)
{
    if(tUnit_ != sub.tUnit_)
        return *this;

    WFClass sum(1, tUnit_);
    for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
        sum.AddSample(samples_.at(iSample) + sub.samples_.at(iSample));

    return sum;
}

//----------subtraction and assignment----------------------------------------------------
WFClass& WFClass::operator-=(const WFClass& sub)
{
    if(tUnit_ == sub.tUnit_)
        for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
            samples_.at(iSample) -= sub.samples_.at(iSample);
    
    return *this;
}

//----------addition and assignmet--------------------------------------------------------
WFClass& WFClass::operator+=(const WFClass& sub)
{
    if(tUnit_ == sub.tUnit_)
        for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
            samples_.at(iSample) += sub.samples_.at(iSample);
    
    return *this;
}
