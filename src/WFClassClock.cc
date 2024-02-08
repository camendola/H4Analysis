#include "interface/WFClassClock.h"
#include "TError.h" 

//**********Constructors******************************************************************

WFClassClock::WFClassClock(int polarity,float tUnit):
    WFClass(polarity, tUnit), extPeriod_(-1), clkPeriod_(-1), tmplFitPeriod_(-1)
{}

//---------estimate the baseline from a template fit of the clock signal
WFBaseline WFClassClock::SubtractBaselineFit(int min, int max)
{
    if(min!=-1 && max==-1)
    {
        bWinMin_=min;
        bWinMax_=max;
    }
    //---compute baseline
    float baseline_=0;
    bRMS_=10; //custom value as starting point...

    //Remove annoying Minuit2 warnings
    gErrorIgnoreLevel = kError;

    //---setup minimization
    fWinMin_=bWinMin_;
    fWinMax_=bWinMax_;
    int maxSample_=bWinMin_;
    int minSample_=bWinMin_;
    for(int iSample=bWinMin_; iSample<bWinMax_; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        if(samples_.at(iSample) > samples_.at(maxSample_)) 
            maxSample_ = iSample;
        if(samples_.at(iSample) < samples_.at(minSample_)) 
            minSample_ = iSample;
    }    

    float baseline_mean=(samples_.at(maxSample_)+samples_.at(minSample_))/2.;
    std::default_random_engine generator;
    auto t0 = (times_.at(bWinMin_)+times_.at(bWinMax_))/2.;
    std::uniform_real_distribution<float> rnd_t0(-2,2);
    ROOT::Math::Functor chi2(this, &WFClassClock::TemplateChi2, 2);
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(1000);
    minimizer->SetTolerance(1e-4);
    minimizer->SetPrintLevel(-1);
    minimizer->SetFunction(chi2);
    minimizer->SetLimitedVariable(0, "deltaV", baseline_mean, 1e-2, baseline_mean*0.7,baseline_mean*1.3);
    minimizer->SetLimitedVariable(1, "deltaT", t0, 1e-3, t0-10, t0+10);
    //---fit
    minimizer->Minimize();
    int tries=0;
    while(minimizer->Status() != 0 && tries < 1000)
      {
	minimizer->SetVariableValue(1, (times_[bWinMin_]+times_[bWinMax_])/2.+rnd_t0(generator));
	minimizer->Minimize();
	++tries;
      }    
    baseline_ = minimizer->X()[0];
    tmplFitAmp_ = minimizer->X()[0];
    tmplFitTime_ = minimizer->X()[1];
    delete minimizer;        

    //---subtract baseline
    for(unsigned int iSample=0; iSample<samples_.size(); ++iSample)
        samples_.at(iSample) = (samples_.at(iSample) - baseline_);    
    //---interpolate baseline
    float A=0, B=0;
    float myChi2= float(TemplateChi2()/(bWinMax_-bWinMin_));
    bRMS_*=TMath::Sqrt(myChi2); //fermi method to estimate uncertainty...
    //restore the initialised fit conditions before return
    tmplFitAmp_ = -1;
    tmplFitTime_ = 0;
    return WFBaseline{baseline_, bRMS_, A, B, myChi2};
}

//**********Getters***********************************************************************

//----------Get time of a clock pulse-----------------------------------------------------
WFFitResults WFClassClock::GetTime(string method, vector<float>& params)
{
    //---CLK
    if(method.find("CLK") != string::npos)
    {
        if(params.size()==0)
            return GetTimeCLK();
        else if(params.size()==2)
            return GetTimeCLK(params[0], params[1]);
        else if(params.size()==3)
            return GetTimeCLK(params[0], params[1], params[2]);
        else if(params.size()==4)
	  return GetTimeCLK(params[0], params[1], params[2],params[3]);
        else
            cout << ">>>ERROR: wrong number of arguments passed for CLK time computation: " << params.size()
                 << ". usage: CLK [min max] [+/- 1]" << endl;

    }
    //---LED
    else if(method.find("LED") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR: to few arguments passed for LED time computation" << endl;
        else if(params.size()<2)
            return GetTimeLE(params[0]);
        else if(params.size()<4)
            return GetTimeLE(params[0], params[1], params[2]);
        else
            return GetTimeLE(params[0], params[1], params[2], params[3], params[4]);

    }
   
    cout << ">>>ERROR, WFClassClock: time reconstruction method <" << method << "> not supported" << endl;
    return WFFitResults{-1, -1000, -1, 0};
}

//----------Get time of a clock pulse-----------------------------------------------------
WFFitResults WFClassClock::GetTimeLE(float thr, int nmFitSamples, int npFitSamples, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{leThr_, -1000, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---compute LED time value 
    float A=0, B=0;
    if(thr != leThr_ || leSample_ == -1)
    {
        //---find first sample above thr on RISING edge
        leThr_ = thr;
        for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
        {
            if(samples_.at(iSample) > leThr_ && samples_.at(iSample) > samples_.at(iSample-1)) 
            {
                leSample_ = iSample;
                break;
            }
        }
        //---interpolate -- A+Bx = amp
        chi2le_ = LinearInterpolation(A, B, leSample_-nmFitSamples, leSample_+npFitSamples);
        leTime_ = (leThr_ - A) / B;
    }

    return WFFitResults{leThr_, leTime_, chi2le_, B};
}

//----------Get time of a clock pulse-----------------------------------------------------
//---Each rising edge zero crossing is measured, the clock period is
//---then extracted from a heuristic fit (cit Marc)
WFFitResults WFClassClock::GetTimeCLK(float wleft, float wright, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{-1, -1000, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);

    double S=0, S2=0, T=0, ST=0; 
    vector<double> clkZeros;
    TGraph gClockEdge(samples_.size(), times_.data(), samples_.data());
    TF1 fClockEdge("fClockEdge", "[0]*tanh((x-[1])/[2])+[3]", times_[min], times_[max]);
    
    for(int iSample=min; iSample<max; ++iSample)
    {
        if(samples_.at(iSample) < 0 && samples_.at(iSample+1) > 0)
        {
            //---edge fit
            auto t0 = times_[iSample];
            if(t0+wleft > times_[min] && t0+wright < times_[max])
            {
                fClockEdge.SetParameters(300., t0, 0.5, 0.);
                gClockEdge.Fit(&fClockEdge, "QRSO", "goff", t0+wleft, t0+wright);
                clkZeros.push_back(fClockEdge.GetParameter(1));            
                //---fill variables
                auto N = clkZeros.size();
                S += N;
                S2 += N*N;
                T += clkZeros.back();
                ST += N*clkZeros.back();
            
                ++iSample;
            }
        }
    }

    //---compute clock period, phase and phase error
    auto N = clkZeros.size();
    clkPeriod_ = (ST-S*T/N)/(S2-S*S/N);
    double phase = (T-S*clkPeriod_)/N;
    TH1F h_phase_err("h_phase_err", "", 75, -1, 1);
    for(unsigned int i=1; i<N; ++i)
    {
        auto t_clk = clkZeros[i]-clkZeros[0];
        h_phase_err.Fill(i*clkPeriod_-t_clk);
        //phase_err += (phase+i*period-t_clk)*(phase+i*period-t_clk);
    }
    //auto phase_err = h_phase_err.Fit("gaus", "QRSO", "goff")->Parameter(2);    
    auto phase_err = h_phase_err.GetRMS();
    
//    return WFFitResults{0, phase, std::sqrt(phase_err/(N-1)), 0};
    return WFFitResults{0, phase, phase_err, -1, 0};
}

//----------template fit to the WF--------------------------------------------------------
//---Template fit of each single clock cycle
WFFitResults WFClassClock::TemplateFit(float ampl_threshold, float offset, int lW, int hW)
{
    //Remove annoying Minuit2 warnings
    gErrorIgnoreLevel = kError;

    if(tmplFitAmp_ == -1)
    {
        // BaselineRMS();
        GetAmpMax();
        bRMS_ = 10.;
        int N=0;
        double S=0, S2=0, T=0, ST=0; 
        vector<double> clkZeros;
        double period = 0;
        //---cycles loop
        for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
        {
            //---single clock cycle fit
            if(samples_.at(iSample) < 0 && samples_.at(iSample+1) > 0)
            {
                auto t0 = (times_.at(iSample+1)+times_.at(iSample))/2.;
                //---set template fit window around zero crossing, [min, max)
                fWinMin_ = iSample - lW;
                fWinMax_ = iSample + hW;
                
                //---setup minimization
                std::default_random_engine generator;
                std::uniform_real_distribution<float> rnd_t0(times_.at(iSample), times_.at(iSample+1));
                ROOT::Math::Functor chi2(this, &WFClassClock::TemplateChi2, 2);
                ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
                minimizer->SetMaxFunctionCalls(100000);
                minimizer->SetMaxIterations(1000);
                minimizer->SetTolerance(1e-2);
                minimizer->SetPrintLevel(-1);
                minimizer->SetFunction(chi2);
                minimizer->SetLimitedVariable(0, "deltaV", 0, 1e-2, -0.1*GetAmpMax(), 0.1*GetAmpMax() );
                minimizer->SetLimitedVariable(1, "deltaT", t0, 1e-3, times_[fWinMin_]-1., times_[fWinMax_]+1);
                //---fit
                minimizer->Minimize();
                int tries=0;
                while(minimizer->Status() != 0 && tries < 1000)
                {
                    minimizer->SetVariableValue(1, rnd_t0(generator));
                    minimizer->Minimize();
                    ++tries;
                }    
                if(tmplFitAmp_ == -1)
                    tmplFitAmp_ = minimizer->X()[0];
                tmplFitTime_ = minimizer->X()[1];
                clkZeros.push_back(minimizer->X()[1]);
                //---fill variables
                ++N;
                S += N;
                S2 += N*N;
                T += clkZeros.back();
                ST += N*clkZeros.back();                
                ++iSample;                
                delete minimizer;        
            }
        }        
        //---compute clock period, phase and phase error
        tmplFitPeriod_ = (ST-S*T/N)/(S2-S*S/N);
        period = (extPeriod_ == -1 ? tmplFitPeriod_ : extPeriod_) ;  
        tmplFitTime_ = (T-S*period)/N;
        TH1F h_phase_err("h_phase_err", "", 100, -0.05, 0.05);
        for(unsigned int i=1; i<N; ++i)
        {
            auto t_clk = clkZeros[i]-clkZeros[0];
            h_phase_err.Fill(i*tmplFitPeriod_-t_clk);
        }
        tmplFitTimeErr_ = h_phase_err.GetRMS();
        //tmplFitTimeErr_ = h_phase_err.Fit("gaus", "QRSO", "goff")->Parameter(2);
    }
    return WFFitResults{tmplFitAmp_, tmplFitTime_, tmplFitTimeErr_, TemplateChi2()/(hW-lW), 0};
}

//---------Add waveform sample to the list of uncalibrated samples------------------------
//---sample is inserted at the end of uncalibSamples_
//---the times vector is filled with the uncalibrated sample time computed from the time unit
void WFClassClock::AddSample(float sample)
{
    uncalibSamples_.push_back(polarity_*sample); 
    times_.push_back( (samples_.size()-1.)*tUnit_ );
    samples_ = uncalibSamples_;
};

//----------Set the fit template----------------------------------------------------------
void WFClassClock::SetTemplate(TH1* templateWF)
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
    //tmplFitTime_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tmplFitTime_ = 0; 
    tmplFitAmp_  = -1;

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

double WFClassClock::TemplateChi2(const double* par)
{
    double chi2 = 0;
    double delta2 = 0;
#ifdef PARALLEL
#pragma omp parallel for reduction(+:chi2)
#endif    
    for(int iSample=fWinMin_; iSample<=fWinMax_; ++iSample)
    {
        if(iSample < 0 || iSample >= int(samples_.size()) 
                || (par && (times_.at(iSample)-par[1] < interpolatorMin_ || times_.at(iSample)-par[1] > interpolatorMax_) ) 
                )

        {
            //cout << ">>>WARNING: template fit out of samples rage (chi2 set to 9999)" << endl;
            //cout << "iSample " << iSample << endl;
            //cout << "polarity " << polarity_<< endl;
            //cout << " interpolatorMin_ "<< interpolatorMin_<< endl;
            //cout << " interpolatorMax_ "<< interpolatorMax_<< endl;
            //cout << " fWinMin_ "<< fWinMin_<< endl;
            //cout << " fWinMax_ "<< fWinMax_<< endl;
            //cout << "( times_.at(iSample) -par[1])" << times_.at(iSample)-par[1] << endl;
            //cout << "par[1] " << par[1] << endl;
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
                delta2 = pow((samples_.at(iSample) - par[0]+interpolator_->Eval(times_[iSample]-par[1])), 2)/err2;
            }
            else
            {
                //auto deriv = tmplFitAmp_*interpolator_->Deriv(times_[iSample]-tmplFitTime_);
                auto err2 = bRMS_*bRMS_;// + pow(tUnit_/sqrt(12)*deriv/2, 2);
                delta2 = pow((samples_.at(iSample) - tmplFitAmp_+interpolator_->Eval(times_[iSample]-tmplFitTime_)), 2)/err2;
            }
            chi2 += delta2;
        }
    }
    return chi2;
}
