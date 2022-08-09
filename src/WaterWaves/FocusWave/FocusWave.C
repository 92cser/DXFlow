/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 Dalian Ocean University 
						    National Marine Environment Monitoring Center 
							Dalian University of Technology
    Copyright (C) 2018-2022 Ocean University of China
	Copyright (C) 2022-     Ningbo University 
-------------------------------------------------------------------------------
License
    This file is part of DXFlow, a toolbox developed based on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "FocusWave.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(FocusWave, 0);
    addToRunTimeSelectionTable(waterWaves, FocusWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

Foam::waterWaveModels::FocusWave::FocusWave
(	
    const IOdictionary& dict
)
:
    SuperimposeWave(dict),
    amplitude_(readScalar(WWParameters_.lookup("FocusAmplitude"))),
    f_start_(readScalar(WWParameters_.lookup("Fre_start"))),
    f_end_(readScalar(WWParameters_.lookup("Fre_end"))),
    f_peak_(WWParameters_.lookupOrDefault<scalar>("Fre_peak", 0.5*(f_start_+f_end_))),
    t_focus_(readScalar(WWParameters_.lookup("FocusTime"))),
    x_focus_(readScalar(WWParameters_.lookup("FocusPosition"))),
    AmplitudeDiscretization_(WWParameters_.lookup("AmplitudeDiscretization")),
    TF_final_super_(NComponents_),
    TF_final_sub_(NComponents_)
{
    updateCoeff();
    if(SecondOrderCorrection_)
    {
        storeTFs();
    }
}

Foam::word Foam::waterWaveModels::FocusWave::name()
{
    return "Focused wave";
}

void Foam::waterWaveModels::FocusWave::updateCoeff()
{
    Cal_TF_array();
    Cal_LkC_array();
    Cal_A_array();
    Cal_Tf_array();
}

void Foam::waterWaveModels::FocusWave::Cal_TF_array()
{
    scalar deltaF = (f_end_-f_start_)/(NComponents_-1);
    forAll(T_array_, i)
    {
            F_array_[i] = f_start_ + i*deltaF;
            T_array_[i] = 1/F_array_[i];
    }
}

void Foam::waterWaveModels::FocusWave::Cal_LkC_array()
{
    forAll(L_array_, i)
    {
        scalar mu0_ = pow(2*pi_,2)*waterDepth_/g_*pow(T_array_[i],-2);
        scalar index0_ = 1.835+1.225*pow(mu0_,1.35);
        scalar mu_ = mu0_*(1+mu0_*exp(-index0_))*pow(tanh(mu0_),-0.5);
        // another iteration to improve precise
        scalar muFinal_ = (pow(mu_,2)+mu0_*pow(cosh(mu_),2))/(mu_+0.5*sinh(2*mu_));

        L_array_[i] = 2*pi_*waterDepth_/muFinal_;
        k_array_[i] = 2*pi_/L_array_[i];
        C_array_[i] = L_array_[i]/T_array_[i];
    }
}

void Foam::waterWaveModels::FocusWave::Cal_A_array() 
{
    if(AmplitudeDiscretization_ == "CWA")
    {
        scalar AveragedAmp = amplitude_/NComponents_;
        forAll(A_array_,i)
        {
             A_array_[i] = AveragedAmp;
        }
    }
    else if(AmplitudeDiscretization_ == "CWS")
    {
        scalar k_sum = 0;
        forAll(k_array_,i)
        {
                k_sum += 1/k_array_[i];
        }
        Info << "k_reverse_sum = " << k_sum << endl;
        forAll(A_array_,i)
        {
                A_array_[i] = amplitude_/k_array_[i]/k_sum;
        }
        // represent wave length
        scalar Lr_ = 2*pi_*amplitude_/k_array_[0]/A_array_[0]/NComponents_;
        Info << "The represented wave length is " << Lr_ << endl;
    }
    else if(AmplitudeDiscretization_ == "NewWaveJS")
    {
        // spectrum elevation factor
        scalar gamma = WWParameters_.lookupOrDefault<scalar>("GammaForSpectrum", 3.3);

        // calculate wave spectrum energy
        for(int i = 0; i < NComponents_; i++)
        {
            WS_[i] = Jonswap_(amplitude_, 1/f_peak_, gamma, F_array_[i]);
        }
        // calculate gross energy
        scalar GrossEnergy = 0;
        scalar deltaF = (f_end_-f_start_)/(NComponents_-1);
        for(int i = 0; i < NComponents_; i++)
        {
            GrossEnergy += WS_[i]*deltaF;
        }
        // calculate amplitude for each component
        for(int i = 0; i < NComponents_; i++)
        {
            A_array_[i] = amplitude_*WS_[i]*deltaF/GrossEnergy;
        }
    }
    else if(AmplitudeDiscretization_ == "NewWavePM")
    {
        // spectrum elevation factor
        scalar gamma = WWParameters_.lookupOrDefault<scalar>("GammaForSpectrum", 3.3);

        // calculate wave spectrum energy
        for(int i = 0; i < NComponents_; i++)
        {
            WS_[i] = PM_(f_peak_, F_array_[i]);
        }
        // calculate gross energy
        scalar GrossEnergy = 0;
        scalar deltaF = (f_end_-f_start_)/(NComponents_-1);
        for(int i = 0; i < NComponents_; i++)
        {
            GrossEnergy += WS_[i]*deltaF;
        }
        // calculate amplitude for each component
        for(int i = 0; i < NComponents_; i++)
        {
            A_array_[i] = amplitude_*WS_[i]*deltaF/GrossEnergy;
        }
    }
    else
    {
        FatalErrorInFunction
        << "Only CWA(constant wave amplitude), "<< nl
        << "	 CWS(constant wave steepness), "<< nl
        << "	 NewWave(JONSWAP or PM spectrum) are valid for focus wave. "
        << exit(FatalError);
    }
    // independant variable for curve fitting formular
    List<scalar> theta(NComponents_, 0.0);

    // enhance the amplitude
    forAll(Ef_array_,i)
    {
        scalar relativeWaterDepth = waterDepth_/L_array_[i];
        scalar nonlinearity = 2*A_array_[i]/waterDepth_;

        // theta value
        theta[i] = (1 - (1 - exp(4*nonlinearity) * nonlinearity)*relativeWaterDepth)*relativeWaterDepth;

        scalar thetaLimit = 0.5/(1 - exp(4*nonlinearity) * nonlinearity);

        if(relativeWaterDepth > thetaLimit)// deep water waves, theta decreases with the increase of D/L
        {
            // deep water approximation: theta is calculated by interpolation
            theta[i] = theta[i-1] + 0.2*(mag(theta[i]-theta[i-1]) + mag(theta[i-1]-theta[i-2]) + mag(theta[i-2]-theta[i-3]) + mag(theta[i-3]-theta[i-4]) + mag(theta[i-4]-theta[i-5]));
        }

        // enhanced factor
        Ef_array_[i] = 183.88*pow3(theta[i]) - 63.26*pow(theta[i],2) + 7.58*theta[i] + 0.72;

        // the minimum value is 1.0, maximum value is 3.115
        Ef_array_[i] = min(max(1.0, Ef_array_[i]), 3.115);

        // store and enhance
        A_EF_array_[i] = A_array_[i];
        A_EF_array_[i] *= Ef_array_[i];
        //Info << "D/L = " << relativeWaterDepth << ", enhanced factor = " << Ef_array_[i] << ", " << theta[i] << endl;
    }
}

void Foam::waterWaveModels::FocusWave::Cal_Tf_array()
{
    forAll(Tf_array_,i)
    {
        Tf_array_[i] = 4.0*sinh(k_array_[i]*waterDepth_)*sinh(k_array_[i]*waterDepth_)/(2.0*k_array_[i]*waterDepth_+sinh(2.0*k_array_[i]*waterDepth_));
    }
}

// For a focus wave, the phase velocities have been
// embedded in the calculation of surface elevation
Foam::scalar Foam::waterWaveModels::FocusWave::Cal_C()
{
    return 1.0;
}

Foam::scalar Foam::waterWaveModels::FocusWave::eta(scalar t)
{
    scalar eta_ = 0.0;
    // array stores first order surface elevation (contain phase velocity)
    List<scalar> FirstOrderEta(NComponents_, 0.0);

    forAll(FirstOrderEta,i)
    {
        FirstOrderEta[i] = A_EF_array_[i]*(cos(-k_array_[i]*x_focus_ - 2*pi_*F_array_[i]*(t-t_focus_)));
        eta_ += C_array_[i]*FirstOrderEta[i];
    }
    return eta_;
}

Foam::scalar Foam::waterWaveModels::FocusWave::eta1st(scalar t)
{
    scalar eta_ = 0.0;
    // array stores first order surface elevation (contain phase velocity)
    List<scalar> FirstOrderEta(NComponents_, 0.0);

    forAll(FirstOrderEta,i)
    {
        FirstOrderEta[i] = A_EF_array_[i]*(cos(-k_array_[i]*x_focus_ - 2*pi_*F_array_[i]*(t-t_focus_)));
        eta_ += FirstOrderEta[i];
    }
    return eta_;
}

Foam::scalar Foam::waterWaveModels::FocusWave::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
    position_ = 0.0;
    scalar secondOrderCorrection = 0.0;
    if(SecondOrderCorrection_)
    {
        forAll(TF_final_super_, i)
        {
            // get the reference
            List<complex>& curr_TF_final_super = TF_final_super_[i];
            List<complex>& curr_TF_final_sub = TF_final_sub_[i];
            //Info << "1" << endl;
            for(label j=i;j<NComponents_;j++)
            {
                // label for the subList
                const label lfs = j - i;
                const complex& curr_TF_super = curr_TF_final_super[lfs];
                const complex& curr_TF_sub = curr_TF_final_sub[lfs];
                waveComponent wave1(2*A_array_[i],waterDepth_,T_array_[i], true);
                waveComponent wave2(2*A_array_[j],waterDepth_,T_array_[j], true);
                SecondOrderWave CombinedWave(wave1,wave2,"lightConstruction");
                secondOrderCorrection += (CombinedWave.X2_super(0, x_focus_, t_focus_, t)*curr_TF_super).Re() + (CombinedWave.X2_sub(0, x_focus_, t_focus_, t)*curr_TF_sub).Re();
            }
        }
    }
    //first order signal
    forAll(A_array_,i)
    {
        scalar phi = k_array_[i]*(-x_focus_)-2*pi_*F_array_[i]*(t-t_focus_);
        position_ += A_array_[i]*sin(-phi)/Tf_array_[i];
    }
    scalar t_soft = min(0.5, 0.1*t_focus_);
    scalar SSF = (t<=t_soft)? pow(t/t_soft,2):1.0;
    return (position_ + secondOrderCorrection)*SSF + Cal_X_corr(t, deltaT, WaterLevel);
}

void Foam::waterWaveModels::FocusWave::PrintWaveProperties()
{
    Info << "-----------------------------" << "\n"
         << "Wave properties of Focus wave" << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Focus Ampitude (m) = " << amplitude_ << "\n"
         << " F_start (Hz) = " << f_start_ << "\n"
         << " F_end (Hz) = " << f_end_ << "\n"
         << " F_peak (Hz) = " << f_peak_ << "\n"
         << " Focus Postion (m) = " << x_focus_ << "\n"
         << " Focus Time (s) = " << t_focus_ << "\n"
         << " Amplitude Discretization = " << AmplitudeDiscretization_ << "\n"
         << " Active absorption = " << activeAbsorption_ << "\n"
         << "-----------------------------" << endl;
        if(AmplitudeDiscretization_=="NewWaveJS")
	{
                Info << " Target energy spectrum (Jonswap):" << endl;
		forAll(F_array_, i)
		{
                        Info << F_array_[i]<< ", " << WS_[i] << endl;
		}
	}
        else if(AmplitudeDiscretization_=="NewWavePM")
        {
                Info << " Target energy spectrum (PM):" << endl;
                forAll(F_array_, i)
                {
                        Info << F_array_[i]<< ", " << WS_[i] << endl;
                }
        }
}

void Foam::waterWaveModels::FocusWave::storeTFs()
{
    Info << "Calculating all the transfer functions for the focused wave..." << endl;
    Info << "The total number of components = " << NComponents_ << endl;
    forAll(TF_final_super_, i)
    {
        // get the reference
        List<complex>& curr_TF_final_super = TF_final_super_[i];
        List<complex>& curr_TF_final_sub = TF_final_sub_[i];
        // set size
        curr_TF_final_super.setSize(NComponents_ - i);
        curr_TF_final_sub.setSize(NComponents_ - i);
        // calculate
        for(label j=i;j<NComponents_;j++)
        {
            waveComponent wave1(2*A_array_[i],waterDepth_,T_array_[i]);
            waveComponent wave2(2*A_array_[j],waterDepth_,T_array_[j]);
            // label for the subList
            const label lfs = j - i;
            // update
            SecondOrderWave CombinedWave(wave1,wave2);
            curr_TF_final_super[lfs] = CombinedWave.TF_final_super();
            curr_TF_final_sub[lfs] = CombinedWave.TF_final_sub();
        }
    }
}

Foam::scalar Foam::waterWaveModels::FocusWave::Cal_X_corr(scalar t, scalar deltaT, scalar WaterLevel)
{
    if(activeAbsorption_)
    {
        if(t > 0.2*t_focus_)
        {
            scalar correctWaterLevel = eta1st(t) - (WaterLevel- waterDepth_);
            scalar U_corr_ = pow(g_/waterDepth_, 0.5)*correctWaterLevel;
            X_corr_ += U_corr_*deltaT;
        }
    }
    return X_corr_;
}
