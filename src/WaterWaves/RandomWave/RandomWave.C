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

#include "RandomWave.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace waterWaveModels
{
    defineTypeNameAndDebug(RandomWave, 0);
    addToRunTimeSelectionTable(waterWaves, RandomWave, nameList);
}// end of namespace waterWaveModels
}// end of namespace Foam

// Constructor
Foam::waterWaveModels::RandomWave::RandomWave
(	
	const IOdictionary& dict
)
:
    SuperimposeWave(dict),
    SpectrumType_(WWParameters_.lookupOrDefault<word>("SpectrumType","Jonswap")),
    waveHeight_sig_(readScalar(WWParameters_.lookup("SignificantWaveHeight"))),
    peakPeriod_(readScalar(WWParameters_.lookup("PeakPeriod"))),
    AmplitudeDiscretization_(WWParameters_.lookup("AmplitudeDiscretization")),
    activeAbsorption_(WWParameters_.lookupOrDefault<bool>("ActiveAbsorption",false)),
    peakFrequency_(1.0/peakPeriod_),
    gamma_(WWParameters_.lookupOrDefault<scalar>("GammaForSpectrum", 3.3)),
    f_start_(0.6*peakFrequency_),
    f_end_(3.0*peakFrequency_),
    RandomPhase_(T_array_),
    X_corr_(0.0),
    TF_final_super_(NComponents_),
    TF_final_sub_(NComponents_)
{
    updateCoeff();
    if(SecondOrderCorrection_)
    {
        storeTFs();
    }
}

void Foam::waterWaveModels::RandomWave::updateCoeff()
{
    Cal_TF_array();
    Cal_LkC_array();
    Cal_EF_array();
    Cal_randomPhase();
    Cal_Tf_array();
}


void Foam::waterWaveModels::RandomWave::Cal_TF_array()
{
    scalar averagedEnergy = 0.0;
    scalar deltaF = (f_end_-f_start_)/(NComponents_-1);
    forAll(T_array_, i)
    {
            F_array_[i] = f_start_ + i*deltaF;
            T_array_[i] = 1/F_array_[i];
    }
    if(AmplitudeDiscretization_ == "AE")  // averaged energy
    {
            // calculate total spectrum energy
            scalar totalEnergy = 0.0;
            forAll(F_array_,i)
            {
                if(SpectrumType_ == "Jonswap")
                {
                    totalEnergy += deltaF*Jonswap_(waveHeight_sig_, peakPeriod_, gamma_, F_array_[i]);
                }
                else if(SpectrumType_ == "PM")
                {
                    totalEnergy += deltaF*PM_(1/peakPeriod_, F_array_[i]);
                }
                else
                {
                        FatalErrorInFunction
                        << "Only Jonswap or PM spectrum is valid for random waves. "<< nl
                        << exit(FatalError);
                }
            }
            // calculate averaged energy
            averagedEnergy = totalEnergy/(NComponents_-1);

            // recalculate freqency array according to averaged energy
            // the start and end frequency is fixed during calculation
            for(int i = 1; i < NComponents_ - 1; i++)
            {
                    scalar currentEnergy_ = 0.0;
                    F_array_[i] = F_array_[i-1];
                    while(currentEnergy_ < averagedEnergy)
                    {
                            F_array_[i] += Vsmall_*deltaF;
                             if(SpectrumType_ == "Jonswap")
                             {
                                 currentEnergy_ = 0.5*(F_array_[i] - F_array_[i-1])
                                                 *(Jonswap_(waveHeight_sig_, peakPeriod_, gamma_, F_array_[i-1]) + Jonswap_(waveHeight_sig_, peakPeriod_, gamma_, F_array_[i]));
                             }
                             else if(SpectrumType_ == "PM")
                             {
                                 currentEnergy_ = 0.5*(F_array_[i] - F_array_[i-1])*(PM_(1/peakPeriod_, F_array_[i-1]) + PM_(1/peakPeriod_, F_array_[i]));
                             }
                    }
            }
            forAll(F_array_, i) // calculate represented frequency
            {
                    if(i > 0)  // from the second item of array
                    {
                            F_array_[i-1] = 0.5*(F_array_[i] + F_array_[i-1]);
                    }
            }
    }
    else if(AmplitudeDiscretization_ == "AF") // averaged frequency
    {}
    else
    {
            FatalErrorInFunction
            << "Only AE(averaged energy), "<< nl
            << "     AF(averaged frequency) are valid for random wave. "
            << exit(FatalError);
    }

    // update spectrum value
    if(SpectrumType_ == "Jonswap")
    {
        forAll(WS_, i)
        {
            WS_[i] = Jonswap_(waveHeight_sig_, peakPeriod_, gamma_, F_array_[i]);
        }
    }
    else if(SpectrumType_ == "PM")
    {
        forAll(WS_, i)
        {
            WS_[i] = PM_(1/peakPeriod_, F_array_[i]);
        }
    }

    // calculate amplitude
    Cal_A_array(deltaF, averagedEnergy);
}

void Foam::waterWaveModels::RandomWave::Cal_LkC_array()
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

void Foam::waterWaveModels::RandomWave::Cal_A_array
(
    scalar& deltaFreqency,
    scalar& deltaEnergy
)
{
    forAll(A_array_, i)
    {
        if(AmplitudeDiscretization_ == "AF")
        {
             A_array_[i] = sqrt(2.0*WS_[i]*deltaFreqency);
        }
        else
        {
             A_array_[i] = sqrt(2.0*deltaEnergy);
        }
    }
}

void Foam::waterWaveModels::RandomWave::Cal_Tf_array()
{
    forAll(Tf_array_,i)
    {
        Tf_array_[i] = 4.0*sinh(k_array_[i]*waterDepth_)*sinh(k_array_[i]*waterDepth_)/(2.0*k_array_[i]*waterDepth_+sinh(2.0*k_array_[i]*waterDepth_));
    }
}

void Foam::waterWaveModels::RandomWave::Cal_randomPhase()
{
    Random Randomephase_(1);
    forAll(RandomPhase_,i)
    {
        RandomPhase_[i] = 2.0*pi_*Randomephase_.sample01<scalar>();
    }
}

void Foam::waterWaveModels::RandomWave::Cal_EF_array()
{
    // independant variable for curve fitting formular
    List<scalar> theta(NComponents_, 0.0);

    // enhance the amplitude
    forAll(Ef_array_,i)
    {
        scalar relativeWaterDepth = waterDepth_/L_array_[i];
        scalar nonlinearity = 2*A_array_[i]/waterDepth_;

        // D/L value corresponds to the maximum theta
        scalar DL_Limit = 0.5/(1 - exp(4*nonlinearity) * nonlinearity);

        if(relativeWaterDepth > DL_Limit)// deep water waves, theta decreases with the increase of D/L
        {
            // deep water approximation: theta is calculated by interpolation
            theta[i] = theta[i-1] + 0.2*(mag(theta[i-1]-theta[i-2]) + mag(theta[i-2]-theta[i-3]) + mag(theta[i-3]-theta[i-4]) + mag(theta[i-4]-theta[i-5])+mag(theta[i-5]-theta[i-6]));
        }
        else
        {
            // shallow- and intermediate-water wave components
            theta[i] = (1 - (1 - exp(4*nonlinearity) * nonlinearity)*relativeWaterDepth)*relativeWaterDepth;
        }

        // enhanced factor
        Ef_array_[i] = 183.88*pow3(theta[i]) - 63.26*pow(theta[i],2) + 7.58*theta[i] + 0.72;

        // the minimum value is 1.0, maximum value is 3.115
        //Ef_array_[i] = min(max(1.0, Ef_array_[i]), 3.115);

        // minimum value is 1.0, no maximum value
        Ef_array_[i] = max(1.0, Ef_array_[i]);

        A_EF_array_[i] = A_array_[i];
        A_EF_array_[i] *= Ef_array_[i];
            /*if(AmplitudeDiscretization_ == "AF") // seperately enhanced
            {
                    // store and enhance
                    A_EF_array_[i] = A_array_[i];
                    A_EF_array_[i] *= Ef_array_[i];
            }*/
    }

	/*if(AmplitudeDiscretization_ == "AE") //averaged enhanced because all the components share indentical amplitude
	{
		scalar averagedEnhancedFactor = 0.0;
		forAll(A_array_, i)
		{
			// store
			A_EF_array_[i] = A_array_[i];
			// calculate total energy
			averagedEnhancedFactor += Ef_array_[i];
		}
		// average
		averagedEnhancedFactor /= (NComponents_-1);
		forAll(A_EF_array_, i)
		{
			// store
			Ef_array_[i] = averagedEnhancedFactor;
			// enhanced
			A_EF_array_[i] *= averagedEnhancedFactor;
		}
	}*/

}

Foam::word Foam::waterWaveModels::RandomWave::name()
{
	return "Random wave";
}

// For a irregular wave, the phase velocities have been
// embedded in the calculation of surface elevation
Foam::scalar Foam::waterWaveModels::RandomWave::Cal_C()
{
	return 1.0;
}

Foam::scalar Foam::waterWaveModels::RandomWave::eta(scalar t)
{
    scalar eta_ = 0.0;
    // array stores first order surface elevation (contain phase velocity)
    List<scalar> FirstOrderEta(NComponents_, 0.0);
    // array stores second order surface elevation (contain phase velocity)
    //List<scalar> SecondOrderEta(NComponents_, 0.0);

    forAll(FirstOrderEta,i)
    {
        FirstOrderEta[i] = A_EF_array_[i]*(cos(-2*pi_*F_array_[i]*t + RandomPhase_[i]));
            // second order for each components
            //SecondOrderEta[i] = 0.5*pow(A_EF_array_[i],2)*HijPos(i,i)*cos(2*(-2*pi_*F_array_[i]*t + RandomPhase_[i]));
            // second order wave-wave interactions
            /*for(int j = i+1; j < NComponents_; j++)
            {
                    scalar phase_sum = RandomPhase_[i] + RandomPhase_[j];
                    scalar phase_sub = RandomPhase_[i] - RandomPhase_[j];
            scalar omegaij_sum = 2*pi_*F_array_[i] + 2*pi_*F_array_[j];
            scalar omegaij_sub = 2*pi_*F_array_[i] - 2*pi_*F_array_[j];
            scalar phi_sum = -omegaij_sum *t + phase_sum;
            scalar phi_sub = -omegaij_sub *t - phase_sub;
            SecondOrderEta[i] += A_EF_array_[i]*A_EF_array_[j]*(HijPos(i, j)*cos(phi_sum) + HijNeg(i, j)*cos(phi_sub));
            }*/
            //eta_ += C_array_[i]*(FirstOrderEta[i] + SecondOrderEta[i]);
        eta_ += C_array_[i]*FirstOrderEta[i];
    }
    // slowly start factor
    scalar SSF_ = (t < peakPeriod_) ? pow(t/peakPeriod_,2.0) : 1.0;
    return eta_*SSF_;
}

// for piston type wavemaker
Foam::scalar Foam::waterWaveModels::RandomWave::displacement(scalar t, scalar deltaT, scalar WaterLevel)
{
    scalar correctWaterLevel = WaterLevel - waterDepth_;
    scalar PistonVelocity_ = 0.0;
    if(activeAbsorption_)
    {
        forAll(F_array_,i)
        {
            PistonVelocity_ += 2*pi_*F_array_[i]/Tf_array_[i]*(2*A_array_[i]*cos(-2*pi_*F_array_[i]*t + RandomPhase_[i]) - correctWaterLevel/(NComponents_-1));
        }
    }
    else
    {
        forAll(F_array_,i)
        {
            PistonVelocity_ += 2*pi_*F_array_[i]/Tf_array_[i]*A_array_[i]*cos(-2*pi_*F_array_[i]*t + RandomPhase_[i]);
        }
    }

    position_ += PistonVelocity_*deltaT;

    // second order correction
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
                secondOrderCorrection += (CombinedWave.X2_super(RandomPhase_[i]+RandomPhase_[j], 0, 0, t)*curr_TF_super).Re()+(CombinedWave.X2_sub(RandomPhase_[i]-RandomPhase_[j], 0, 0, t)*curr_TF_sub).Re();
            }
        }
    }
    return position_ + secondOrderCorrection;
}

void Foam::waterWaveModels::RandomWave::PrintWaveProperties()
{
    Info << "------------------------------" << "\n"
         << "Wave properties of random wave" << "\n"
         << " Wave spectrum = " << SpectrumType_ << "\n"
         << " Water Depth (m) = " << waterDepth_ << "\n"
         << " Significant Wave Height (m) = " << waveHeight_sig_ << "\n"
         << " Peak Wave Period (s) = " << peakPeriod_ << "\n"
         << "------------------------------" << endl;

    if(activeAbsorption_)
    {
        Info << "    Active absorption will be used after the third peak wave period" << endl;
    }
    else
    {
        Info << "    Active absorption is not adopted" << endl;
    }
    Info << "    The frequency, amplitude, wave number, random phase, enhanced factor and spectrum value of each component are: " << endl;
    forAll(F_array_, i)
    {
        Info << "    " << F_array_[i] << ", " <<  A_array_[i] << ", " << k_array_[i] << ", " << RandomPhase_[i] << ", " << Ef_array_[i] << ", " << WS_[i] <<  endl;
    }
}

void Foam::waterWaveModels::RandomWave::storeTFs()
{
    Info << "Calculating all the transfer functions for the random wave..." << endl;
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
