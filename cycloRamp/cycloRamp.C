/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "cycloRamp.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

#include "transformField.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(cycloRamp, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        cycloRamp,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::cycloRamp::cycloRamp
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::cycloRamp::~cycloRamp()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::cycloRamp::transformation() const
{
    scalar t = time_.value();
    //Calculate azimuth (rotation of outer rotor zone) with ramp.
    vector azimuth(0,0,0);
    vector basePhi=t*OMEGA0_;
    vector dPhiDot=OMEGAF_-OMEGA0_;
    azimuth+=basePhi;
    if(t>TF_){
	azimuth+=(TF_-T0_)*dPhiDot/2+(t-TF_)*OMEGAF_;
    } else if (t>T0_) {
	vector phiDotDiff=(t-T0_)/(TF_-T0_)*dPhiDot;
	azimuth+=(t-T0_)*phiDotDiff/2;
    }
    scalar magAzimuth = azimuth.z();
    //End ramp azimuth calculation

    scalar amplitude = degToRad(max_-min_)/2;
    //scalar average = degToRad(max_+min_)/2;
    scalar pitch = amplitude*cos(magAzimuth+degToRad(phi_-tilt_))-amplitude*cos(degToRad(phi_-tilt_));

    //This assumes that the rotation vector is exclusively in the z-direction.
    point bladeCofr;
    bladeCofr.x()=COFR_.x()+radius_*cos(magAzimuth+degToRad(phi_+90));
    bladeCofr.y()=COFR_.y()+radius_*sin(magAzimuth+degToRad(phi_+90));
    bladeCofr.z()=COFR_.z();

    quaternion R(0,0,azimuth.z());
    septernion TR(septernion(COFR_)*R*septernion(-COFR_));
    quaternion R2(0,0,-pitch);
    TR=septernion(bladeCofr)*R2*septernion(-bladeCofr)*TR;

    Info<< "solidBodyMotionFunctions::cycloRamp::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::cycloRamp::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("tilt") >> tilt_;
    SBMFCoeffs_.lookup("max") >> max_;
    SBMFCoeffs_.lookup("min") >> min_;
    SBMFCoeffs_.lookup("radius") >> radius_;
    SBMFCoeffs_.lookup("phi") >> phi_;
    SBMFCoeffs_.lookup("OMEGA0") >> OMEGA0_;
    SBMFCoeffs_.lookup("OMEGAF") >> OMEGAF_;
    SBMFCoeffs_.lookup("T0") >> T0_;
    SBMFCoeffs_.lookup("TF") >> TF_;
    SBMFCoeffs_.lookup("COFR") >> COFR_;

    return true;
}

// ************************************************************************* //
