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

#include "nestedRotationMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

#include "transformField.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(nestedRotationMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        nestedRotationMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::nestedRotationMotion::nestedRotationMotion
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

Foam::solidBodyMotionFunctions::nestedRotationMotion::~nestedRotationMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::nestedRotationMotion::transformation() const
{
    scalar t = time_.value();

    vector omegaHouse
    (
        t*degToRad(radialVelocityHouse_.x()),
        t*degToRad(radialVelocityHouse_.y()),
        t*degToRad(radialVelocityHouse_.z())
    );
    vector omegaNested
    (
        t*degToRad(radialVelocityNested_.x()),
        t*degToRad(radialVelocityNested_.y()),
        t*degToRad(radialVelocityNested_.z())
    );

    scalar magOmegaHouse = mag(omegaHouse);
    scalar magOmegaNested = mag(omegaNested);

    quaternion R(omegaHouse/magOmegaHouse, magOmegaHouse);//construct a rotation quaternion from direction and magnitude.

    //This assumes that the rotation vector is exclusively in the z-direction.
    scalar nestedAngle=omegaHouse.z()+degToRad(nestedOmega0_);
    point CofGNested;
    CofGNested.x()=CofGHouse_.x()+nestedRadius_*cos(nestedAngle);
    CofGNested.y()=CofGHouse_.y()+nestedRadius_*sin(nestedAngle);
    CofGNested.z()=CofGHouse_.z();

    septernion TR(	septernion(CofGHouse_)*R*septernion(-CofGHouse_)	);
if(magOmegaNested>0){
    quaternion R2(omegaNested/magOmegaNested, magOmegaNested);
//    septernion TR(	septernion(CofGNested)*R2*septernion(-CofGNested)*
//			septernion(CofGHouse_)*R*septernion(-CofGHouse_)	);
    TR=septernion(CofGNested)*R2*septernion(-CofGNested)*TR;
}

    Info<< "solidBodyMotionFunctions::nestedRotationMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::nestedRotationMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("nestedOmega0") >> nestedOmega0_;
    SBMFCoeffs_.lookup("nestedRadius") >> nestedRadius_;
    SBMFCoeffs_.lookup("CofGHouse") >> CofGHouse_;
    SBMFCoeffs_.lookup("radialVelocityNested") >> radialVelocityNested_;
    SBMFCoeffs_.lookup("radialVelocityHouse") >> radialVelocityHouse_;

    return true;
}

// ************************************************************************* //
