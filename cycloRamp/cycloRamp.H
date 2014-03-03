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

Class
    Foam::solidBodyMotionFunctions::cycloRamp

Description
    Allows for ramped nested AMI movement; currently only allows for a rotation vector completely in z-direction.

SourceFiles
    cycloRamp.C

\*---------------------------------------------------------------------------*/

#ifndef cycloRamp_H
#define cycloRamp_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class cycloRamp Declaration
\*---------------------------------------------------------------------------*/

class cycloRamp
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point COFR_;
	vector OMEGA0_;
	vector OMEGAF_;
	scalar T0_;
	scalar TF_;
	scalar radius_;
	scalar tilt_;
	scalar max_;
	scalar min_;
	scalar phi_;

    // Private Member Functions

        //- Disallow copy construct
        cycloRamp(const cycloRamp&);

        //- Disallow default bitwise assignment
        void operator=(const cycloRamp&);


public:

    //- Runtime type information
    TypeName("cycloRamp");


    // Constructors

        //- Construct from components
        cycloRamp
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~cycloRamp();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
