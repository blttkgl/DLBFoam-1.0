/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM-Aalto library, derived from OpenFOAM.

    https://github.com/blttkgl/OpenFOAM-Aalto

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

#include "mixtureFractionRefMapper.H"

bool Foam::mixtureFractionRefMapper::check_if_refcell(const ChemistryProblem& problem)
{
    // TODO: These two lines should be called ONCE 
    auto beta_of = mixture_fraction_.get_beta();
    auto alpha   = mixture_fraction_.get_alpha();

    scalar beta = 0.0;
    forAll(problem.c, iField)
    {
        beta += alpha[iField] * problem.c[iField];
    }

    scalar Z = (beta - beta_of[0]) / (beta_of[1] - beta_of[0]);


    if(refCellFound_ && Z < tolerance_ && abs(problem.Ti-globalReference.Ti) < temperature_tolerance_ )
    {
        return true;
    }
    else if(!refCellFound_ && Z < tolerance_)
    {
        return true;
    }
    else
    {
        return false;
    }
    

}

bool Foam::mixtureFractionRefMapper::shouldMap(const ChemistryProblem& problem)
{
    return check_if_refcell(problem);
}


Foam::ChemistryProblem Foam::mixtureFractionRefMapper::getGlobalRef(const ChemistryProblem& problem)
{
    auto refMapList = LoadBalancerBase::allGather(problem);

    // Return the cell with the maximum temperature as the globalRef
    auto filter = [](const ChemistryProblem& lhs, const ChemistryProblem& rhs){
        return lhs.Ti < rhs.Ti;};

    globalReference = *std::max_element(refMapList.begin(), refMapList.end(), filter);
    refCellFound_ = true;

    return(globalReference);
}    


