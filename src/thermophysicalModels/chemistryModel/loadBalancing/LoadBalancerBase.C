/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "LoadBalancerBase.H"

namespace Foam
{

bool LoadBalancerBase::active() const
{
    return true;
}

void LoadBalancerBase::setState(const BalancerState& state)
{

    state_ = state;

    if(!validState())
    {

        printState();

        throw error("Invalid balance state in ChemistryLoadBalancing.");
    }
}

ChemistryLoad LoadBalancerBase::computeLoad(
    const DynamicList<ChemistryProblem>& problems)
{
    auto lambda = [](scalar sum, const ChemistryProblem& rhs) {
        return sum + rhs.cpuTime;
    };
    scalar sum =
        std::accumulate(problems.begin(), problems.end(), scalar(0), lambda);
    return ChemistryLoad(Pstream::myProcNo(), sum);
}

scalar
LoadBalancerBase::getMean(const DynamicList<ChemistryLoad>& loads)
{

    auto op = [](scalar sum, const ChemistryLoad& load) {
        return sum + load.value;
    };
    return std::accumulate(loads.begin(), loads.end(), 0.0, op) / loads.size();
}
ChemistryLoad
LoadBalancerBase::getMin(const DynamicList<ChemistryLoad>& vec)
{

    auto comp = [](const ChemistryLoad& lhs, const ChemistryLoad& rhs) {
        return lhs.value < rhs.value;
    };
    return *std::min_element(vec.begin(), vec.end(), comp);
}

ChemistryLoad
LoadBalancerBase::getMax(const DynamicList<ChemistryLoad>& vec)
{
    auto comp = [](const ChemistryLoad& lhs, const ChemistryLoad& rhs) {
        return lhs.value < rhs.value;
    };
    return *std::max_element(vec.begin(), vec.end(), comp);
}

label LoadBalancerBase::rankToLoadIdx(
    const DynamicList<ChemistryLoad>& loads, label rank) const
{

    for(label i = 0; i < loads.size(); ++i)
    {
        if(loads[i].rank == rank)
            return i;
    }

    throw "Could not find my rank from loads.";
}

bool LoadBalancerBase::validState() const
{

    /*
    auto sources      = state_.sources;
    auto destinations = state_.destinations;

    if (sources.size() < 1) return false;
    if (destinations.size() < 1) return false;
    if (sources.size() > size_t(Pstream::nProcs())) return false;
    if (destinations.size() > size_t(Pstream::nProcs())) return false;
    if (sources[0] != Pstream::myProcNo())
        return false; // the first index must always be this process
    if (sources[0] != destinations[0]) return false;
    */
    // TODO: enable
    return true;
}

void LoadBalancerBase::printState() const
{

    // receiver
    if(state_.sources.size() > 1)
    {
        Pout << "Receiver rank: " << Pstream::myProcNo() << " receives from "
             << vectorToString(state_.sources) << " own problems "
             << vectorToString(state_.nProblems) << endl;
    }
    // sender
    else if(state_.destinations.size() > 1)
    {
        Pout << "Sender rank: " << Pstream::myProcNo()
             << " sends to: " << vectorToString(state_.destinations)
             << " counts: " << vectorToString(state_.nProblems)
             << " remaining problems:  " << state_.nProblems.back() << endl;
    }

    else
    {
        Pout << "Rank: " << Pstream::myProcNo()
             << " does not take part in balancing.";
    }
}

} // namespace Foam