/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "LoadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

template <class ReactionThermo, class ThermoType>
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    LoadBalancedChemistryModel(const ReactionThermo& thermo)
    : StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
      cpuTimes_(this->mesh().cells().size(), 0.0), balancer_(LoadBalancer()),
      mapper_(createMapper(this->thermo()))
{

    Info << "Running with a load balanced" << endl;
    cpuSolveFile_ = logFile("cpu_solve.out");
    cpuSolveFile_() << "time"
                    << "    "
                    << "get_problem"
                    << "    "
                    << "updateState"
                    << "    "
                    << "balance"
                    << "    "
                    << "solveBuffer"
                    << "    "
                    << "unbalance"
                    << "    "
                    << "rank ID" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template <class ReactionThermo, class ThermoType>
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    ~LoadBalancedChemistryModel()
{
}

template <class ReactionThermo, class ThermoType>
mixtureFractionRefMapper
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::createMapper(
    const ReactionThermo& thermo)
{

    const IOdictionary chemistryDict_tmp(IOobject(
        thermo.phasePropertyName("chemistryProperties"),
        thermo.db().time().constant(),
        thermo.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false));
    return mixtureFractionRefMapper(chemistryDict_tmp, thermo.composition());
    // return balancer_ptr(new bulutLoadBalancing());
}

template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
scalar LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve(
    const DeltaTType& deltaT)
{

    using problem_buffer_t  = LoadBalancerBase::buffer_t<ChemistryProblem>;
    using solution_buffer_t = LoadBalancerBase::buffer_t<ChemistrySolution>;

    BasicChemistryModel<ReactionThermo>::correct();

    if(!this->chemistry_)
    {
        return great;
    }

    clockTime timer;

    timer.timeIncrement();
    DynamicList<ChemistryProblem> all_problems = getProblems(deltaT);
    scalar                        get_problem  = timer.timeIncrement();

    timer.timeIncrement();
    balancer_.updateState(all_problems);
    scalar updateState = timer.timeIncrement();

    balancer_.printState();

    timer.timeIncrement();
    problem_buffer_t balanced_problems = balancer_.balance(all_problems);
    scalar           balance           = timer.timeIncrement();

    timer.timeIncrement();
    solution_buffer_t balanced_solutions = solveBuffer(balanced_problems);
    scalar            solveBuffer        = timer.timeIncrement();

    timer.timeIncrement();
    solution_buffer_t my_solutions = balancer_.unbalance(balanced_solutions);
    scalar            unbalance    = timer.timeIncrement();

    cpuSolveFile_() << this->time().timeOutputValue() << "    " << get_problem
                    << "    " << updateState << "    " << balance << "    "
                    << solveBuffer << "    " << unbalance << "    "
                    << Pstream::myProcNo() << endl;

    Pstream::waitRequests();

    return updateReactionRates(my_solutions);
}

template <class ReactionThermo, class ThermoType>
void LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveSingle(
    ChemistryProblem& problem, ChemistrySolution& solution) const
{

    // TODO: Make this work so that ChemistryProblem is a const &

    scalar timeLeft = problem.deltaT;
    // const scalarList c0       = prob.c;
    scalarField c0 = problem.c;

    // Timer begin
    clockTime time;
    time.timeIncrement();

    const label arbitrary = 0;
    // Calculate the chemical source terms
    while(timeLeft > small)
    {
        scalar dt = timeLeft;
        // this->solve(prob.c, prob.Ti, prob.pi, dt, prob.deltaTChem);
        this->solve(
            problem.pi,
            problem.Ti,
            problem.c,
            arbitrary,
            dt,
            problem.deltaTChem);
        timeLeft -= dt;
    }

    solution.c_increment = (problem.c - c0) / problem.deltaT;
    solution.deltaTChem  = min(problem.deltaTChem, this->deltaTChemMax_);
    // Timer end
    solution.cpuTime = time.timeIncrement();
    solution.cellid  = problem.cellid;
    solution.rhoi    = problem.rhoi;
}

template <class ReactionThermo, class ThermoType>
scalar
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRates(
    const LoadBalancerBase::buffer_t<ChemistrySolution>& solutions)
{

    scalar deltaTMin = great;

    for(const auto& array : solutions)
    {
        for(const auto& solution : array)
        {

            for(label j = 0; j < this->nSpecie_; j++)
            {
                this->RR_[j][solution.cellid] =
                    computeReactionRate(j, solution);
            }

            deltaTMin = min(solution.deltaTChem, deltaTMin);
            this->deltaTChem_[solution.cellid] =
                min(solution.deltaTChem, this->deltaTChemMax_);
            this->cpuTimes_[solution.cellid] = solution.cpuTime;
        }
    }

    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
scalar LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve(
    const scalarField& deltaT)
{
    return this->solve<scalarField>(deltaT);
}

template <class ReactionThermo, class ThermoType>
scalar LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve(
    const scalar deltaT)
{
    // Don't allow the time-step to change more than a factor of 2
    return min(
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT);
}

template <class ReactionThermo, class ThermoType>
LoadBalancerBase::buffer_t<ChemistrySolution>
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveBuffer(
    LoadBalancerBase::buffer_t<ChemistryProblem>& problems) const
{

    // allocate the solutions buffer
    LoadBalancerBase::buffer_t<ChemistrySolution> solutions;
    for(label i = 0; i < problems.size(); ++i)
    {
        DynamicList<ChemistrySolution> sublist(
            problems[i].size(), ChemistrySolution(this->nSpecie_));
        solutions.append(sublist);
    }

    for(label i = 0; i < solutions.size(); ++i)
    {
        for(label j = 0; j < solutions[i].size(); ++j)
        {

            solveSingle(problems[i][j], solutions[i][j]);
        }
    }

    return solutions;
}

template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
DynamicList<ChemistryProblem>
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getProblems(
    const DeltaTType& deltaT)
{

    // TODO: Add refcell and Treact as conditions to get problems.

    const scalarField&            T = this->thermo().T();
    const scalarField&            p = this->thermo().p();
    tmp<volScalarField>           trho(this->thermo().rho());
    const scalarField&            rho          = trho();
    bool                          refCellFound = false;
    DynamicList<ChemistryProblem> chem_problems;
    ChemistrySolution             ref_solution(this->nSpecie_);

    // TODO: reserve
    // chem_problems.reserve(T.size());

    forAll(p, celli)
    {

        for(label i = 0; i < this->nSpecie_; i++)
        {
            this->c_[i] = computeConcentration(rho[celli], i, celli);
        }

        scalar Ti = T[celli];

        if(Ti > this->Treact())
        {

            // Create the problem
            ChemistryProblem problem;
            problem.c          = this->c_;
            problem.Ti         = T[celli];
            problem.pi         = p[celli];
            problem.rhoi       = rho[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.deltaT     = deltaT[celli];
            problem.cpuTime    = this->cpuTimes_[celli];
            problem.cellid     = celli;

            // First reference cell is found and solved, following reference cells 
            // are mapped from the first reference cell found
            if(mapper_.active() && mapper_.shouldMap(getMassFraction(problem)))
            {

                if(!refCellFound)
                {
                    solveSingle(problem, ref_solution);
                    updateReactionRate(ref_solution, ref_solution.cellid);
                    refCellFound                         = true;
                    this->cpuTimes_[ref_solution.cellid] = ref_solution.cpuTime;
                }

                else
                {
                    updateReactionRate(ref_solution, celli);
                    this->cpuTimes_[celli] = ref_solution.cpuTime;
                }
            }
            else
            {
                chem_problems.append(problem);
            }
        }
        else
        {
            for(label i = 0; i < this->nSpecie(); i++)
            {
                this->RR_[i][celli] = 0;
            }
        }
    }
    return chem_problems;
}

template <class ReactionThermo, class ThermoType>
void LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRate(
    const ChemistrySolution& solution, const label& i)
{

    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][i] = computeReactionRate(j, solution);
    }
    this->deltaTChem_[i] = min(solution.deltaTChem, this->deltaTChemMax_);
}

template <class ReactionThermo, class ThermoType>
scalar
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::computeConcentration(
    const scalar& rho, const label& i, const label& celli) const
{

    return (rho * this->Y_[i][celli] / this->specieThermos_[i].W());
}

template <class ReactionThermo, class ThermoType>
scalar
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::computeReactionRate(
    const label& j, const ChemistrySolution& solution) const
{

    return (solution.c_increment[j] * this->specieThermos_[j].W());
}

template <class ReactionThermo, class ThermoType>
scalarField
LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getMassFraction(
    const ChemistryProblem& problem) const
{

    scalarField tmp(this->nSpecie_);
    for(label i = 0; i < this->nSpecie_; i++)
    {
        tmp[i] = this->Y_[i][problem.cellid];
    }

    return (tmp);
}

} // namespace Foam