#pragma once
#include <stdlib.h>     //for using the function sleep
namespace Foam {


void sanity_check(volScalarField& p, volScalarField& rho, PtrList<volScalarField>& Y, psiReactionThermo& thermo) {

    basicChemistryModel* model1 = new ode<StandardChemistryModel<psiReactionThermo, ThermoType>>(thermo);
    basicChemistryModel* model2 = new ode<LoadBalancedChemistryModel<psiReactionThermo, ThermoType>>(thermo);
        


        
    //basicChemistryModel* model2 = new ode<StandardChemistryModel<psiReactionThermo, ThermoType>>(thermo);

    

    for (size_t i = 0; i < 2; ++i){
    set_master_heavy(rho, thermo);
        auto asd = model1->solve(detail::FLOWTIMESTEP);
    set_master_heavy(rho, thermo);
        auto asd2 = model2->solve(detail::FLOWTIMESTEP);
    }



    for (label i = 0; i < model1->nSpecie(); ++i){
        
        
        const scalarField& r1 = model1->RR(i);
        const scalarField& r2 = model2->RR(i);

        scalar sum = 0.0;
        forAll(r1, celli) {

            double a = r1[celli];
            double b = r2[celli];
            double diff = std::abs(a-b);
            if (diff > 1E-9){
                Pout << "Load balanced model and standard model give different answers! Difference: " << diff << endl; 
                
                std::exit(1);
            }
            
        

        }


    }
    

    delete model1; delete model2;

}

}