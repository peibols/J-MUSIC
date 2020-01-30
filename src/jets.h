#pragma once
#ifndef SRC_JETS_H
#define SRC_JETS_H_

#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "data.h"
#include "grid.h"
#include "Parton.h"
#include "hydro_source.h"
#include "eos.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "pretty_ostream.h"

using namespace Pythia8;

struct collision {
    double x,y,z,t;
    double px,py,pz,en;
};

class Jets {

 private:
    const InitData &DATA;
    vector<Parton> parton_list;
    vector<Parton> final_parton_list;
    vector<Parton> hadron_list;
    pretty_ostream music_message;

    Pythia pythia;	//Parton Shower
    Pythia hpythia;	//Hadroniztion
    Rndm randi;
 
    std::vector<std::shared_ptr<jet>> binary_list;

    vector<double> injected_momentum;

    double event_weight;
    double event_cross;

 public:
    Jets(const InitData &DATA_in);
    void GetBinaries();
    void InitJets(hydro_source &hydro_source_terms);
    void InitTestJets();
    void TreeInit();
    bool TreeDoer(vector<double> pos_coll, double total_cross, int icoll);
    bool EvolveJets(double tau, SCGrid &arena_current, const EOS &eos, hydro_source &hydro_source_terms, bool after_hydro_flag);
    void DoEloss(Parton &parton, double tau, SCGrid &arena_current, const EOS &eos, hydro_source &hydro_source_terms);
    void FinalPartons();
    void InitLund();
    void HadronizeJets();
    double GetFluidEnergy(double x, double y, double rap, SCGrid &arena_current);
    void GetFluidFlow(double x, double y, double rap, SCGrid &arena_current, double  *v_flow);
    void GetFluidGrid(int nx, double x, int &ix, double &dx, int ny, double y, int &iy, double &dy, int neta, double rap, int &ieta, double &deta);

};

#endif // SRC_JETS_H


