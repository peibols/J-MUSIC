#include "jets.h"

using namespace std;

void Jets::GetFluidGrid(int nx, double x, int &ix, double &dx, int ny, double y, int &iy, double &dy, int neta, double rap, int &ieta, double &deta) {
    
    double delta_x=DATA.delta_x;
    double delta_y=DATA.delta_y;
    double delta_rap=DATA.delta_eta;
    
    ix=int(x/delta_x)+nx/2;
    dx=(x-double(ix-nx/2)*delta_x)/delta_x;

    iy=int(y/delta_y)+ny/2;
    dy=(y-double(iy-ny/2)*delta_y)/delta_y;
     
    ieta=int(rap/delta_rap)+neta/2;
    deta=(rap-double(ieta-neta/2)*delta_rap)/delta_rap;
    
}
double Jets::GetFluidEnergy(double x, double y, double rap, SCGrid &arena_current) {
    	
    int ix, iy, ieta;
    double dx, dy, deta;

    int nx=arena_current.nX();
    int ny=arena_current.nY();
    int neta=arena_current.nEta();
    
    GetFluidGrid(nx, x, ix, dx, ny, y, iy, dy, neta, rap, ieta, deta);

    double energy=0.;

    if (ix < 0 || ix >= nx-1 || iy < 0 || iy >= ny-1 || ieta < 0 || ieta >= neta-1) return energy;
    
    energy += arena_current  (ix      , iy      , ieta      ).epsilon * (1.-dx) * (1.-dy) * (1.-deta);
    energy += arena_current  (ix      , iy + 1  , ieta      ).epsilon * (1.-dx) * dy      * (1.-deta);
    energy += arena_current  (ix + 1  , iy      , ieta      ).epsilon * dx      * (1.-dy) * (1.-deta);
    energy += arena_current  (ix + 1  , iy + 1  , ieta      ).epsilon * dx      * dy      * (1.-deta);
    energy += arena_current  (ix      , iy      , ieta + 1  ).epsilon * (1.-dx) * (1.-dy) * deta     ;
    energy += arena_current  (ix      , iy + 1  , ieta + 1  ).epsilon * (1.-dx) * dy      * deta     ;
    energy += arena_current  (ix + 1  , iy      , ieta + 1  ).epsilon * dx      * (1.-dy) * deta     ;
    energy += arena_current  (ix + 1  , iy + 1  , ieta + 1  ).epsilon * dx      * dy      * deta     ;

    return energy;
}
void Jets::GetFluidFlow(double x, double y, double rap, SCGrid &arena_current, double *v_flow) {
    	
    int ix, iy, ieta;
    double dx, dy, deta;
    
    int nx=arena_current.nX();
    int ny=arena_current.nY();
    int neta=arena_current.nEta();
    
    GetFluidGrid(nx, x, ix, dx, ny, y, iy, dy, neta, rap, ieta, deta);

    for (unsigned a=0; a<4; a++) v_flow[a]=0.;

    if (ix < 0 || ix >= nx-1 || iy < 0 || iy >= ny-1 || ieta < 0 || ieta >= neta-1) return;
    
    for (int com=1; com<=3; com++) {
        v_flow[com-1] += arena_current  (ix      , iy      , ieta      ).u[com] * (1.-dx) * (1.-dy) * (1.-deta);
        v_flow[com-1] += arena_current  (ix      , iy + 1  , ieta      ).u[com] * (1.-dx) * dy      * (1.-deta);
        v_flow[com-1] += arena_current  (ix + 1  , iy      , ieta      ).u[com] * dx      * (1.-dy) * (1.-deta);
        v_flow[com-1] += arena_current  (ix + 1  , iy + 1  , ieta      ).u[com] * dx      * dy      * (1.-deta);
        v_flow[com-1] += arena_current  (ix      , iy      , ieta + 1  ).u[com] * (1.-dx) * (1.-dy) * deta     ;
        v_flow[com-1] += arena_current  (ix      , iy + 1  , ieta + 1  ).u[com] * (1.-dx) * dy      * deta     ;
        v_flow[com-1] += arena_current  (ix + 1  , iy      , ieta + 1  ).u[com] * dx      * (1.-dy) * deta     ;
        v_flow[com-1] += arena_current  (ix + 1  , iy + 1  , ieta + 1  ).u[com] * dx      * dy      * deta     ;
    }

    v_flow[3]=std::sqrt(1.+v_flow[0]*v_flow[0]+v_flow[1]*v_flow[1]+v_flow[2]*v_flow[2]);

}
void Jets::GetBinaries() {
        
    music_message.info(DATA.binsName);
    ifstream bins_file(DATA.binsName.c_str());
        
    string dummy;
    getline(bins_file,dummy);

    binary_list.clear();

    while (!bins_file.eof()) {
        std::shared_ptr<jet> new_collision(new jet);
        int dum;
	bins_file >> dum >> new_collision->x_perp >> new_collision->y_perp;
        if (bins_file.eof()) break;
	new_collision->eta_source=0.;
	new_collision->tau_form=0.;
	new_collision->dpxdtau=0.;
	new_collision->dpydtau=0.;
	new_collision->dpzdtau=0.;
	new_collision->dEdtau=0.;

	binary_list.push_back(new_collision);
    }
    music_message << "# Binary Collisions=" << binary_list.size();		
    music_message.flush("info");
    music_message << "Using Initial_profile=" << DATA.binsName;
    music_message.flush("info");
}

