#include "jets.h"

using namespace std;

void Jets::GetFluidGrid(int nx, double x, int &ix, double &dx, int ny, double y, int &iy, double &dy, int neta, double rap, int &ieta, double &deta) {
    
    double delta_x=DATA.delta_x;
    double delta_y=DATA.delta_y;
    double delta_rap=DATA.delta_eta;
    
    //std::cout << " JET delta_x = " << delta_x << " nx=  " << nx << std::endl;

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
   
    //cout << " energy sin mas= " << arena_current  (ix      , iy      , ieta      ).epsilon << endl;

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
    if (bins_file.fail()) {
      cout << " No Binaries file! " << endl;
      exit(1);
    }
        
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
void Jets::get_smooth_xy_point(double &x, double &y) {
 
  int g_maxx;
  int g_maxy;
  double g_deltax;
  double g_deltay;
 
  ifstream initial("./u_field.dat");
  if (initial.fail()) { cout << " No initial file in smooth glauber! " << endl; exit(1); }
 
  // Read profile 
  string s, sub;
  getline(initial,s);
  istringstream isi(s);
  string scrap;
  isi >> scrap >> scrap >> scrap >> scrap >> scrap >> scrap >> g_maxx >> scrap >> g_maxy >> scrap >> scrap >> scrap >> g_deltax >> scrap >> g_deltay;
  cout << " g_maxx= " << g_maxx << " g_maxy= " << g_maxy << " g_deltax= " << g_deltax << " g_deltay= " << g_deltay << endl;
  vector<vector<double>> glaub(g_maxx, vector<double>(g_maxy));

  double crap, xc, yc, edens;
  double maxedens=0.;
  int ix=0, iy=0;
  do {
    getline(initial,s);
    if (s.empty()) continue;
    istringstream iss(s);
    iss >> crap >> xc >> yc >> edens;
    //cout << " edens= " << edens << endl;
    if (edens>maxedens) maxedens=edens;
    glaub[ix][iy]=edens;
    iy+=1;
    if (iy==g_maxy) ix+=1, iy=0;
  } while(ix<g_maxx);

  // Normalise
  for (int a=0; a<200; a++) {
    for (int b=0; b<200; b++) {
      glaub[a][b]/=maxedens;
    }
  }

  // Get x and y
  inline double gGlaub(double x, double y,
                double g_deltax, double g_deltay,
                int g_maxx, int g_maxy,
                vector<vector<double>> glaub);


  srand(DATA.Jets_seed+2);
  while (true) {
    double rho=sqrt(150.*((double)rand()/RAND_MAX));
    double phi=2.*3.141592654*((double)rand()/RAND_MAX);
    x=rho*cos(phi);
    y=rho*sin(phi);
    double P = (double)rand()/RAND_MAX;
    if (gGlaub(x,y,g_deltax,g_deltay,g_maxx,g_maxy,glaub)>P) break;
  }

  // DEBUG
  //x=0.148;
  //y=-1.027;

  cout << "X= " << x << " Y= " << y << endl;
  return;

}
double gGlaub(double x, double y,
		double g_deltax, double g_deltay,
		int g_maxx, int g_maxy,
		vector<vector<double>> glaub)
{
  double gdens=0.;

  int ix, dx, iy, dy;

  if (x>=0.) {
    ix = int(x/g_deltax)+(g_maxx)/2;
    dx = (x - double(ix-(g_maxx)/2)*g_deltax)/g_deltax;
  }
  else {
    ix = int(x/g_deltax)+(g_maxx)/2-1;
    dx = (x - double(ix-(g_maxx)/2)*g_deltax)/g_deltax;
  }

  if (y>=0.) {
    iy = int(y/g_deltay)+(g_maxy)/2;
    dy = (y - double(iy-(g_maxy)/2)*g_deltay)/g_deltay;
  }
  else {
    iy = int(y/g_deltay)+(g_maxy)/2-1;
    dy = (y - double(iy-(g_maxy)/2)*g_deltay)/g_deltay;
  }

  if (ix<0 || ix>=g_maxx-1 || iy<0 || iy>=g_maxy-1) return gdens;
  gdens=glaub[ix][iy]*(1.-dx)*(1.-dy);
  gdens+=glaub[ix][iy+1]*(1.-dx)*dy;
  gdens+=glaub[ix+1][iy]*dx*(1.-dy);
  gdens+=glaub[ix+1][iy+1]*dx*dy;

  if (gdens>1.) cout << " gGlaub not properly normalised: gdens = " << gdens << endl;
  return gdens;
}














