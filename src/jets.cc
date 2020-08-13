#include "jets.h"
#include "vector"
#include "vector_operators.h" 
#include <stdlib.h>
#include <time.h>

using namespace std;

Jets::Jets(const InitData &DATA_in) : DATA(DATA_in) 
{
  for (unsigned a=0; a<4; a++) injected_momentum.push_back(0.);

  surface_in_binary = DATA.freeze_surface_in_binary; 

  boost_invariant = DATA.boost_invariant;

  bulk_deltaf_kind = 1;

  // Tolerance for medium hadronization
  x_tol = 2.*DATA.delta_x;
  y_tol = 2.*DATA.delta_y;
  tau_tol = 2.*DATA.delta_tau;
  eta_tol = 2.*DATA.delta_eta;

  hadfile.open("hadrons_list.dat");
  negafile.open("sampled_partons_list.dat");
}

void Jets::InitJets(hydro_source &hydro_source_terms) {

    srand (time(NULL));

    //total cross section[mb] 
    double total_cross=70.; //for 2.76 TeV

    //Clear parton list
    parton_list.clear();
    
    // Does PYTHIA tree (int nhat, int seed, vector<Parton> parton_list)
    TreeInit();

    //Loop over #Coll
    int njets=0;
    int Ncoll=binary_list.size();
   // Ncoll=100;
    if (DATA.one_hard_collision==1) {
       vector<double> coll_pos;
       int icoll;
       if (DATA.smooth_glauber==1) {
	icoll = 0;
	double x_perp, y_perp;
        get_smooth_xy_point(x_perp, y_perp); 
	vector<double> smooth_coll_pos {x_perp, y_perp, 0., 0.};
        coll_pos = smooth_coll_pos;
      }
      else {
        icoll = rand() % Ncoll;
        vector<double> point_coll_pos {binary_list[icoll]->x_perp,binary_list[icoll]->y_perp,binary_list[icoll]->eta_source,binary_list[icoll]->tau_form}; 
        coll_pos = point_coll_pos;
      }
      bool put_jet=0;
      while (!put_jet) put_jet=TreeDoer(coll_pos, total_cross, icoll);
      njets+=1;
    }
    else {
      for (int icoll=0; icoll<Ncoll; icoll++) {
          //cout << " icoll= " << icoll  << endl;
  	  vector<double> coll_pos {binary_list[icoll]->x_perp,binary_list[icoll]->y_perp,binary_list[icoll]->eta_source,binary_list[icoll]->tau_form}; 
          bool put_jet=TreeDoer(coll_pos, total_cross, icoll);
          if (put_jet) {
	      njets+=1;
              //hydro_source_terms.update_sources(binary_list[icoll]);
	  }
      }
    }
    cout << " Njets= " << njets << endl;
}
void Jets::InitTestJets() {

    parton_list.clear();

    srand (DATA.Jets_seed+2);

    if (DATA.smooth_glauber==0) {
      double test_px=8.768*cos(M_PI/4.);
      double test_py=8.768*sin(M_PI/4.);
      double test_pz=0.;
      double test_en=8.768;
      int test_id=1;
      //Smooth tests
      parton_list.push_back ( Parton ( test_px, test_py, test_pz, test_en, 0., 0., 0, -1, -1, test_id, "test", 0, 0, true ) );
      vector<double> test_pos {0.,0.,0.,0.};
      parton_list[parton_list.size()-1].vSetPos(test_pos);
      parton_list[parton_list.size()-1].vSetPosIn(test_pos);
    }
    else {
      double en=10.;
      double phi=(double)rand() / RAND_MAX * 2. * M_PI;
      ofstream jetphi("jetphi.dat");
      jetphi << phi << endl;
      jetphi.close();
      double px=en*cos(phi);
      double py=en*sin(phi);
      double pz=0.;
      int id=1;
      parton_list.push_back ( Parton ( px, py, pz, en, 0., 0., 0, -1, -1, id, "test", 0, 0, true ) );
      double x_perp, y_perp;
      get_smooth_xy_point(x_perp, y_perp);
      vector<double> smooth_coll_pos {x_perp, y_perp, 0., 0.};
      parton_list[parton_list.size()-1].vSetPos(smooth_coll_pos);
      parton_list[parton_list.size()-1].vSetPosIn(smooth_coll_pos);
    }
}

void Jets::FinalPartons() {
    final_parton_list.clear();
    int initial_partons=0;
    for (int iparton=0; iparton<parton_list.size(); iparton++) {
        Parton parton=parton_list[iparton];
	int d1=parton.GetD1();
	if (d1!=-1) continue;
	if (parton.GetOrig()=="rem") continue;

	initial_partons+=1;
	//Check if it has been dealt with. If not, means a mother was absorbed
	if (parton.vGetPos()[3]==-1.) {
            if (parton.vGetP()[3]==0.) {cout << " PROBLEM, not activated and still 0 energy" << endl;}
	    else {
	        int m1=parton.GetMom();
		if (m1==-1) cout << " NO MOM and still NOT DEALT WITH " << endl;
		bool found_the_mother=0;
		while (!found_the_mother) {
		    double ener_mom=parton_list[m1].vGetP()[3];
		    if (ener_mom==0.) {
			found_the_mother=1;
			parton.vSetP(parton.vGetP()*0.);
		    }
		    else if (m1==-1) {
                        cout << " COULD NOT FIND THE MOM !" << endl;
		    }
		    m1=parton_list[m1].GetMom();
		}
	    }
	}
	else {
	    if (parton.vGetP()[3]!=0.) {
	        final_parton_list.push_back(parton);
		//cout << " Lambda of this parton= " << parton.vGetP()[3]/parton.vGetPIn()[3] << endl;
		//
		// Medium Hadronization
                if (DATA.single_parton==0) {
		  SampleSurface(parton);
		  if (fin_and_therm_parton_list.size()>0) HadronizeTherm();
		}
	    }
	}
    }

    // Corona hadronization
    if (corona_parton_list.size()>0) HadronizeCorona();

    hadfile.close();
    negafile.close();

    //End of Parton Loop 
    cout << " Final partons list size= " << final_parton_list.size() << " and before quenching= " << initial_partons << endl;

    int nxmax=400;
    double hstep_x=0.05;
    int nymax=400;
    double hstep_y=0.05;
    double creat_hist[nxmax][nymax]={0.,0.};
    double creat_n[nxmax][nymax]={0.,0.};

    ofstream partons_file("final_partons_file.dat");
    for (int iparton=0; iparton<final_parton_list.size(); iparton++) {
        vector<double> pfin=final_parton_list[iparton].vGetP();
	vector<double> xfin=final_parton_list[iparton].vGetPos();
	partons_file << xfin[0] << " " << xfin[1] << " " << xfin[2] << " " << xfin[3] << " "
		     << pfin[0] << " " << pfin[1] << " " << pfin[2] << " " << pfin[3] << " "
		     << final_parton_list[iparton].GetId() << endl;
        //Fill creation position hist
	vector<double> x_in=final_parton_list[iparton].vGetPosIn();
	int nx,ny;
	if (x_in[0]!=x_in[0]) cout << " proooblem in hist ";
	if (x_in[0]>0.) nx=int(x_in[0]/hstep_x)+nxmax/2;
	else nx=int(x_in[0]/hstep_x)-1+nxmax/2;
	if (x_in[1]>0.) ny=int(x_in[1]/hstep_y)+nymax/2;
	else ny=int(x_in[1]/hstep_y)-1+nymax/2;
//	double heta=1./2.*log(()/());
	//int neta=x_in[0]/hstep_x;
	if (nx>nxmax-1) nx=nxmax-1;
	if (nx<0) nx=0;
	if (ny>nymax-1) ny=nymax-1;
	if (ny<0) ny=0;
	//cout << " x= " << x_in[0] << " y= " << x_in[1] << endl;
	if (x_in[0]<-10.) cout << " nx==" << nx << endl;
	//cout << " Lambda when hist= " << scientific << setprecision(8) << final_parton_list[iparton].vGetP()[3]/final_parton_list[iparton].vGetPIn()[3] << endl;
	//final_parton_list[iparton].display();
	creat_hist[nx][ny]+=final_parton_list[iparton].vGetP()[3]/final_parton_list[iparton].vGetPIn()[3];
	creat_n[nx][ny]+=1.;
    }

    ofstream creat_file("creation_final_energy_weighted_density.dat");
    for (unsigned a=0; a<nxmax; a++) {
        for (unsigned b=0; b<nymax; b++) {
	    creat_file << hstep_x*double(a)-nxmax/2.*hstep_x << " " << hstep_y*double(b)-nymax/2.*hstep_y << " " << creat_hist[a][b]/max(creat_n[a][b],1.) << endl;
	}
	creat_file << endl;
    }
    partons_file.close();
}
