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
  x_tol = 0.1;
  y_tol = 0.1;
  tau_tol = 0.1;
  eta_tol = 0.1;

  hadfile.open("hadrons_list.dat");
  negafile.open("sampled_partons_list.dat");
}

void Jets::InitJets(hydro_source &hydro_source_terms) {

    srand (time(NULL));

    //total cross section[mb]
    double total_cross=70.;

    //Clear parton list
    parton_list.clear();
    
    // Does PYTHIA tree (int nhat, int seed, vector<Parton> parton_list)
    TreeInit();

    //Loop over #Coll
    int njets=0;
    int Ncoll=binary_list.size();
   // Ncoll=100;
    if (DATA.one_hard_collision==1) {
      int icoll = rand() % Ncoll;
      vector<double> coll_pos {binary_list[icoll]->x_perp,binary_list[icoll]->y_perp,binary_list[icoll]->eta_source,binary_list[icoll]->tau_form}; 
      bool put_jet=TreeDoer(coll_pos, total_cross, icoll);
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

    double test_px=25.;
    double test_py=0.;
    double test_pz=0.;
    double test_en=25.;
    int test_id=1;
    //Smooth tests
    parton_list.push_back ( Parton ( test_px, test_py, test_pz, test_en, 0., 0., 0, -1, -1, test_id, "test", 0, 0, true ) );
    vector<double> test_pos {0.,0.,0.,0.};
    parton_list[parton_list.size()-1].vSetPos(test_pos);
    parton_list[parton_list.size()-1].vSetPosIn(test_pos);
}
bool Jets::EvolveJets(double tau, SCGrid &arena_current, const EOS &eos, hydro_source &hydro_source_terms, bool after_hydro_flag) {

    //cout << "EvolveJets:: with tau= " << tau << " and parton_list.size= " << parton_list.size() << endl;

    //Remove previous sources
    if (tau>=DATA.tau0) {
      cout << " Removing sources at tau= " << tau << endl; hydro_source_terms.remove_sources(tau);
    }

    bool something_happens=0;
    //Loop over parton list
    for (int iparton=0; iparton<parton_list.size(); iparton++) {
        //Skip if remnant
	if (parton_list[iparton].GetOrig()=="rem") continue;
	//Skip if completely quenched
	if (parton_list[iparton].vGetP()[3]==0.) continue;
	// If it is active
	if (parton_list[iparton].vGetPos()[3] != -1.) { 
	
	    if (parton_list[iparton].GetTauForm()<0. && after_hydro_flag==1) continue;
	    something_happens=1;

		
	    //if (after_hydro==1) cout << " at tau= " << tau << " splitTime= " << parton_list[iparton].splitTime() << " tau_form= "<< parton_list[iparton].GetTauForm() << " d1= " << parton_list[iparton].GetD1() << " mom= " << parton_list[iparton].GetMom() << endl;

	    //Set its splitTime if it is initial parton
            if (parton_list[iparton].splitTime()==-1. && parton_list[iparton].GetTauForm()!=-1.) parton_list[iparton].SetsplitTime();
	    
	    //cout << " splitTime = " << parton_list[iparton].splitTime() << " and tau form= " << parton_list[iparton].GetTauForm() << endl;

	    double sEta = parton_list[iparton].GetsEta();

	    //Compute energy loss and advance position
	    vector<double> p_bef=parton_list[iparton].vGetP();
	    DoEloss(parton_list[iparton], tau, arena_current, eos, hydro_source_terms);
	    vector<double> p_aft=parton_list[iparton].vGetP();
	    
	    //Compute source part
	    //if (tau+DATA.tau0>parton_list[iparton].GetTau()) cout << " Tau outside= " << parton_list[iparton].GetTau() << endl;

	    std::shared_ptr<jet> sourceterm(new jet);

	    sourceterm->tau_form=parton_list[iparton].GetTau();
	    sourceterm->x_perp=parton_list[iparton].vGetPos()[0];
	    sourceterm->y_perp=parton_list[iparton].vGetPos()[1];
	    if (p_aft[3]!=0.) sourceterm->eta_source=parton_list[iparton].GetsEta();
	    else sourceterm->eta_source=sEta;

	    sourceterm->dpxdtau=p_bef[0]-p_aft[0];
	    sourceterm->dpydtau=p_bef[1]-p_aft[1];
	    sourceterm->dpzdtau=p_bef[2]-p_aft[2];
	    sourceterm->dEdtau=p_bef[3]-p_aft[3];

	    if (p_bef[3]-p_aft[3]!=0.) hydro_source_terms.update_sources(sourceterm);
	    if (p_bef[3]-p_aft[3]!=0.) {
              for (unsigned a=0; a<4; a++) {
                injected_momentum[a]+=p_bef[a]-p_aft[a];
	      }
	    }
	    
	    //if (p_bef[3]-p_aft[3]!=0.) cout << " delta En= " << p_bef[3]-p_aft[3] << endl; 
	    //if (tau>=DATA.tau0) cout << setprecision(6) << "parton " << iparton << " deltaEn= " << p_bef[3]-p_aft[3] <<  " x= " << parton_list[iparton].vGetPos()[0] << " y= " << parton_list[iparton].vGetPos()[1] << endl;
	    //if (p_aft[3]==0.) cout << " ABSORBED !!! \n \n";
	    
	    double lambda=parton_list[iparton].vGetP()[3]/parton_list[iparton].vGetPIn()[3];
	    parton_list[iparton].SetLambda(lambda);	
	    
	    //Parton splits
	    if (parton_list[iparton].splitTime()<tau && parton_list[iparton].GetTauForm()!=-1.) {
            
		//cout << " \n \n Parton Split \n \n" << endl;
		
		int d1=parton_list[iparton].GetD1();
	        int d2=parton_list[iparton].GetD2();
	       
                parton_list[d1].set_hyper_point(parton_list[iparton].hyper_point());
                parton_list[d2].set_hyper_point(parton_list[iparton].hyper_point());

                parton_list[d1].set_last_temp(parton_list[iparton].last_temp());		
                parton_list[d2].set_last_temp(parton_list[iparton].last_temp());		

		parton_list[d1].vSetPos(parton_list[iparton].vGetPos());
		parton_list[d2].vSetPos(parton_list[iparton].vGetPos());
		
		parton_list[d1].vSetPosIn(parton_list[iparton].vGetPos());
		parton_list[d2].vSetPosIn(parton_list[iparton].vGetPos());
		
		parton_list[d1].vSetP(parton_list[d1].vGetP()*lambda);
		parton_list[d2].vSetP(parton_list[d2].vGetP()*lambda);
		
		parton_list[d1].vSetPinh(parton_list[d1].vGetP());
		parton_list[d2].vSetPinh(parton_list[d2].vGetP());
	
		parton_list[d1].SetsplitTime();
		parton_list[d2].SetsplitTime();

		//Set splitted parton as inactive
		vector<double> finpos {-1.,-1.,-1.,-1.};
		parton_list[iparton].vSetPos(finpos);
            
	    } 
	}

    } //parton list loop
              
    for (unsigned a=0; a<4; a++) {
      cout << " Injected " << a << " = " << injected_momentum[a] << endl;
    }

    return something_happens;
}
void Jets::DoEloss(Parton &parton, double tau, SCGrid &arena_current, const EOS &eos, hydro_source &hydro_source_terms) {
    
    double Tc=0.145;
    double kappa=0.41;

    double d_tau=DATA.delta_tau;

    //Check if done
    //if (parton.splitTime()<tau && parton.splitTime()!=-1.) {
      //  d_tau = tau-parton.splitTime(); 
      //  delta_t=d_tau;
    //}
    
    //Extract local fluid properties
    double x=parton.vGetPos()[0];
    double y=parton.vGetPos()[1];
    double rap=parton.GetsEta();

    double f_temp;
    double u0,u1,u2,u3;

    if (tau>=DATA.tau0) {
        double f_energy=GetFluidEnergy(x, y, rap, arena_current);
        const double rhob_center = 0.;
        f_temp=eos.get_temperature(f_energy, rhob_center)*hbarc;
        if (DATA.single_parton) cout << setprecision(6) << " local f_energy= " << f_energy << " f_temp= " << f_temp << endl; 
   
        double *v_flow = new double [4];
        GetFluidFlow(x, y, rap, arena_current, v_flow);
    
        u0=v_flow[0];
        u1=v_flow[1];
        u2=v_flow[3]*sinh(rap)+v_flow[2]*cosh(rap);
        u3=v_flow[3]*cosh(rap)+v_flow[2]*sinh(rap);
      
        delete v_flow;    
    }
    else {
        f_temp=0.;
        u0=0., u1=0., u2=0.;
        u3=1.;
    }

    vector<double> v;
    v.push_back(u0/u3);
    v.push_back(u1/u3);
    v.push_back(u2/u3);
    v.push_back(1.);

    vector<double> w;
    for (unsigned a=0; a<4; a++) w.push_back(parton.vGetP()[a]/parton.vGetP()[3]);

    double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
    double w2=pow(w[0],2.)+pow(w[1],2.)+pow(w[2],2.);
    double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
    if (v2>=1.) v2=0.999999999, cout << " V2 >= 1 \n";
    double lore=1./sqrt(1.-v2);

    //Color Factor
    double CF;
    if (parton.GetId()==21) CF=pow(9./4.,1./3.);        //If gluon, color charge dependence is ratio of casimirs to power 1/3
    else if (fabs(parton.GetId())<=6) CF=1.;
    else CF=0.;

    double new_tau=tau+d_tau; 
    double new_time=new_tau*cosh(rap);

    double delta_t=new_time-parton.vGetPos()[3];
    if (delta_t<0.) {cout << new_time << " " << parton.vGetPos()[3] << endl; exit(0); }

    //double prev_f_dist=parton.GetFluidDist();

    double delta_f_dist=delta_t*sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
    if (f_temp>=Tc) parton.AddFluidDist(delta_f_dist);

    double f_dist=parton.GetFluidDist();
    double ei=parton.vGetPinh()[3];

    vector<double> p_now=parton.vGetP();

    //Update momentum
    double quench=1.;
    if (parton.vGetP()[3]>0. && f_temp>=Tc && CF>0.) {
        if (kappa!=0.) {
            double Efs=ei*lore*(1.-vscalw);
            double tstop=0.2*pow(Efs,1./3.)/(2.*pow(f_temp,4./3.)*kappa)/CF;
            double beta=tstop/f_dist;
	    //double gamma=tstop/prev_f_dist;
	    //cout << " tstop= " << tstop << " f_dist= " << f_dist << " Efs= " << Efs << endl;
            if (beta>1.)
            {
                double intpiece=Efs*delta_t*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)));
                //double beta_intpiece=2.*Efs/3.141592/beta/beta*
			//(sqrt(beta*beta-1.)-beta*beta*atanh(1/sqrt(beta*beta-1.)));
                //double gamma_intpiece=2.*Efs/3.141592/gamma/gamma*
			//(sqrt(gamma*gamma-1.)-gamma*gamma*atanh(1/sqrt(gamma*gamma-1.)));
	        //if (prev_f_dist==0.) gamma_intpiece=0.; 
		//double intpiece=beta_intpiece-gamma_intpiece;
		//double quenched_energy = parton.vGetP()[3]*(1.+intpiece);
		
		//double quench=1.+intpiece;
		quench=(p_now[3]-intpiece)/p_now[3];
		if (quench<0.) quench=0.;
		parton.vSetP(parton.vGetP()*quench);                
            }
            else
            {
                vector<double> pzero (4, 0.);
	        parton.vSetP(pzero);
            }
        }
    }

    //Store position of hyper-surface crossing
    if (f_temp < Tc && parton.last_temp() >= Tc) {
      parton.set_hyper_point( x, y, rap, tau );
    }
    parton.set_last_temp(f_temp);

    if (tau>=DATA.tau0 && DATA.single_parton==1) {
    //if (tau>=DATA.tau0) {  
      cout << setprecision(6) << " px= " << parton.vGetP()[0] << " py= " << parton.vGetP()[1] << " pz= " << parton.vGetP()[2] << " en= " << parton.vGetP()[3] << " quench= " << quench << endl;
      cout << setprecision(6) << " delta en in= " << p_now[3]-parton.vGetP()[3] << endl;
      //cout << " px= " << p_now[0]*0.6 << " py= " << p_now[1]*0.6 << " pz= " << p_now[2]*0.6 << " en= " << p_now[3]*0.6 << " TOY " << endl;
    }

    //if (tau>=DATA.tau0) parton.vSetP(p_now);

    //Update position
    vector<double> p_parton;
    if (parton.vGetP()[3]!=0.) p_parton=parton.vGetP();
    else p_parton=p_now;

    vector<double> new_x = parton.vGetPos() +  (p_parton/p_parton[3]) * (new_time-parton.vGetPos()[3]);
    new_x[2]=new_tau*sinh(rap);  //should use new rap in case there is broadening, but watch out with completely quenched partons!

    parton.vSetPos(new_x);

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
		  HadronizeTherm();
		}
	    }
	}
    }
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
