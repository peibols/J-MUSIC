#include "jets.h"
#include "vector"
#include "vector_operators.h" 
#include <stdlib.h>
#include <time.h>

using namespace std;

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

            if (p_bef[3]-p_aft[3]!=0. && DATA.single_parton==1) {
              sourcefile << (*(sourceterm)).tau_form << " "
		      	 << (*(sourceterm)).x_perp << " "
			 << (*(sourceterm)).y_perp << " "
			 << (*(sourceterm)).eta_source << " "
			 << (*(sourceterm)).dpxdtau << " "
			 << (*(sourceterm)).dpydtau << " "
			 << (*(sourceterm)).dpzdtau << " "
			 << (*(sourceterm)).dEdtau << endl;
            }

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
      if (tau>=DATA.tau0) cout << " Injected " << a << " = " << injected_momentum[a] << endl;
    }

    return something_happens;
}
void Jets::DoEloss(Parton &parton, double tau, SCGrid &arena_current, const EOS &eos, hydro_source &hydro_source_terms) {
    
    double Tc=0.145;

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

    if (tau>=DATA.tau0 || DATA.prehydro_quenching==1) {
        double f_energy=GetFluidEnergy(x, y, rap, arena_current);
        const double rhob_center = 0.;
        f_temp=eos.get_temperature(f_energy, rhob_center)*hbarc;
        if (DATA.single_parton) {
	  cout << setprecision(6) << " local f_energy= " << f_energy << " f_temp= " << f_temp << endl; 
          cout << " x now= " << x << " y now= " << y << " eta now= " << rap << endl;
	}

	/*
	// CHECK
	for (int ix=-5; ix<=5; ix++) {
          for (int iy=-5; iy<=5; iy++) {
            cout << "MY INTERP= " << ix << " " << iy << " " << GetFluidEnergy(ix, iy, rap, arena_current) << endl;
	  }
        }
	int neta=arena_current.nEta();
	for (int ix=0; ix<300; ix++) {
          for (int iy=0; iy<300; iy++) {
            cout << "RAW= " << ix << " " << iy << " " << arena_current(ix,iy,neta/2).epsilon << endl;
	  }
        }
	*/

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
    if (parton.GetId()==21) {
      if (eloss_model==0) CF=pow(9./4.,1./3.);        //If gluon, color charge dependence is ratio of casimirs to power 1/3
      else CF=9./4.;
    }
    else if (fabs(parton.GetId())<=6) CF=1.;
    else CF=0.;

    //if (parton.GetId()==23) {
      //cout << " here it is the Z, with E= " << parton.vGetP()[3] << endl;
    //}

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
        if (kappa!=0. && eloss_model==0) { //Strong Coupling
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
	if (kappa!=0. && eloss_model==1) {  //Radiative
          double intpiece=(delta_t/0.2)*kappa*CF*f_temp*f_temp*f_temp*(f_dist/0.2);
          quench=(p_now[3]-intpiece)/p_now[3];
          if (quench<0.) quench=0.;
	  parton.vSetP(parton.vGetP()*quench);
	}
	if (kappa!=0. && eloss_model==2) {  //Collisional
          double intpiece=(delta_t/0.2)*kappa*CF*f_temp*f_temp;
          quench=(p_now[3]-intpiece)/p_now[3];
          if (quench<0.) quench=0.;
	  parton.vSetP(parton.vGetP()*quench);
	}
    }

    if (tau>=DATA.tau0) {
      //Store position of hyper-surface crossing
      if (f_temp < Tc && parton.last_temp() >= Tc) {
        parton.set_hyper_point( x, y, rap, tau );
      }
      parton.set_last_temp(f_temp);
    }

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

