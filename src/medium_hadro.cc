#include "./jets.h"

using namespace std;

void Jets::ReadFreezeOutSurface() {
    music_message.info("reading freeze-out surface");
    
    ostringstream surfdat_stream;
    surfdat_stream << "./surface.dat";

    // new counting, mac compatible ...
    if (surface_in_binary) {
        NCells = get_number_of_lines_of_binary_surface_file(
                                                    surfdat_stream.str());
    } else {
        NCells = get_number_of_lines_of_text_surface_file(
                                                    surfdat_stream.str());
    }
    music_message << "NCells = " << NCells;
    music_message.flush("info");

    ifstream surfdat;
    if (surface_in_binary) {
        surfdat.open(surfdat_stream.str().c_str(), std::ios::binary);
    } else {
        surfdat.open(surfdat_stream.str().c_str());
    }
    // Now allocate memory: array of surfaceElements with length NCells
    //surface = (SurfaceElement *) malloc((NCells)*sizeof(SurfaceElement));
    int i = 0;
    while (i < NCells) {
        SurfaceElement temp_cell;
        if (surface_in_binary) {
            float array[32];
            for (int ii = 0; ii < 32; ii++) {
                float temp = 0.;
                surfdat.read((char*)&temp, sizeof(float));
                array[ii] = temp;
            }
            temp_cell.x[0] = array[0];
            temp_cell.x[1] = array[1];
            temp_cell.x[2] = array[2];
            temp_cell.x[3] = array[3];
            if (boost_invariant) {
                temp_cell.x[3] = 0.0;
            }

            temp_cell.s[0] = array[4];
            temp_cell.s[1] = array[5];
            temp_cell.s[2] = array[6];
            temp_cell.s[3] = array[7];

            temp_cell.u[0] = array[8];
            temp_cell.u[1] = array[9];
            temp_cell.u[2] = array[10];
            temp_cell.u[3] = array[11];

            temp_cell.epsilon_f            = array[12];
            temp_cell.T_f                  = array[13];
            temp_cell.mu_B                 = array[14];
            temp_cell.eps_plus_p_over_T_FO = array[15];
            
            temp_cell.W[0][0] = array[16];
            temp_cell.W[0][1] = array[17];
            temp_cell.W[0][2] = array[18];
            temp_cell.W[0][3] = array[19];
            temp_cell.W[1][1] = array[20];
            temp_cell.W[1][2] = array[21];
            temp_cell.W[1][3] = array[22];
            temp_cell.W[2][2] = array[23];
            temp_cell.W[2][3] = array[24];
            temp_cell.W[3][3] = array[25];

            temp_cell.pi_b  = array[26];
            temp_cell.rho_B = array[27];

            temp_cell.q[0] = array[28];
            temp_cell.q[1] = array[29];
            temp_cell.q[2] = array[30];
            temp_cell.q[3] = array[31];
        } else {
	    std::string ss;
	    getline(surfdat,ss);
	    istringstream ssurfdat(ss);
            // position in (tau, x, y, eta)
            ssurfdat >> temp_cell.x[0] >> temp_cell.x[1]
                    >> temp_cell.x[2] >> temp_cell.x[3];
            
            // hypersurface vector in (tau, x, y, eta)
            ssurfdat >> temp_cell.s[0] >> temp_cell.s[1]
                    >> temp_cell.s[2] >> temp_cell.s[3];
            
            // flow velocity in (tau, x, y, eta)
            ssurfdat >> temp_cell.u[0] >> temp_cell.u[1]
                    >> temp_cell.u[2] >> temp_cell.u[3];

            ssurfdat >> temp_cell.epsilon_f >> temp_cell.T_f
                    >> temp_cell.mu_B >> temp_cell.eps_plus_p_over_T_FO;

            // freeze-out Wmunu
            ssurfdat >> temp_cell.W[0][0] >> temp_cell.W[0][1]
                    >> temp_cell.W[0][2] >> temp_cell.W[0][3]
                    >> temp_cell.W[1][1] >> temp_cell.W[1][2]
                    >> temp_cell.W[1][3] >> temp_cell.W[2][2]
                    >> temp_cell.W[2][3] >> temp_cell.W[3][3];
            if (DATA.turn_on_bulk) {
                ssurfdat >> temp_cell.pi_b;
            } else {
                temp_cell.pi_b = 0.;
            }
            if (DATA.turn_on_rhob) {
                ssurfdat >> temp_cell.rho_B;
            } else {
                temp_cell.rho_B = 0.;
            }
            if (DATA.turn_on_diff) {
                ssurfdat >> temp_cell.q[0] >> temp_cell.q[1]
                        >> temp_cell.q[2] >> temp_cell.q[3];
            } else {
                temp_cell.q[0] = 0.;
                temp_cell.q[1] = 0.;
                temp_cell.q[2] = 0.;
                temp_cell.q[3] = 0.;
            }
        }
        temp_cell.sinh_eta_s = sinh(temp_cell.x[3]);
        temp_cell.cosh_eta_s = cosh(temp_cell.x[3]);

        if (temp_cell.epsilon_f < 0)  {
	    cout << " lets seeeee= " << temp_cell.epsilon_f << endl;
            music_message.error("epsilon_f < 0.!");
            exit(1);
        }
        if (temp_cell.T_f < 0) {
            music_message.error("T_f < 0.!");
            exit(1);
        }
        surface.push_back(temp_cell);
        i++;
    }
    surfdat.close();
}


int Jets::get_number_of_lines_of_binary_surface_file(string filename) {
    std::ifstream surface_file(filename.c_str(), std::ios::binary);
    int count = 0;
    float temp = 0.;
    while(surface_file) {
        surface_file.read((char*) &temp, sizeof(float));
        count++;
    }
    int counted_line = count/32;
    surface_file.close();
    return(counted_line);
}


int Jets::get_number_of_lines_of_text_surface_file(string filename) {
    std::ifstream surface_file(filename.c_str(), std::ios::binary);
    int counted_lines = 0;
    std::string temp_line;
    while (std::getline(surface_file, temp_line)) {
        ++counted_lines;
    }
    surface_file.close();
    return(counted_lines);
}

void Jets::SampleSurface(Parton parton)
{
    fin_and_therm_parton_list.clear();

    double tau_parton = parton.hyper_point()[3];
    double eta_s_parton = parton.hyper_point()[2];
    double x_parton = parton.hyper_point()[0];
    double y_parton = parton.hyper_point()[1];

    if (tau_parton == -1) {
      cout << " Unassigned hyper-surface point! " << endl;
      return;
    }
    
    int p_col = parton.GetCol();
    int p_acol = parton.GetAcol();

    if ( p_col==0 && p_acol==0 ) {
      fin_and_therm_parton_list.push_back( parton );
      return;
    }

    double ptmax = DATA.max_pt;
    double ptmin = DATA.min_pt;
    int iptmax = DATA.pt_steps+1;  // Number of points is iptmax + 1
    
    double y_minus_eta_cut = 4.0;
    double etamax = DATA.max_pseudorapidity;
    int ietamax = DATA.pseudo_steps + 1;  // pseudo_steps is number of steps.
    double deltaeta = 0;
    if (ietamax > 1) {
        deltaeta = 2.*etamax/DATA.pseudo_steps;
    }
    
    int iphimax = DATA.phi_steps;  // number of points
    double deltaphi = 2*PI/iphimax;

    double spectrum[ietamax][iptmax][iphimax]={{{0.}}};
    double norm=0.;

    //double m = particleList[j].mass;
    //int deg = particleList[j].degeneracy;

    double mu_PCE = 0.; // not sure?
    double baryon = 1./3.;
    double m = 0.33;  //just sampling u and d quarks for now
    int sign = -1;

    // caching 
    double *cos_phi = new double[iphimax];
    double *sin_phi = new double[iphimax];
    for (int iphi = 0; iphi < iphimax; iphi++) {
        double phi_local = deltaphi*iphi;
        cos_phi[iphi] = cos(phi_local);
        sin_phi[iphi] = sin(phi_local);
    }
    double* pt_array = new double[iptmax];
    for (int ipt = 0; ipt < iptmax; ipt++) {
        double pt =  (ptmin + (ptmax - ptmin)
                              *pow(static_cast<double>(ipt), 2.)
                              /pow(static_cast<double>(iptmax - 1), 2.));
        pt_array[ipt] = pt;
    }

    // Look for nearby cell
    int the_cell = -1000;
    for (int icell = 0; icell < NCells; icell++) {
      double tau_cell     = surface[icell].x[0];
      double eta_s_cell   = surface[icell].x[3];
      double x_cell       = surface[icell].x[1];
      double y_cell       = surface[icell].x[2];

      if ( fabs(tau_cell - tau_parton) < tau_tol &&
	   fabs(eta_s_cell - eta_s_parton) < eta_tol &&
           fabs(x_cell - x_parton) < x_tol &&	   
           fabs(y_cell - y_parton) < y_tol ) {
        the_cell = icell;
	break;
      }
    }

    if ( the_cell == -1000 ) {
      cout << " Could not find nearby cell! " << endl;
      return;
    }

    // Cell properties
    int icell = the_cell;
                
    double tau        = surface[icell].x[0];
    double eta_s      = surface[icell].x[3];
    double cosh_eta_s = surface[icell].cosh_eta_s;
    double sinh_eta_s = surface[icell].sinh_eta_s;

    double T   = surface[icell].T_f*hbarc;  // GeV
    double muB = surface[icell].mu_B*hbarc;  // GeV
    double mu  = baryon*muB;  // GeV
    if (DATA.whichEOS>=3 && DATA.whichEOS < 10) {
	    // for PCE use the previously computed mu
	    // at the freeze-out energy density
	    mu += mu_PCE;  // GeV
    }
    
    double sigma_mu[4];
    double u_flow[4];
    for (int ii = 0; ii < 4; ii++) {
      sigma_mu[ii] = surface[icell].s[ii];
      u_flow[ii] = surface[icell].u[ii];
    }

    // Just ideal cooper-frye for now (improve!)
    double W00 = 0.0;  
    double W01 = 0.0;
    double W02 = 0.0;
    double W03 = 0.0;
    double W11 = 0.0;
    double W12 = 0.0;
    double W13 = 0.0;
    double W22 = 0.0;
    double W23 = 0.0;
    double W33 = 0.0;
    int flag_shear_deltaf = 0;
   
    double Pi_bulk = 0.0;
    int flag_bulk_deltaf = 0;

    double qmu_0 = 0.0;
    double qmu_1 = 0.0;
    double qmu_2 = 0.0;
    double qmu_3 = 0.0;
    double deltaf_qmu_coeff = 1.0;
    double deltaf_qmu_coeff_14mom_DV = 0.0;
    double deltaf_qmu_coeff_14mom_BV = 0.0;
    int flag_qmu_deltaf = 0;
    
    double rhoB = 0.0;
    
    double eps_plus_P_over_T = surface[icell].eps_plus_p_over_T_FO;
    double prefactor_shear = 1./(2.*eps_plus_P_over_T*T*T*T)*hbarc; 
    // fm^4/GeV^2
    double prefactor_qmu = rhoB/(eps_plus_P_over_T*T);   // 1/GeV    
    //

    double alpha=0.;
    int ieta;
    for (ieta = 0; ieta < ietamax; ieta++) {
      double eta = -etamax + ieta*deltaeta;
        
      double* rapidity = new double [iptmax];
      double* cosh_y = new double [iptmax];
      double* sinh_y = new double [iptmax];
      
      for (int ipt = 0; ipt < iptmax; ipt++) {
        double pt = pt_array[ipt];
        double y_local;
        // rapidity as a function of pseudorapidity:
        if (DATA.pseudofreeze == 1) {
          y_local = Rap(eta, pt, m);
        } else {
          y_local = eta;
        }
        rapidity[ipt] = y_local;
        cosh_y[ipt] = cosh(y_local);
        sinh_y[ipt] = sinh(y_local);
      }
      
      // Build spectrum      
      for (int ipt = 0; ipt < iptmax; ipt++) {
        double y = rapidity[ipt];    
      	if (fabs(y - eta_s) < y_minus_eta_cut) {
          double pt = pt_array[ipt];
          double cosh_y_local = cosh_y[ipt];  
	  double sinh_y_local = sinh_y[ipt];
	  double mt = sqrt(m*m + pt*pt);     // all in GeV

	  double ptau = mt*(cosh_y_local*cosh_eta_s
				- sinh_y_local*sinh_eta_s); 
	  double peta = mt*(sinh_y_local*cosh_eta_s
		        	- cosh_y_local*sinh_eta_s); 

          for (int iphi = 0; iphi < iphimax; iphi++) {
            double px = pt*cos_phi[iphi];
    	    double py = pt*sin_phi[iphi];
    
	    double sum;
            // compute p^mu*dSigma_mu [fm^3*GeV]
                    
	    double pdSigma = tau*(ptau*sigma_mu[0]
		  		+ px*sigma_mu[1]
	      			+ py*sigma_mu[2]
	       		        + peta/tau*sigma_mu[3]);
	    double E = (ptau*u_flow[0] - px*u_flow[1]
			- py*u_flow[2] - peta*u_flow[3]);
                  
            // this is the equilibrium f, f_0:        
	    double f = 1./(exp(1./T*(E - mu)) + sign);
                     
            double delta_f_shear = 0.0;
            double delta_f_bulk = 0.0;
                            
            // delta f for qmu
            double qmufactor = 0.0;
            double delta_f_qmu = 0.0;

            double max_ratio = 1.0;
            double total_deltaf = (delta_f_shear
                                 + delta_f_bulk
                                 + delta_f_qmu);

            if (fabs(total_deltaf)/f > max_ratio) {
              total_deltaf *= f/fabs(total_deltaf);
            }
            sum = (f + total_deltaf)*pdSigma;
                     
            if (sum > 10000) {
              music_message << "sum>10000 in summation. sum = "
                            << sum 
                            << ", f=" << f << ", deltaf="
                            << delta_f_shear
                            << ", pdSigma=" << pdSigma << ", T=" << T 
                            << ", E=" << E << ", mu=" << mu;
              music_message.flush("warning");
            }
            if (f < 0.) {
              music_message << " f_eq < 0.! f_eq = " << f
                            << ", T = " << T << " GeV, mu = "
                            << mu << " GeV, E = "
                            << E << " GeV";
              music_message.flush("error");
            }
            
	    spectrum[ieta][ipt][iphi] += sum;
      	    norm += sum;                  
	  }
        }
      }
    
      delete[] rapidity;
      delete[] cosh_y;
      delete[] sinh_y;

    } // eta loop

    // Fill final parton
    fin_and_therm_parton_list.push_back( parton );

    for (int icol=0; icol<2; icol++) {
	    
      int the_col;
      if ( icol==0 ) the_col = p_col;
      else the_col = p_acol;
      
      if ( the_col==0) continue;

      int id;
      int t_col, t_acol;

      if ( rand() > 0.5 ) id = 1;	// u quark
      else id = 2;			// d quark
      
      if ( icol==0 ) id *= -1, t_col = 0, t_acol = the_col;
      else t_col = the_col, t_acol = 0;
      

      // MonteCarlo
      double nrand = 1.;
      double dist = 0.;
      int ipt, iphi;
      do {
        ieta = int(ietamax * rand());
        ipt = int(iptmax * rand());
        iphi = int(iphimax * rand());

        dist = spectrum[ieta][ipt][iphi] / norm;
        nrand = rand();
      } while ( nrand > dist );

      double pt = pt_array[ipt];
      double mt = sqrt(pt*pt + m*m);
      double px = pt*cos_phi[iphi];
      double py = pt*sin_phi[iphi];
      double eta = -etamax + ieta*deltaeta;
      double pz = mt*sinh(eta);
      double en = mt*cosh(eta);

      // Fill thermal parton
      fin_and_therm_parton_list.push_back ( Parton ( px, py, pz, en, 0., m, 0, -1, -1, id, "therm", t_col, t_acol, true ) );

      negafile << px << " " << py << " " << pz << " " << en << " " << id << " "
	          x_cell << " " << y_cell << " " << eta_s_cell << " " << tau_cel << endl;  

    }

    // Clean up
    delete[] cos_phi;
    delete[] sin_phi;
    delete[] pt_array;
    
    return;

}

double Jets::Rap(double eta, double pt, double m) {
    double y = log((sqrt(m*m + pt*pt*cosh(eta)*cosh(eta)) + pt*sinh(eta))
                   /sqrt(m*m+pt*pt));
    return y;
}

