// Copyright 2017 Chun Shen

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <memory>

#include "hydro_source.h"
#include "util.h"
#include <gsl/gsl_sf_bessel.h>

using std::string;
using std::cout;
using std::endl;

hydro_source::hydro_source(const InitData &DATA_in) :
    DATA(DATA_in) {
    source_tau_max = 0.0;
    source_tau_min = 100.0;
    if (DATA.Initial_profile == 13) {  // MC-Glauber-LEXUS
        sigma_tau            = 0.1;
        sigma_x              = 0.5;
        sigma_eta            = 0.5;
        volume               = DATA.delta_x*DATA.delta_y*DATA.delta_eta;
        string_dump_mode     = DATA.string_dump_mode;
        string_quench_factor = DATA.string_quench_factor;
        parton_quench_factor = 1.0;    // no diffusion current from the source
        read_in_QCD_strings_and_partons();
    }
    if (DATA.Initial_profile == 30) {  // AMPT
        sigma_tau = 0.1;
        sigma_x   = 0.5;
        sigma_eta = 0.2;
        volume = DATA.delta_x*DATA.delta_y*DATA.delta_eta;
        parton_quench_factor = 1.0;
        read_in_AMPT_partons();
    }

    // Options for reading jet energy deposition source terms - MS
    if (DATA.jet_medium_response_flag == 1 || DATA.jet_medium_response_flag == 2 || DATA.jet_medium_response_flag == 3) {
	jet_sigma_tau = DATA.jet_sigma_tau;
	jet_sigma_x   = DATA.jet_sigma_x;
	jet_sigma_eta = DATA.jet_sigma_eta;
	jet_volume    = DATA.delta_x*DATA.delta_y*DATA.delta_eta;
	if (DATA.jet_medium_response_flag == 1) read_in_jet_energy_loss();
	if (DATA.jet_medium_response_flag == 2) read_in_jet_energy_loss_2();
    }

    //Causal diffusion parameters
    tau_delay=0.5;
    time_relax=0.4;
    d_diff=0.08;
    width_delta=1.0;
    c_diff=sqrt(d_diff/time_relax);
    gamma_relax = 0.5/time_relax;
}

hydro_source::~hydro_source() {
    if (DATA.Initial_profile == 13) {
        QCD_strings_list.clear();
    }
    if (DATA.Initial_profile == 30) {
        parton_list.clear();
    }
    if (DATA.jet_medium_response_flag != 0 ) jet_energy_list.clear();
}

//! This function reads in the spatal information of the strings and partons
//! which are produced from the MC-Glauber-LEXUS model
void hydro_source::read_in_QCD_strings_and_partons() {
    string QCD_strings_filename = DATA.initName;
    string partons_filename = DATA.initName_rhob;
    music_message << "read in QCD strings list from " << QCD_strings_filename
                  << " and partons list from " << partons_filename;
    music_message.flush("info");
    string text_string;

    std::ifstream QCD_strings_file(QCD_strings_filename.c_str());
    if (!QCD_strings_file) {
        music_message << "hydro_source::read_in_QCD_strings_and_partons: "
                      << "can not open QCD strings file: "
                      << QCD_strings_filename;
        music_message.flush("error");
        exit(1);
    }
    getline(QCD_strings_file, text_string);  // read the header
    // now we read in data
    getline(QCD_strings_file, text_string);
    while (!QCD_strings_file.eof()) {
        std::stringstream text_stream(text_string);
        std::shared_ptr<QCD_string> new_string(new QCD_string);
        text_stream >> new_string->norm >> new_string->delta_E
                    >> new_string->tau_form
                    >> new_string->tau_0 >> new_string->eta_s_0
                    >> new_string->x_perp >> new_string->y_perp
                    >> new_string->eta_s_left >> new_string->eta_s_right
                    >> new_string->y_l >> new_string->y_r
                    >> new_string->frac_l >> new_string->frac_r
                    >> new_string->y_l_i;
        if (!text_stream.eof()) {
            // read in the last element
            text_stream >> new_string->y_r_i;
        } else {
            // the string is too short
            music_message << "read_in_QCD_strings_and_partons: "
                          << "the format of file"
                          << QCD_strings_filename << "is wrong~";
            music_message.flush("error");
            exit(1);
        }
        if (!text_stream.eof()) {
            // the string is too long
            music_message << "read_in_QCD_strings_and_partons: "
                          << "the format of file"
                          << QCD_strings_filename << "is wrong~";
            music_message.flush("error");
            exit(1);
        }

        double temp_factor1 = (new_string->tau_0*new_string->tau_0
                               - new_string->tau_form*new_string->tau_form);
        double temp_factor2 = (new_string->tau_0
                        *cosh(new_string->eta_s_left - new_string->eta_s_0));
        double temp_factor3 = (new_string->tau_0
                    *cosh(new_string->eta_s_right - new_string->eta_s_0));
        double tau_end_left_local = (
            temp_factor2 + sqrt(temp_factor2*temp_factor2 - temp_factor1));
        double tau_end_right_local = (
            temp_factor3 + sqrt(temp_factor3*temp_factor3 - temp_factor1));
        new_string->tau_end_left = tau_end_left_local;
        new_string->tau_end_right = tau_end_right_local;
        if (new_string->eta_s_left > new_string->eta_s_0) {
            new_string->tau_start = tau_end_left_local;
            new_string->eta_s_start = new_string->eta_s_left;
        } else if (new_string->eta_s_right < new_string->eta_s_0) {
            new_string->tau_start = tau_end_right_local;
            new_string->eta_s_start = new_string->eta_s_right;
        } else {
            new_string->tau_start = new_string->tau_0 + new_string->tau_form;
            new_string->eta_s_start = new_string->eta_s_0;
        }

        // read in one string properly
        QCD_strings_list.push_back(new_string);

        // record the proper time of the first and last string sources
        double source_tau = new_string->tau_form;
        if (new_string->tau_end_left > new_string->tau_end_right) {
            source_tau = new_string->tau_end_left;
        } else {
            source_tau = new_string->tau_end_right;
        }

        if (source_tau_max < source_tau) {
            source_tau_max = source_tau;
        }

        if (source_tau_min > (new_string->tau_0 + new_string->tau_form)) {
            source_tau_min = new_string->tau_0 + new_string->tau_form;
        }
        getline(QCD_strings_file, text_string);
    }
    QCD_strings_file.close();
    music_message << "hydro_source: tau_min = " << source_tau_min << " fm/c.";
    music_message.flush("info");
    music_message << "hydro_source: tau_max = " << source_tau_max << " fm/c.";
    music_message.flush("info");
    
    double total_baryon_number = 0;
    for (auto const& it : QCD_strings_list)
        total_baryon_number += it->frac_l + it->frac_r;
    music_message << "total baryon number = " << total_baryon_number;
    music_message.flush("info");
    compute_norm_for_strings();
}


//! This function reads in the partons information from the AMPT model
void hydro_source::read_in_AMPT_partons() {
    parton_list.clear();
    string AMPT_filename = DATA.initName_AMPT;
    music_message << "hydro_source: "
                  << "read in AMPT parton list from " << AMPT_filename;
    music_message.flush("info");

    string text_string;
    std::ifstream AMPT_file(AMPT_filename.c_str());
    if (!AMPT_file) {
        music_message << "hydro_source::read_in_AMPT_partons: "
                      << "can not open the AMPT file: " << AMPT_filename;
        music_message.flush("error");
        exit(1);
    }

    int n_partons = 0;
    int event_id  = 0;
    int dummy;
    getline(AMPT_file, text_string);
    std::stringstream text_stream1(text_string);
    text_stream1 >> event_id >> dummy >> n_partons;

    // now we read in data
    for (int ipart = 0; ipart < n_partons; ipart++) {
        getline(AMPT_file, text_string);
        std::stringstream text_stream(text_string);
        std::shared_ptr<parton> new_parton(new parton);
        double t_local, z_local, pz_local;
        int pid;
        text_stream >> pid >> new_parton->px >> new_parton->py >> pz_local
                    >> new_parton->mass
                    >> new_parton->x >> new_parton->y >> z_local >> t_local;
        if (t_local < z_local) continue;
        if (std::abs(pid) > 3) continue;

        // Now the parton is inside the light cone
        new_parton->E = sqrt(  new_parton->mass*new_parton->mass
                            + new_parton->px*new_parton->px
                            + new_parton->py*new_parton->py
                            + pz_local*pz_local);
        new_parton->tau      = sqrt(t_local*t_local - z_local*z_local);
        new_parton->eta_s    = 0.5*log( (t_local + z_local)
                                      /(t_local - z_local + 1e-15));
        new_parton->rapidity = 0.5*log( (new_parton->E + pz_local)
                                      /(new_parton->E - pz_local));
        double u_perp = (sqrt(  new_parton->px*new_parton->px
                              + new_parton->py*new_parton->py)
                         /new_parton->mass);
        new_parton->rapidity_perp = asinh(u_perp);
        if (pid == 1) {
            // d quark
            new_parton->baryon_number   =  1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge = -1./3.;
        } else if (pid == -1) {
            // anti-d quark
            new_parton->baryon_number   = -1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge =  1./3.;
        } else if (pid == 2) {
            // u quark
            new_parton->baryon_number   =  1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge =  2./3.;
        } else if (pid == -2) {
            // anti-u quark
            new_parton->baryon_number   = -1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge = -2./3.;
        } else if (pid == 3) {
            // s quark
            new_parton->baryon_number   =  1./3.;
            new_parton->strangness      = -1.0  ;
            new_parton->electric_charge = -1./3.;
        } else if (pid == -3) {
            // anti-s quark
            new_parton->baryon_number   = -1./3.;
            new_parton->strangness      =  1.0  ;
            new_parton->electric_charge =  1./3.;
        } else {
            cout << "pid = " << pid << endl;
        }
        parton_list.push_back(new_parton);
        if (source_tau_max < new_parton->tau) {
            source_tau_max = new_parton->tau;
        }
        if (source_tau_min > new_parton->tau) {
            source_tau_min = new_parton->tau;
        }
    }
    AMPT_file.close();
    music_message << "hydro_source:: read in " << parton_list.size() << "/"
                  << n_partons << " partons.";
    music_message.flush("info");
    music_message << "hydro_source:: tau_min = " << source_tau_min << " fm.";
    music_message.flush("info");
    music_message << "hydro_source:: tau_max = " << source_tau_max << " fm.";
    music_message.flush("info");
}


void hydro_source::prepare_list_for_current_tau_frame(double tau_local) {
    double dtau = DATA.delta_tau;
    QCD_strings_list_current_tau.clear();
    QCD_strings_baryon_list_current_tau.clear();
    parton_list_current_tau.clear();
    if (DATA.Initial_profile == 13) {
        for (auto &it: QCD_strings_list) {
            if ((   it->tau_end_left >= (tau_local - 1./2.*dtau)
                 && it->tau_end_left <  (tau_local + 3./2.*dtau))
                || (   it->tau_end_right >= (tau_local - 1./2.*dtau)
                    && it->tau_end_right <  (tau_local + 3./2.*dtau))) {
                QCD_strings_baryon_list_current_tau.push_back(it);
            }
            if (   it->tau_start <= tau_local + 3./2.*dtau
                && std::max(it->tau_end_left, it->tau_end_right) >= tau_local - dtau/2.) {
                QCD_strings_list_current_tau.push_back(it);
            }
        }
        music_message << "hydro_source: tau = " << tau_local << " fm."
                      << " number of strings for energy density: "
                      << QCD_strings_list_current_tau.size()
                      << " number of strings for net baryon density: "
                      << QCD_strings_baryon_list_current_tau.size();
        music_message.flush("info");
    } else if (DATA.Initial_profile == 30) {
        for (auto &it: parton_list) {
            double tau_dis = it->tau - tau_local;
            if (tau_dis > 0. && tau_dis < dtau) {
                parton_list_current_tau.push_back(it);
            }
        }
        music_message << "hydro_source: tau = " << tau_local
                      << " number of source: "
                      << parton_list_current_tau.size();
        music_message.flush("info");
    }
}

void hydro_source::get_hydro_energy_source(
    double tau, double x, double y, double eta_s, 
    FlowVec &u_mu, EnergyFlowVec &j_mu) {
    j_mu = {0};
    // flow velocity
    const double gamma_perp_flow = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow     = acosh(gamma_perp_flow);
    const double y_long_flow     = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sin_phi_flow    = u_mu[1]/gamma_perp_flow;
    const double cos_phi_flow    = u_mu[2]/gamma_perp_flow;
    const double dtau            = DATA.delta_tau;

    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_tau  = 1./dtau;
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double n_sigma_skip   = 5.;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;
    const double sfactor        = DATA.sFactor/hbarc;
    const double exp_tau = 1./tau;
    if (DATA.Initial_profile == 13) {
        for (auto const&it: QCD_strings_list_current_tau) {
            // energy source from strings
            const double tau_0     = it->tau_0;
            const double delta_tau = it->tau_form;
            
            double x_dis = x - it->x_perp;
            if (std::abs(x_dis) > skip_dis_x) continue;

            double y_dis = y - it->y_perp;
            if (std::abs(y_dis) > skip_dis_x) continue;

            // calculate the crossed string segments in the eta direction
            // normally, there will be two segments
            // [eta_L_next, eta_L] and [eta_R, eta_R_next]
            // the envelop profile for a segment [eta_L, eta_R] is
            // f(eta) = 0.5*(- Erf((eta_L - eta)/sigma)
            //               + Erf((eta_R - eta)/sigma))
            double eta_s_shift = 0.0;
            double tau_L = tau - dtau/2.;
            if (tau_L > tau_0 + delta_tau) {
                eta_s_shift = acosh((tau_L*tau_L + tau_0*tau_0
                                        - delta_tau*delta_tau)
                                       /(2.*tau_L*tau_0 + 1e-10));
            }
            double eta_s_L = std::min(it->eta_s_right, it->eta_s_0 - eta_s_shift);
            double eta_s_R = std::max(it->eta_s_left,  it->eta_s_0 + eta_s_shift);

            double eta_s_next_shift = 0.0;
            double tau_next = tau + dtau/2.;
            if (tau_next > tau_0 + delta_tau) {
                eta_s_next_shift = acosh((tau_next*tau_next + tau_0*tau_0
                                          - delta_tau*delta_tau)
                                         /(2.*tau_next*tau_0 + 1e-10));
            }
            double eta_s_L_next = std::max(it->eta_s_left,  it->eta_s_0 - eta_s_next_shift);
            double eta_s_R_next = std::min(it->eta_s_right, it->eta_s_0 + eta_s_next_shift);

            bool flag_left = true;  // the left string segment is valid
            if (eta_s_L_next > eta_s_L) flag_left = false;
            
            bool flag_right = true;  // the right string segment is valid
            if (eta_s_R_next < eta_s_R) flag_right = false;

            double exp_eta_s = 0.;
            if (flag_left) {
                if (   eta_s > eta_s_L_next - skip_dis_eta 
                    && eta_s < eta_s_L + skip_dis_eta) {
                    exp_eta_s += 0.5*(- erf((eta_s_L_next - eta_s)/sigma_eta)
                                      + erf((eta_s_L - eta_s)/sigma_eta));
                }
            }
            if (flag_right) {
                if (   eta_s > eta_s_R - skip_dis_eta 
                    && eta_s < eta_s_R_next + skip_dis_eta) {
                    exp_eta_s += 0.5*(- erf((eta_s_R - eta_s)/sigma_eta)
                                      + erf((eta_s_R_next - eta_s)/sigma_eta));
                }
            }

            double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                    /(sigma_x*sigma_x));

            double e_frac = 1.0;
            if (eta_s < it->eta_s_left) {
                e_frac = it->frac_l;
            } else if (eta_s < it->eta_s_right) {
                e_frac = (it->frac_l
                          + (it->frac_r - it->frac_l)
                            /(it->eta_s_right - it->eta_s_left)
                            *(eta_s - it->eta_s_left));
            } else {
                e_frac = it->frac_r;
            }
            double e_local = e_frac*exp_tau*exp_xperp*exp_eta_s;
            e_local *= it->norm*sfactor;  // 1/fm^4
            double y_string = (
                    it->y_l + (it->y_r - it->y_l)
                                /(it->eta_s_right - it->eta_s_left)
                                *(eta_s - it->eta_s_left));
            double y_dump = ((1. - string_quench_factor)*y_string
                             + string_quench_factor*y_long_flow);
            double y_dump_perp = string_quench_factor*y_perp_flow;
            double cosh_long = cosh(y_dump - eta_s);
            double sinh_long = sinh(y_dump - eta_s);
            double cosh_perp = 1.0;
            double sinh_perp = 0.0;
            if (std::abs(y_dump_perp) > 1e-6) {
                cosh_perp = cosh(y_dump_perp);
                sinh_perp = sinh(y_dump_perp);
            }
            j_mu[0] += e_local*cosh_long*cosh_perp;
            j_mu[1] += e_local*sinh_perp*cos_phi_flow;
            j_mu[2] += e_local*sinh_perp*sin_phi_flow;
            j_mu[3] += e_local*sinh_long*cosh_perp;
        }

        for (auto const&it: QCD_strings_baryon_list_current_tau) {
            // add baryon energy at the string ends
            bool flag_left = false;
            if (   it->tau_end_left >= tau - dtau/2.
                && it->tau_end_left <  tau + dtau/2.) {
                flag_left = true;
            }

            bool flag_right = false;
            if (   it->tau_end_right >= tau - dtau/2.
                && it->tau_end_right <  tau + dtau/2.) {
                flag_right = true;
            }
            
            double x_dis = x - it->x_perp;
            if (std::abs(x_dis) > skip_dis_x) continue;
            
            double y_dis = y - it->y_perp;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double exp_eta_s_left = 0.0;
            if (flag_left) {
                double eta_dis_left = std::abs(eta_s - it->eta_s_left);
                if (eta_dis_left < skip_dis_eta) {
                    exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                          /(sigma_eta*sigma_eta)));
                }
            }
            double exp_eta_s_right = 0.0;
            if (flag_right) {
                double eta_dis_right = std::abs(eta_s - it->eta_s_right);
                if (eta_dis_right < skip_dis_eta) {
                    exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                           /(sigma_eta*sigma_eta)));
                }
            }
            double exp_factors = exp_tau*(
                      exp_eta_s_left*it->frac_l*it->E_baryon_norm_L
                    + exp_eta_s_right*it->frac_r*it->E_baryon_norm_R);
            double e_baryon_local = 0.0;
            if (exp_factors > 0) {
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                e_baryon_local = exp_xperp*exp_factors;
            }
            j_mu[0] += e_baryon_local*sfactor;
        }
        double prefactors = prefactor_tau*prefactor_prep*prefactor_etas;
        j_mu[0] *= prefactors;
        j_mu[1] *= prefactors;
        j_mu[2] *= prefactors;
        j_mu[3] *= prefactors;
    } else if (DATA.Initial_profile == 30) {
        // AMPT parton sources
        double n_sigma_skip = 5.;
        double tau_dis_max = tau - source_tau_max;
        if (tau_dis_max < n_sigma_skip*sigma_tau) {
            for (auto &it: parton_list_current_tau) {
                double x_dis = x - it->x;
                if (std::abs(x_dis) > skip_dis_x) continue;

                double y_dis = y - it->y;
                if (std::abs(y_dis) > skip_dis_x) continue;

                double eta_s_dis = eta_s - it->eta_s;
                if (std::abs(eta_s_dis) > skip_dis_eta) continue;

                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double exp_eta_s = (
                        exp(-eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));

                double f_smear = exp_tau*exp_xperp*exp_eta_s;
                double p_perp_sq = it->px*it->px + it->py*it->py;
                double m_perp = sqrt(it->mass*it->mass + p_perp_sq);
                j_mu[0] += m_perp*cosh(it->rapidity - eta_s)*f_smear;
                j_mu[1] += it->px*f_smear;
                j_mu[2] += it->py*f_smear;
                j_mu[3] += m_perp*sinh(it->rapidity - eta_s)*f_smear;
            }
            double norm = DATA.sFactor/hbarc;     // 1/fm^4
            double prefactor = norm*prefactor_tau*prefactor_prep*prefactor_etas;
            j_mu[0] *= prefactor;
            j_mu[1] *= prefactor;
            j_mu[2] *= prefactor;
            j_mu[3] *= prefactor;
        }
    }
}

double hydro_source::get_hydro_rhob_source(double tau, double x, double y,
                                           double eta_s, FlowVec &u_mu) {
    double res = 0.;

    // flow velocity
    const double gamma_perp_flow  = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow      = acosh(gamma_perp_flow);
    const double y_long_flow      = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sinh_y_perp_flow = sinh(y_perp_flow);
    const double sin_phi_flow     = u_mu[1]/gamma_perp_flow;
    const double cos_phi_flow     = u_mu[2]/gamma_perp_flow;
    const double dtau             = DATA.delta_tau;

    const double exp_tau        = 1.0/tau;
    const double n_sigma_skip   = 5.;
    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double prefactor_tau  = 1./dtau;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;
    if (DATA.Initial_profile == 13) {
        for (auto &it: QCD_strings_baryon_list_current_tau) {
            // skip the evaluation if the strings is too far away in the
            // space-time grid
            // dumping energy into the medium from the active strings
            //double tau_dis_left = fabs(tau - it->tau_end_left);
            //double tau_dis_right = fabs(tau - it->tau_end_right);
            int flag_left = 0;
            if (   it->tau_end_left >= tau - dtau/2.
                && it->tau_end_left <  tau + dtau/2.) {
                flag_left = 1;
            }

            int flag_right = 0;
            if (   it->tau_end_right >= tau - dtau/2.
                && it->tau_end_right <  tau + dtau/2.) {
                flag_right = 1;
            }

            if (flag_left == 0 && flag_right == 0) continue;

            double x_dis = x - it->x_perp;
            if (std::abs(x_dis) > skip_dis_x) continue;
            
            double y_dis = y - it->y_perp;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double exp_eta_s_left = 0.0;
            if (flag_left == 1) {
                double eta_dis_left = std::abs(eta_s - it->eta_s_left);
                if (eta_dis_left < skip_dis_eta) {
                    exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                          /(sigma_eta*sigma_eta)));
                }
            }

            double exp_eta_s_right = 0.0;
            if (flag_right == 1) {
                double eta_dis_right = std::abs(eta_s - it->eta_s_right);
                if (eta_dis_right < skip_dis_eta) {
                    exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                           /(sigma_eta*sigma_eta)));
                }
            }
            
            double exp_factors = exp_tau*(
                    exp_eta_s_left*it->frac_l + exp_eta_s_right*it->frac_r);
            if (exp_factors > 0) {
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double fsmear = exp_xperp*exp_factors;
                double rapidity_local = (
                    (  exp_eta_s_left*it->frac_l*it->y_l
                     + exp_eta_s_right*it->frac_r*it->y_r)
                    /(  exp_eta_s_left*it->frac_l
                      + exp_eta_s_right*it->frac_r));
                double y_dump = ((1. - parton_quench_factor)*rapidity_local
                                 + parton_quench_factor*y_long_flow);
                double y_dump_perp = parton_quench_factor*y_perp_flow;
                double p_dot_u = 1.;
                if (parton_quench_factor < 1.) {
                    p_dot_u = (  u_mu[0]*cosh(y_dump)*cosh(y_dump_perp)
                               - u_mu[1]*sinh(y_dump_perp)*cos_phi_flow
                               - u_mu[2]*sinh(y_dump_perp)*sin_phi_flow
                               - u_mu[3]*sinh(y_dump)*cosh(y_dump_perp));
                }
                res += p_dot_u*fsmear;
            }
        }
        res *= prefactor_tau*prefactor_prep*prefactor_etas;
    } else if (DATA.Initial_profile == 30) {
        double tau_dis_max = tau - source_tau_max;
        if (tau_dis_max < n_sigma_skip*sigma_tau) {
            for (auto &it: parton_list_current_tau) {
                // skip the evaluation if the strings is too far away in the
                // space-time grid
                double x_dis = x - it->x;
                if (std::abs(x_dis) > skip_dis_x) continue;

                double y_dis = y - it->y;
                if (std::abs(y_dis) > skip_dis_x) continue;

                double eta_s_dis = eta_s - it->eta_s;
                if (std::abs(eta_s_dis) > skip_dis_eta) continue;

                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double exp_eta_s = (
                        exp(-eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));
                double f_smear = exp_tau*exp_xperp*exp_eta_s;
                double y_dump = ((1. - parton_quench_factor)*it->rapidity
                                 + parton_quench_factor*y_long_flow);
                double y_dump_perp = ((1. - parton_quench_factor)*it->rapidity_perp
                                      + parton_quench_factor*y_perp_flow);
                double p_dot_u = (u_mu[0]
                    - tanh(y_dump_perp)*sinh_y_perp_flow/cosh(y_dump - eta_s)
                    - tanh(y_dump - eta_s)*u_mu[3]);
                res += p_dot_u*f_smear;
            }
            res *= prefactor_tau*prefactor_prep*prefactor_etas;
        } }
    return(res);
}

void hydro_source::get_hydro_energy_source_before_tau(
    double tau, double x, double y, double eta_s, double *j_mu) {
    FlowVec u                   = {0};
    u[0]                        = 1.0;
    EnergyFlowVec j_mu_one_step = {0};

    double tau0 = 0.0;
    double dtau = DATA.delta_tau;
    int n_tau_steps = static_cast<int>((tau - tau0)/dtau);
    for (int i = 0; i < n_tau_steps; i++) {
        j_mu_one_step = {0};
        const double tau_local = tau0 + (i + 0.5)*dtau;
        get_hydro_energy_source(tau_local, x, y, eta_s, u, j_mu_one_step);
        for (int j = 0; j < 4; j++) {
            j_mu[j] += tau_local*j_mu_one_step[j]*dtau;
        }
    }
    for (int j = 0; j < 4; j++) {
        j_mu[j] /= tau;
    }
}

double hydro_source::get_hydro_rhob_source_before_tau(
        double tau, double x, double y, double eta_s) {
    FlowVec u = {0};
    u[0] = 1.0;

    double res  = 0.;
    double tau0 = 0.0;
    double dtau = DATA.delta_tau;

    int n_tau_steps = static_cast<int>((tau - tau0)/dtau);
    for (int i = 0; i < n_tau_steps; i++) {
        const double tau_local = tau0 + (i + 0.5)*dtau;
        const double res_local = get_hydro_rhob_source(tau_local, x, y, eta_s, u);
        res += tau_local*res_local*dtau;
    }

    return res/tau;
}

void hydro_source::compute_norm_for_strings() {
    const int neta              = 500;
    const double eta_range      = 6.;
    const double deta           = 2.*eta_range/(neta - 1);
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);

    double E_string_total   = 0.0;
    double E_baryon_total = 0.0;
    for (auto &it: QCD_strings_list) {
        double E_string_norm   = 0.;
        double E_baryon_L_norm = 0.;
        double E_baryon_R_norm = 0.;
        for (int ieta = 0; ieta < neta; ieta++) {
            double eta_local = - eta_range + ieta*deta;
            double f_eta = (it->frac_l
                + (it->frac_l - it->frac_r)/(it->eta_s_left - it->eta_s_right)
                  *(eta_local - it->eta_s_left));
            double y_eta = (it->y_l
                + (it->y_l - it->y_r)/(it->eta_s_left - it->eta_s_right)
                  *(eta_local - it->eta_s_left));

            double expon_left  = (it->eta_s_left - eta_local)/sigma_eta;
            double expon_right = (it->eta_s_right - eta_local)/sigma_eta;
            double e_eta = 0.5*(- erf(expon_left) + erf(expon_right));
            E_string_norm += f_eta*e_eta*cosh(y_eta);

            double e_baryon_L = exp(-expon_left*expon_left);
            double e_baryon_R = exp(-expon_right*expon_right);
            E_baryon_L_norm += e_baryon_L*cosh(eta_local);
            E_baryon_R_norm += e_baryon_R*cosh(eta_local);
        }
        E_string_norm   *= prefactor_etas*deta;
        double E_string = (  it->frac_l*cosh(it->y_l_i)
                           + it->frac_r*cosh(it->y_r_i)
                           - it->frac_l*cosh(it->y_l)
                           - it->frac_r*cosh(it->y_r));
        it->norm = E_string/E_string_norm;
        E_string_total += E_string;

        E_baryon_L_norm *= it->frac_l*prefactor_etas*deta;
        E_baryon_R_norm *= it->frac_r*prefactor_etas*deta;
        double E_baryon_L   = it->frac_l*cosh(it->y_l);
        double E_baryon_R   = it->frac_r*cosh(it->y_r);
        it->E_baryon_norm_L = E_baryon_L/E_baryon_L_norm;
        it->E_baryon_norm_R = E_baryon_R/E_baryon_R_norm;
        E_baryon_total += E_baryon_L + E_baryon_R;
    }
    music_message << "E_total = "
                  << (E_string_total + E_baryon_total)*DATA.sFactor << " GeV. "
                  << "E_string_total = " << E_string_total*DATA.sFactor
                  << " GeV" << ", E_baryon_total = "
                  << E_baryon_total*DATA.sFactor << " GeV.";
    music_message.flush("info");
}

// jet response read
void hydro_source::get_jet_energy_source(
    double time, double time_next, double x, double y, double eta_s, double *j_mu, double epsilon) {
    // clean up j_mu
    for (int i = 0; i < 4; i++) {
        j_mu[i] = 0.0;
    }
//    if (DATA_ptr->jet_response_flag == 1) {
        double n_sigma_skip = 5.;
            double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
            double prefactor_tau = 1./(DATA.delta_tau);
            double prefactor_eta = 1./(sqrt(M_PI)*sigma_eta*time);
	    double energy_check = 0.0;
            for (auto &it: jet_energy_list) {
		if((*it).tau_form < time) {
			music_message << "Some problem  " << time << "  " << (*it).tau_form;
			music_message.flush("error");
		}
		else{
                
                if ((*it).tau_form > time_next) {
		    break;
                }
                double x_dis = x - (*it).x_perp;
                if (fabs(x_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double y_dis = y - (*it).y_perp;
                if (fabs(y_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double eta_dis = eta_s - (*it).eta_source;
                if (fabs(eta_dis) > n_sigma_skip*sigma_eta) {
                    continue;
                }
                double exp_tau = 1.0;
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double exp_eta_s = exp(-(eta_dis*eta_dis)
                                        /(sigma_eta*sigma_eta));
                double e_local = exp_tau*exp_xperp*exp_eta_s;
		double j_in[4];

		j_in[0] = ((*it).dEdtau*cosh((*it).eta_source) - (*it).dpzdtau*sinh((*it).eta_source));
		j_in[1] = (*it).dpxdtau;
		j_in[2] = (*it).dpydtau;
		j_in[3] = (-(*it).dEdtau*sinh((*it).eta_source) + (*it).dpzdtau*cosh((*it).eta_source));
	// Lorentz transform to lab frame from fluid rest frame
	if(j_in[0] < 1.e-5 && j_in[1] < 1.e-5 && j_in[2] < 1.e-5 && j_in[3] < 1.e-5){
		for (int i = 0; i < 4; i++) {
	        	j_mu[i] += 0.0;
    		}
	}else{
	double u_mu[4];
	u_mu[0] = 1./sqrt(1. - ((*it).vx)*((*it).vx) - ((*it).vy)*((*it).vy) - ((*it).vz)*((*it).vz));
	u_mu[1] = ((*it).vx)*u_mu[0];
	u_mu[2] = ((*it).vy)*u_mu[0];
	u_mu[3] = ((*it).vz)*u_mu[0];

	double temp_u0 = (u_mu[0]*cosh((*it).eta_source) - u_mu[3]*sinh((*it).eta_source));
	double temp_u3 = (-u_mu[0]*sinh((*it).eta_source) + u_mu[3]*cosh((*it).eta_source));
	u_mu[0] = temp_u0;
	u_mu[3] = temp_u3;
	double j_out[4];
	Util::LorentzBoost4Vector(j_in,u_mu,j_out);

		for(int ii = 0; ii < 4; ii++){	
            		j_mu[ii] += j_out[ii]*e_local*prefactor_tau*prefactor_prep*prefactor_eta/hbarc;
		}

//			energy_check += j_out[0]*e_local*prefactor_prep*prefactor_eta/hbarc;
	}

	}
	}
	for(int ii = 0; ii < 4; ii++){	
		if(j_mu[ii] < 1.e-5) j_mu[ii] = 0.;
	}
return;
}

void hydro_source::RemoveSources(double time_now){
//	vector<jet> jetsource = jet_energy_list;
	jet_energy_list.erase(
		std::remove_if(jet_energy_list.begin(),jet_energy_list.end(), predicate_time(time_now)),
		jet_energy_list.end());
}


void hydro_source::CountSources(){
//	vector<jet> jetsource = jet_energy_list;
	int my_count = std::count_if(jet_energy_list.begin(),jet_energy_list.end(), predicate_count());
	music_message << "Number of source terms = " << my_count;	
	music_message.flush("info");
}


void hydro_source::read_in_jet_energy_loss() {
    string jet_filename = DATA.jetName;
    music_message << "read in jet energy loss from " << jet_filename;
    music_message.flush("info");
    string text_string;

    std::ifstream jet_file(jet_filename.c_str());
    if (!jet_file) {
        music_message << "Error:hydro_source::read_in_jet_energy_loss: "
                      << "can not open jet energy loss file: " << jet_filename;
	music_message.flush("error");
        exit(1);
    }
//    getline(jet_file, text_string);  // read the header
    // now we read in data
//    getline(jet_file, text_string);
    while (!jet_file.eof()) {
	std::stringstream text_stream(text_string);
        std::shared_ptr<jet> new_jetsource(new jet);
        text_stream >> new_jetsource->flag 
                    >> new_jetsource->tau_form
                    >> new_jetsource->x_perp >> new_jetsource->y_perp
                    >> new_jetsource->eta_source 
		    >> new_jetsource->dEdtau
                    >> new_jetsource->dpxdtau
		    >> new_jetsource->dpydtau
		    >> new_jetsource->dpzdtau
		    >> new_jetsource->vx
		    >> new_jetsource->vy
		    >> new_jetsource->vz;
        jet_energy_list.push_back(new_jetsource);
        getline(jet_file, text_string);
    }
    jet_file.close();

	std::sort(jet_energy_list.begin(), jet_energy_list.end(), sortPar);
//    cout << "hydro_source: tau_max = " << source_tau_max << endl;
    
}

bool hydro_source::sortPar(std::shared_ptr<jet>& x, std::shared_ptr<jet>& y)
     {
	return (x->tau_form < y->tau_form);
     }

// jet response read for causal diffusion
void hydro_source::get_causal_jet_energy_source(
    double time, double time_next, double x, double y, double eta_s, double *j_mu) {
    
    double dtau = DATA.delta_tau;
    // clean up j_mu
    for (int i = 0; i < 4; i++) {
        j_mu[i] = 0.0;
    }
            
    //if (jet_energy_list.size()==0 && time > DATA.tau0) cout << " Jet Energy List size= " << jet_energy_list.size() << endl;
    //if (time > DATA.tau0 && int((time-DATA.tau0)/dtau) != jet_energy_list.size()) cout << "int((time-DATA.tau0)/dtau)= " << int((time-DATA.tau0)/dtau) << " Jet Energy List size= " << jet_energy_list.size() << endl;
    //cout << " Jet Energy list size= " << jet_energy_list.size() << endl;
    //exit(0);
	    
    for (auto &it: jet_energy_list) {

	double tau_j = (*it).tau_form;    
	double x_j = (*it).x_perp;    
	double y_j = (*it).y_perp;    
	double eta_j = (*it).eta_source;    

	double ds2
	    = time*time + tau_j*tau_j
	    - 2.0*time*tau_j*cosh(eta_s-eta_j)
	    - (x-x_j)*(x-x_j)
	    - (y-y_j)*(y-y_j);

	if ( time>=tau_j && ds2 >= 0. ) {
	  
	    double t_j = tau_j * cosh(eta_j);
	    double t = time * cosh(eta_s);

	    double z_j = tau_j * sinh(eta_j);
	    double z = time * sinh(eta_s);

	    if ( time - 0.5*dtau <= tau_j + tau_delay &&
	        time + 0.5*dtau > tau_j + tau_delay ) {

	        double t_delay = tau_delay * cosh(eta_s);
		double delta_r2 = (x-x_j)*(x-x_j)+(y-y_j)*(y-y_j)+(z-z_j)*(z-z_j);

		double kernel = causal_diffusion_kernel(t_delay, sqrt(delta_r2))/dtau;
	  	
		//get jmu_i
	        double j_in[4];
		j_in[0] = kernel*((*it).dEdtau*cosh(eta_s) - (*it).dpzdtau*sinh(eta_s));
		j_in[1] = kernel*(*it).dpxdtau;
		j_in[2] = kernel*(*it).dpydtau;
		j_in[3] = kernel*(-(*it).dEdtau*sinh(eta_s) + (*it).dpzdtau*cosh(eta_s));
		//
                for(int ii = 0; ii < 4; ii++){	
		    j_mu[ii] += j_in[ii];
	        }
	    }
	}
    }
	    
    return;
}

double hydro_source::causal_diffusion_kernel(
		        double t, double r) const {

    double smooth = causal_diffusion_smooth(t, r);
    double delta = causal_diffusion_delta(t, r);
    
    return smooth+delta;

}

double hydro_source::causal_diffusion_smooth(double t, double r) const {
    

    if( r < c_diff*t ){
        
        
        double u = sqrt( c_diff*c_diff*t*t - r*r );
        double x = gamma_relax*u/c_diff; // unitless

        double i1 = gsl_sf_bessel_I1(x);
        double i2 = gsl_sf_bessel_In(2,x);

        return (exp(-gamma_relax*t)/(20.*M_PI))*(2.*gamma_relax*gamma_relax/c_diff)*(i1/(c_diff*u) + 4.*t*i2/u/u);
    }else{
        return 0.0;
    }
    
}

double hydro_source::causal_diffusion_delta(double t, double r) const {

    double r_w = width_delta;
    if( c_diff*t <= width_delta ){
        r_w = c_diff*t;
    }
    
    if( r >= c_diff*t - r_w &&
        r < c_diff*t ){
        
        return (exp(-gamma_relax*t)/(20.*M_PI))*(8. - 3.*exp(-gamma_relax*t) + 2.*gamma_relax*t +4.*gamma_relax*gamma_relax*t*t )/r/r/r_w;
        
    }else{
        return 0.0;
    }
    
}

// jet response read for Dani's Source terms
void hydro_source::get_jet_energy_source_2(
    double time, double time_next, double x, double y, double eta_s, double *j_mu) {
    // clean up j_mu
    for (int i = 0; i < 4; i++) {
        j_mu[i] = 0.0;
    }
            
            double n_sigma_skip = 5.;
            double prefactor_prep = 1./(M_PI*jet_sigma_x*jet_sigma_x);
            double prefactor_tau = 1./(DATA.delta_tau);
            double prefactor_eta = 1./(sqrt(M_PI)*jet_sigma_eta*time);
	    double energy_check = 0.0;
	    //if (jet_energy_list.size()==0 && time > DATA.tau0) cout << " Jet Energy List size= " << jet_energy_list.size() << endl;
	    
	    for (auto &it: jet_energy_list) {
		if((*it).tau_form < time-0.00001 && fabs(eta_s) < 0.01) {
			cout << std::scientific << time << "  " << (*it).tau_form;
			//music_message.flush("error");
		}
		else {
              	    //cout << " a source " << endl; 
		    if (fabs((*it).tau_form-time)>1.e-5) {
		        cout << " dif in time= " << (*it).tau_form-time << "\n \n";   	
		        exit(0);
		    }                    
		    if ((*it).tau_form > time_next) {
		      cout << " WTFFF = " << (*it).tau_form << " vs time_next " << time_next << endl;  
		      break;
                    }
                    double x_dis = x - (*it).x_perp;
                    if (fabs(x_dis) > n_sigma_skip*jet_sigma_x) {
                        continue;
                    }
                    double y_dis = y - (*it).y_perp;
                    if (fabs(y_dis) > n_sigma_skip*jet_sigma_x) {
                        continue;
                    }
                    double eta_dis = eta_s - (*it).eta_source;
                    if (fabs(eta_dis) > n_sigma_skip*jet_sigma_eta) {
                        continue;
                    }
		    //cout << " WE GOT SOMETHING !!! " << (*it).dEdtau << "\n \n";
                    double exp_tau = 1.0;
                    double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(jet_sigma_x*jet_sigma_x));
                    double exp_eta_s = exp(-(eta_dis*eta_dis)
                                        /(jet_sigma_eta*jet_sigma_eta));
                    double e_local = exp_tau*exp_xperp*exp_eta_s;
		    double j_in[4];
		    //j_in[0] = ((*it).dEdtau*cosh((*it).eta_source) - (*it).dpzdtau*sinh((*it).eta_source));
		    //j_in[1] = (*it).dpxdtau;
		    //j_in[2] = (*it).dpydtau;
		    //j_in[3] = (-(*it).dEdtau*sinh((*it).eta_source) + (*it).dpzdtau*cosh((*it).eta_source));
		    j_in[0] = ((*it).dEdtau*cosh(eta_s) - (*it).dpzdtau*sinh(eta_s));
		    j_in[1] = (*it).dpxdtau;
		    j_in[2] = (*it).dpydtau;
		    j_in[3] = (-(*it).dEdtau*sinh(eta_s) + (*it).dpzdtau*cosh(eta_s));

		    for (int ii = 0; ii < 4; ii++){	
            	        j_mu[ii] += j_in[ii]*e_local*prefactor_tau*prefactor_prep*prefactor_eta/hbarc;
			//cout << " at time= " << std::setprecision(6) << time << " J mu= " << j_mu[0] << " " << j_mu[1] << " " << j_mu[2] << " " << j_mu[3] << endl;
		        //cout << " at time= " << std::setprecision(6) << time << " J in= " << j_in[0] << " " << j_in[1] << " " << j_in[2] << " " << j_in[3] << endl;
		    }

                }
            }
	
	    for(int ii = 0; ii < 4; ii++){	
		if(fabs(j_mu[ii]) < 1.e-16) j_mu[ii] = 0.;
	    }
return;
}
void hydro_source::read_in_jet_energy_loss_2() {
    string jet_filename = DATA.jetName;
    cout << "read in jet energy loss from " << jet_filename << endl;
    string text_string;

    std::ifstream jet_file(jet_filename.c_str());
    if (!jet_file) {
        music_message << "Error:hydro_source::read_in_jet_energy_loss: "
                      << "can not open jet energy loss file: " << jet_filename;
        music_message.flush("error");
	exit(1);
    }
    // now we read in data
    getline(jet_file, text_string);
    getline(jet_file, text_string);
    getline(jet_file, text_string);
    getline(jet_file, text_string);
    while (!jet_file.eof()) {
	std::stringstream text_stream(text_string);
        std::shared_ptr<jet> new_jetsource(new jet);
        text_stream >> new_jetsource->tau_form
                    >> new_jetsource->x_perp >> new_jetsource->y_perp
                    >> new_jetsource->eta_source 
		    >> new_jetsource->dEdtau
                    >> new_jetsource->dpxdtau
		    >> new_jetsource->dpydtau
		    >> new_jetsource->dpzdtau;
        jet_energy_list.push_back(new_jetsource);
        getline(jet_file, text_string);
    }
    jet_file.close();
//    cout << "hydro_source: tau_max = " << source_tau_max << endl;
    std::sort(jet_energy_list.begin(), jet_energy_list.end(), sortPar);
}
void hydro_source::update_sources(std::shared_ptr<jet> newsource) {
    jet_energy_list.push_back(newsource);
    //cout << " Jet Energy List size= " << jet_energy_list.size() << endl;
    //for (unsigned int i=0; i<jet_energy_list.size(); i++) {
        //cout << "#" << i << " tau= " << jet_energy_list[i]->tau_form << endl;
    //}
}
void hydro_source::remove_sources(double time) {
    if (!DATA.causal_diffusion) jet_energy_list.clear();
    else {
	double dtau = DATA.delta_tau;
	auto it = jet_energy_list.begin();
	while (it != jet_energy_list.end())
	{
	    double tau_j=(*it)->tau_form;
            if ( time - 0.5*dtau > tau_j + tau_delay ||
                time + 0.5*dtau >= tau_j + tau_delay ) it = jet_energy_list.erase(it);
	    else ++it; 
	}

    }
}
