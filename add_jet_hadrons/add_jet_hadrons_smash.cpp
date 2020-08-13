//Adds hadrons produced using lund string model on minijet partons
//from J-MUSIC to OSCAR hadrons list
//To compile the code:  g++ add_jet_hadrons.cpp

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>

typedef struct {
    int decay_Npart;
    double branching_ratio;
    int decay_part[5];
} decay_channel_info;

typedef struct {
    int monval;     // Monte Carlo number according PDG
    std::string name;
    double mass;
    double width;
    int gspin;      // spin degeneracy
    int baryon;
    int strange;
    int charm;
    int bottom;
    int gisospin;   // isospin degeneracy
    int charge;
    int decays;     // amount of decays listed for this resonance
    int stable;     // defines whether this particle is considered as stable
    std::vector<decay_channel_info*> decay_channels;
    int sign;                   // Bose-Einstein or Dirac-Fermi statistics
} particle_info;

typedef struct {
    double px;
    double py;
    double pz;
    double energy;
    double mass;
    double x;
    double y;
    double z;
    double t;
    int part_id;
    int charge;
} hadron_info;

int main() {
    int particle_list_size, num_of_jet_hadrons;
    std::vector<hadron_info> jet_hadrons;
    std::vector<particle_info> particle;
    void append_jet_hadrons_to_OSCAR(std::vector<hadron_info> &jet_hadrons, int num_of_jet_hadrons, std::vector<particle_info> particle, int particle_list_size);
    int read_jet_hadrons_list(std::vector<hadron_info> &jet_hadrons, std::vector<particle_info> &particle, int particle_list_size);
    int read_resonances_list(std::vector<particle_info> &particle);

    particle_list_size = read_resonances_list(particle);
    num_of_jet_hadrons = read_jet_hadrons_list(jet_hadrons, particle, particle_list_size);

    append_jet_hadrons_to_OSCAR(jet_hadrons, num_of_jet_hadrons, particle, particle_list_size);

}

void append_jet_hadrons_to_OSCAR(std::vector<hadron_info> &jet_hadrons, int num_of_jet_hadrons, std::vector<particle_info> particle, int particle_list_size) {
    std::cout << " -- Appending hadrons from J-MUSIC to soft hadrons...";
    std::string line;
    std::ifstream softhad("OSCAR_soft.dat");
    if (softhad.fail()) { std::cout << " no SoftHad File, exiting " << std::endl; exit(1); }
    std::string OSCAR_output_filename = "OSCAR.DAT";
    std::ofstream oscar(OSCAR_output_filename.c_str());

    int nevent=1;

    std::getline(softhad, line);
    std::getline(softhad, line);
    std::getline(softhad, line);
    
    oscar << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge" << std::endl;
    oscar << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none" << std::endl;

    int event_num, num_of_particles;
    double dummy = 0.;
    while (!softhad.eof()) {
	std::getline(softhad, line);
	if (softhad.eof()) break;
	std::istringstream iss(line);
        iss >> event_num >> num_of_particles >> dummy >> dummy;
	if (softhad.eof()) break;
	int total_number_of_particles = num_of_particles + num_of_jet_hadrons;
	oscar << "# event " << nevent << " out " << total_number_of_particles << std::endl;
	for (int idx = 0; idx < num_of_particles; idx++) {
            std::getline(softhad, line);
	    std::istringstream pss(line);
	    int ipart, id;
	    double px, py, pz, p0, mass, x, y, z, t;
	    pss >> ipart >> id >> px >> py >> pz >> p0 >> mass >>  x >> y >> z >> t;
	    int charge = -1000;
	    for (int idy = 0; idy <= particle_list_size; idy++) {
	      if (id == particle[idy].monval) {
                charge = particle[idy].charge;
	      }
	    }
	    if (charge == -1000) {
              std::cout << "WARNING: Unrecognized soft particle with pdg= " << id << std::endl;
	      continue;
	    }
	    oscar << t << " " << x << " " << y << " " << z << " "
		  << mass << " " << p0 << " " << px << " " << py << " " << pz << " "
		  << id << " " << ipart << " " << charge << std::endl;
	}
	for (int ipart = 0; ipart < num_of_jet_hadrons; ipart++) {
            int partnum = ipart + num_of_particles;
	    oscar << jet_hadrons[ipart].t << " " << jet_hadrons[ipart].x << " " << jet_hadrons[ipart].y << " " << jet_hadrons[ipart].z << " "
		  << jet_hadrons[ipart].mass << " " << jet_hadrons[ipart].energy << " " << jet_hadrons[ipart].px << " " << jet_hadrons[ipart].py << " " << jet_hadrons[ipart].pz << " "
		  << jet_hadrons[ipart].part_id << " " << partnum << " " << jet_hadrons[ipart].charge << std::endl;
	}
	oscar << "# event " << nevent << " end 0 impact   0.000" << std::endl;
	std::cout << "Done event= " << nevent << std::endl;
	nevent++;
    }
    softhad.close();
    oscar.close();
    std::cout << "done" << std::endl;
    return;
}

int read_jet_hadrons_list(std::vector<hadron_info> &jet_hadrons, std::vector<particle_info> &particle, int particle_list_size) {
    std::cout << " -- Reading in hadrons from J-MUSIC..." << std::endl;
    std::ifstream hadfile( "hadrons_list.dat");
    if (hadfile.fail()) { std::cout << " no HadronsList File, exiting " << std::endl; exit(1); }
    double px, py, pz, energy, x, y, z, t;
    int part_id;
    std::ofstream trigfile("trigger_particle.dat");
    while (!hadfile.eof()) {
        hadfile >> px;
        hadfile >> py;
        hadfile >> pz;
        hadfile >> energy;
        hadfile >> x;
        hadfile >> y;
        hadfile >> z;
        hadfile >> t;
        hadfile >> part_id;
	if (hadfile.eof()) break;
	if (part_id==23 || (part_id==22 && sqrt(px*px+py*py)>30.)) {
	  trigfile << px << py << pz << energy << part_id << std::endl;
	  continue;
	}
	for (int idx = 0; idx <= particle_list_size; idx++) {
	    if (part_id == particle[idx].monval) {
	        hadron_info hadron_i;
                hadron_i.px = px;
                hadron_i.py = py;
                hadron_i.pz = pz;
                hadron_i.energy = energy;
                hadron_i.x = x;
                hadron_i.y = y;
                hadron_i.z = z;
                hadron_i.t = t;
                hadron_i.part_id = part_id;
                hadron_i.mass = particle[idx].mass;
	        hadron_i.charge = particle[idx].charge;
		jet_hadrons.push_back(hadron_i);
		break;
	    } else if (idx == particle_list_size) {
		std::cout << "WARNING:: Unidentified particle with id " << part_id << std::endl;
	        
	    }
        }
    }
    hadfile.close();
    trigfile.close();
    std::cout << "done" << std::endl;
    return(jet_hadrons.size());
}

//Function to read resonance decay table
//Taken from iSS by Chun Shen and Zhi Qiu 
int read_resonances_list(std::vector<particle_info> &particle) {
    double eps = 1e-15;
    std::cout << " -- Read in particle resonance decay table...";
    std::ifstream resofile("pdg-SMASH.dat");
    if (resofile.fail()) { std::cout << " no Reso File, exiting " << std::endl; exit(1); }
    int local_i = 0;
    int dummy_int;
    while (!resofile.eof()) {
	particle_info particle_i;

        resofile >> particle_i.monval;
    	if (resofile.eof()) break;
        resofile >> particle_i.name;
        resofile >> particle_i.mass;
        resofile >> particle_i.width;
        resofile >> particle_i.gspin;        //spin degeneracy
        resofile >> particle_i.baryon;
        resofile >> particle_i.strange;
        resofile >> particle_i.charm;
        resofile >> particle_i.bottom;
        resofile >> particle_i.gisospin;     //isospin degeneracy
        resofile >> particle_i.charge;
        resofile >> particle_i.decays;
	if (resofile.eof()) break;
        for (int j = 0; j < particle_i.decays; j++) {
            decay_channel_info *temp_decay_channel = new decay_channel_info;
            resofile >> dummy_int;
            resofile >> temp_decay_channel->decay_Npart;
            resofile >> temp_decay_channel->branching_ratio;
            resofile >> temp_decay_channel->decay_part[0];
            resofile >> temp_decay_channel->decay_part[1];
            resofile >> temp_decay_channel->decay_part[2];
            resofile >> temp_decay_channel->decay_part[3];
            resofile >> temp_decay_channel->decay_part[4];
            particle_i.decay_channels.push_back(temp_decay_channel);
        }

        //decide whether particle is stable under strong interactions
        if (particle_i.decay_channels[0]->decay_Npart == 1) {
            particle_i.stable = 1;
        } else {
            particle_i.stable = 0;
        }

        if (!resofile.eof()) {
            particle.push_back(particle_i);
        } else {
            particle_i.baryon = 0;
        }

        //add anti-particle entry
        if (particle_i.baryon == 1) {
            local_i++;
            particle_info particle_j;
            particle_j.monval = -particle_i.monval;
	    std::ostringstream antiname;
            antiname << "Anti-" << particle_i.name;
            particle_j.name = antiname.str();
            particle_j.mass = particle_i.mass;
            particle_j.width = particle_i.width;
            particle_j.gspin = particle_i.gspin;
            particle_j.baryon = -particle_i.baryon;
            particle_j.strange = -particle_i.strange;
            particle_j.charm = -particle_i.charm;
            particle_j.bottom = -particle_i.bottom;
            particle_j.gisospin = particle_i.gisospin;
            particle_j.charge = -particle_i.charge;
            particle_j.decays = particle_i.decays;
            particle_j.stable = particle_i.stable;
            for (int j = 0; j < particle_j.decays; j++) {
                decay_channel_info *temp_anti_decay_channel = (
                                                    new decay_channel_info);
                temp_anti_decay_channel->decay_Npart = (
                        particle_i.decay_channels[j]->decay_Npart);
                temp_anti_decay_channel->branching_ratio =
                        particle_i.decay_channels[j]->branching_ratio;
                for (int k = 0; k < 5; k++) {
                    int decay_part_monval = (
                            particle_i.decay_channels[j]->decay_part[k]);
                    if (decay_part_monval == 0) {
                        // a null entry
                        temp_anti_decay_channel->decay_part[k] = 0;
                    } else {
                        // find the index for decay particle in the
                        // current resonance table
                        int idx;
                        for (idx = 0; idx < local_i; idx++) {
                            if (particle[idx].monval == decay_part_monval) {
                                break;
                            }
                        }
                        double temp_br = (
                            particle_i.decay_channels[j]->branching_ratio);
                        if (idx == local_i && particle_i.stable == 0
                            && temp_br > eps) {
			    std::cout << "Can not find decay particle index for "
                                      << "anti-baryon!" << std::endl;
			    std::cout << "particle monval : "
                                      << decay_part_monval << std::endl;
                            exit(1);
                        }
                        if (particle[idx].baryon == 0
                            && particle[idx].charge == 0
                            && particle[idx].strange == 0) {
                            temp_anti_decay_channel->decay_part[k] = (
                                particle_i.decay_channels[j]->decay_part[k]);
                        } else {
                            temp_anti_decay_channel->decay_part[k] = (
                                - particle_i.decay_channels[j]->decay_part[k]);
                        }
                    }
                }
                particle_j.decay_channels.push_back(temp_anti_decay_channel);
            }
            particle.push_back(particle_j);
        }
        local_i++;   // Add one to the counting variable "i" for the meson/baryon
    }
    for (auto &particle_i: particle) {
        if (particle_i.baryon == 0) {
            particle_i.sign = -1;
        } else {
            particle_i.sign=1;
        }
    }
    resofile.close();
    std::cout << "done." << std::endl;
    return(particle.size());
}

