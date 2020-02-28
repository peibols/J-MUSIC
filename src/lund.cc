#include "jets.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

using namespace Pythia8;

void Jets::InitLund() { 
    // No event record printout.
    hpythia.readString("Next:numberShowInfo = 0");
    hpythia.readString("Next:numberShowProcess = 0");
    hpythia.readString("Next:numberShowEvent = 0");
    
    hpythia.readString("ProcessLevel:all = off");
   
    hpythia.readString("StringFragmentation:TraceColours = on"); 
    hpythia.readString("Fragmentation:setVertices = on");

    // Handle hadron decays to match limited PDG list from UrQMD
    hpythia.readString("HadronLevel:Decay = on");
    HandleDecays();

    // And initialize
    hpythia.init();

    hadron_list.clear();
}

void Jets::HandleDecays() {

    ifstream urqmd_pdg("./EOS/pdg-urqmd_v3.3+.dat");
    if (urqmd_pdg.fail()) {
      cout << "No PDG file for UrQMD!" << endl;
      exit(0);
    }

    string s;
    while (true) {
      getline(urqmd_pdg,s);
      if (urqmd_pdg.eof()) break;
      istringstream iss(s);
      int pdg_id;
      iss >> pdg_id;
      if (!hpythia.particleData.isParticle(pdg_id)) continue;
      ostringstream oss;
      oss << pdg_id << ":mayDecay = off";
      hpythia.readString(oss.str());
    }

}

void Jets::HadronizeTherm() {

    int max_attempts=5;

    int attempts=0;
    while (true) {	
      Event& event      = hpythia.event;
      ParticleData& pdt = hpythia.particleData;

      event.reset();

      vector<Parton> pIn = fin_and_therm_parton_list;
      //cout << " fin_and_therm_parton_list size= " << fin_and_therm_parton_list.size() << endl;

      vector<double> hyper_point;

      for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
      {  
          vector<double> temp_hyper_point = pIn[ipart].hyper_point();
          if (temp_hyper_point[3]!=-1.) hyper_point = temp_hyper_point;	   

  	  int ide=pIn[ipart].GetId();
          double px=pIn[ipart].vGetP()[0];
          double py=pIn[ipart].vGetP()[1];
          double pz=pIn[ipart].vGetP()[2];
          double ee=pIn[ipart].vGetP()[3];
          double mm=pdt.m0(int(ide));
  	  int col=pIn[ipart].GetCol();
	  int acol=pIn[ipart].GetAcol();
          ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
          if (col==0 && acol==0 && (ide==21 || abs(ide)<=6)) {
              cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
              exit(0);
          }
          event.append(int(ide),23,col,acol,px,py,pz,ee,mm);
      }
    
      if (hyper_point.size()==0) {
        cout << " Hyper size is 0 in Lund! " << endl;
        return;
      }

      double fin_pos[4];
      fin_pos[0] = hyper_point[0];
      fin_pos[1] = hyper_point[1];
      fin_pos[2] = hyper_point[3]*sinh(hyper_point[2]);
      fin_pos[3] = hyper_point[3]*cosh(hyper_point[2]);

      //cout << " fin_pos= " << fin_pos[0] << " " << fin_pos[1] << " " << fin_pos[2] << " " << fin_pos[3] << " " << endl;

      if (!hpythia.next()) {
        attempts++;
	if (attempts>=max_attempts) break;	
	else continue;
      }
      //hadfile <<"weight " << event_weight << " cross " << event_cross << endl;
    
      for (unsigned int ipart=0; ipart < event.size(); ++ipart)
      {
          if (event[ipart].isFinal())
          {
              int ide=event[ipart].id();
              vector<double> p {event[ipart].px(),event[ipart].py(),event[ipart].pz(),event[ipart].e()};
              hadron_list.push_back ( Parton ( p, 0., event[ipart].m(), 0, -1, -1, ide, "hadron", 0, 0, true ) );
    
              double had_pos[4];
	      had_pos[0] = fin_pos[0] + event[ipart].xProd()*MM2FM;	    
	      had_pos[1] = fin_pos[1] + event[ipart].yProd()*MM2FM;	    
	      had_pos[2] = fin_pos[2] + event[ipart].zProd()*MM2FM;	    
	      had_pos[3] = fin_pos[3] + event[ipart].tProd()*MM2FM;
	      //cout << " xprod= " << event[ipart].xProd()*MM2FM
		//   << " yprod= " << event[ipart].yProd()*MM2FM
		//   << " zprod= " << event[ipart].zProd()*MM2FM
		//   << " tprod= " << event[ipart].tProd()*MM2FM << endl;

              hadron_list[hadron_list.size()-1].SetPos(had_pos[0],had_pos[1],had_pos[2],had_pos[3]);

	      //Print on output file
              hadfile << event[ipart].px() << " " << event[ipart].py() << " " << event[ipart].pz() << " " << event[ipart].e() << " " 
		    << had_pos[0] << " " << had_pos[1] << " " << had_pos[2] << " " << had_pos[3] << " "
		    << event[ipart].id() << endl;
          }  
      }

      break;

    }

}

void Jets::HadronizeCorona() {
    
    Event& event      = hpythia.event;
    ParticleData& pdt = hpythia.particleData;

    double Lambda_QCD=0.2;
    double rempx=0.2;
    double rempy=0.2;
    double p_fake=1380.;
    double rempz=p_fake;
    double reme=std::sqrt(std::pow(rempx,2.)+std::pow(rempy,2.)+std::pow(rempz,2.));
    
    // Hadronize all showers together
    vector<Parton> pIn = corona_parton_list;
    cout << " Corona list size= " << pIn.size() << endl;

    // Check whether event is empty
    if (pIn.size()==0) return;

    int col[pIn.size()+2], acol[pIn.size()+2], isdone[pIn.size()+2];
    memset( col, 0, (pIn.size()+2)*sizeof(int) ), memset( acol, 0, (pIn.size()+2)*sizeof(int) ), memset( isdone, 0, (pIn.size()+2)*sizeof(int) );
  
    // Find number of quarks
    int nquarks=0;
    int isquark[pIn.size()+2];
    memset( isquark, 0, (pIn.size()+2)*sizeof(int) );
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
        if (abs(pIn[ipart].GetId())<=6) {
            isquark[nquarks]=ipart;
            nquarks+=1;  
        }
    }
  
    // Find number of strings
    int nstrings=max(int(double(nquarks)/2.+0.6),1);
    
    // If there are no quarks, need to attach two of them
    int istring=0;
    int one_end[nstrings], two_end[nstrings];
    if (nquarks==0) { // Only attach remnants if event is not empty
        // First quark
        pIn.push_back ( Parton ( rempx, rempy, rempz, reme, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true ) );
        isquark[nquarks]=pIn.size()-1;
        nquarks+=1;
        isdone[pIn.size()-1]=1;
        one_end[0]=pIn.size()-1;
        // Second quark
        pIn.push_back ( Parton ( rempx, rempy, -rempz, reme, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true ) );
        isquark[nquarks]=pIn.size()-1;
        nquarks+=1;
        isdone[pIn.size()-1]=1;
        two_end[istring]=pIn.size()-1;
    }

    // Assign ends of strings (order matters in this algo)
    for(unsigned int iquark=0; iquark<nquarks; iquark++) {
        if (isdone[isquark[iquark]]==0) {
            isdone[isquark[iquark]]=1;
            one_end[istring]=isquark[iquark];
            double min_delR=1000000.;
            int partner=-2;
            for(unsigned int jquark=0; jquark<nquarks; jquark++) {  
                if (iquark==jquark) continue;
                int d_jquark=isquark[jquark];
                if (isdone[d_jquark]==0) {
                    double delR = pIn[isquark[iquark]].delta_R(pIn[d_jquark]);
	            if (delR<min_delR) min_delR=delR, partner=jquark;
                }
            }
            if (partner!=-2) {
                isdone[isquark[partner]]=1;
                two_end[istring]=isquark[partner];
                istring+=1;
            }
            else {
                pIn.push_back ( Parton ( rempx, rempy, rempz, reme, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true ) );
                isquark[nquarks]=pIn.size()-1;
                nquarks+=1;
                isdone[pIn.size()-1]=1;
                two_end[istring]=pIn.size()-1;
                cout << "Attached quark remnant flying down +Pz beam" << endl;
            }
        }
    }

    // Assign gluons to a certain string
    int my_string[pIn.size()];
    memset( my_string, 0, pIn.size()*sizeof(int) );
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
        if (pIn[ipart].GetId()==21) {
            double min_delR=100000.;
            for (int ns=0; ns<nstrings; ns++)
            {         
                int fq=one_end[ns];
                int sq=two_end[ns];
          	double f_delR = pIn[ipart].delta_R(pIn[fq]);
		double s_delR = pIn[ipart].delta_R(pIn[sq]);
          	double delR=(f_delR+s_delR)/2.;
               if (delR<min_delR) my_string[ipart]=ns, min_delR=delR;
            }
        }
    }

    // Build up chain using gluons assigned to each string, in a closest pair order
    int lab_col=102;
    for (int ns=0; ns<nstrings; ns++)
    {
        event.reset();

	//cout << " NEXT STRING " << endl;
        // Energy weighted starting position of the hadronization process
	double string_pos[4]={0.};
	double sum_pt=0.;
    
        // First end
        int tquark=one_end[ns];
        if (pIn[tquark].GetId()>0) col[tquark]=lab_col;
        else acol[tquark]=lab_col;
        lab_col+=1;
        int link=tquark;

	// Position of parton
	vector<double> parton_pos=pIn[tquark].hyper_point();
	if (parton_pos[3]==-1.) parton_pos=pIn[tquark].vGetPos();
	else {
          double temp_z=parton_pos[3]*sinh(parton_pos[2]);
	  double temp_t=parton_pos[3]*cosh(parton_pos[2]);
	  parton_pos[2]=temp_z;
	  parton_pos[3]=temp_t;
        }
	if (parton_pos[3]!=-1.) { 
          string_pos[0]+=parton_pos[0]*pIn[tquark].GetPt();
          string_pos[1]+=parton_pos[1]*pIn[tquark].GetPt();
          string_pos[2]+=parton_pos[2]*pIn[tquark].GetPt();
          string_pos[3]+=parton_pos[3]*pIn[tquark].GetPt();
	  //cout << " parton_pos= " << parton_pos[0] << " " <<  parton_pos[1] << " " << parton_pos[2] << " " << parton_pos[3] << " " << endl;
          sum_pt+=pIn[tquark].GetPt();
        }
        else {
          cout << " First quark Unassigned position! id= " << pIn[tquark].GetOrig() << endl;
	}

	// Feed into PYTHIA
        int ide=pIn[tquark].GetId();
        double px=pIn[tquark].vGetP()[0];
        double py=pIn[tquark].vGetP()[1];
        double pz=pIn[tquark].vGetP()[2];
        double ee=pIn[tquark].vGetP()[3];
        double mm=pdt.m0(int(ide));
        ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
        if (col[tquark]==0 && acol[tquark]==0 && (ide==21 || abs(ide)<=6)) {
            cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
            exit(0);
        }
        event.append(int(ide),23,col[tquark],acol[tquark],px,py,pz,ee,mm);
        
	int changes=1;
        do {
            changes=0;
            double min_delR=100000.;
            int next_link=0;
            for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
            {
                if (pIn[ipart].GetId()==21 && isdone[ipart]==0 && my_string[ipart]==ns)
                {
                    changes=1;
                    double delR = pIn[link].delta_R(pIn[ipart]);
                    if (delR<min_delR) min_delR=delR, next_link=ipart;
                }
            }
            if (changes==1)
            {
		// Position of parton
	        vector<double> parton_pos=pIn[next_link].hyper_point();
	        if (parton_pos[3]==-1.) parton_pos=pIn[next_link].vGetPos();
	        else {
                  double temp_z=parton_pos[3]*sinh(parton_pos[2]);
	          double temp_t=parton_pos[3]*cosh(parton_pos[2]);
	          parton_pos[2]=temp_z;
	          parton_pos[3]=temp_t;
                }
	        if (parton_pos[3]!=-1.) { 
                  string_pos[0]+=parton_pos[0]*pIn[next_link].GetPt();
                  string_pos[1]+=parton_pos[1]*pIn[next_link].GetPt();
                  string_pos[2]+=parton_pos[2]*pIn[next_link].GetPt();
                  string_pos[3]+=parton_pos[3]*pIn[next_link].GetPt();
                  sum_pt+=pIn[next_link].GetPt();
	          //cout << " parton_pos= " << parton_pos[0] << " " <<  parton_pos[1] << " " << parton_pos[2] << " " << parton_pos[3] << " " << endl;
                }
                else {
                  cout << " Gluon Unassigned position! id= " << pIn[next_link].GetOrig() << endl;
	        }
	       
                isdone[next_link]=1;
                if (col[link]==lab_col-1) col[next_link]=lab_col, acol[next_link]=lab_col-1;
                else col[next_link]=lab_col-1, acol[next_link]=lab_col;
                lab_col+=1;
                link=next_link;
		
		// Feed into PYTHIA
                ide=pIn[next_link].GetId();
                px=pIn[next_link].vGetP()[0];
                py=pIn[next_link].vGetP()[1];
                pz=pIn[next_link].vGetP()[2];
                ee=pIn[next_link].vGetP()[3];
                mm=pdt.m0(int(ide));
                ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
                if (col[next_link]==0 && acol[next_link]==0 && (ide==21 || abs(ide)<=6)) {
                  cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
                  exit(0);
                }
                event.append(int(ide),23,col[next_link],acol[next_link],px,py,pz,ee,mm);
	    }
        } while (changes==1);
        // Attach second end
        if (col[link]==lab_col-1) col[two_end[ns]]=0, acol[two_end[ns]]=lab_col-1;
        else col[two_end[ns]]=lab_col-1, acol[two_end[ns]]=0;
	        
	parton_pos=pIn[two_end[ns]].hyper_point();
	if (parton_pos[3]==-1.) parton_pos=pIn[two_end[ns]].vGetPos();
	else {
	  double temp_z=parton_pos[3]*sinh(parton_pos[2]);
	  double temp_t=parton_pos[3]*cosh(parton_pos[2]);
          parton_pos[2]=temp_z;
          parton_pos[3]=temp_t;
	}
        if (parton_pos[3]!=-1.) {
	  string_pos[0]+=parton_pos[0]*pIn[two_end[ns]].GetPt();
	  string_pos[1]+=parton_pos[1]*pIn[two_end[ns]].GetPt();
	  string_pos[2]+=parton_pos[2]*pIn[two_end[ns]].GetPt();
	  string_pos[3]+=parton_pos[3]*pIn[two_end[ns]].GetPt();
	  sum_pt+=pIn[two_end[ns]].GetPt();
	  //cout << " parton_pos= " << parton_pos[0] << " " <<  parton_pos[1] << " " << parton_pos[2] << " " << parton_pos[3] << " " << endl;
	}
        else {
	  cout << " Second quark Unassigned position! id= " << pIn[two_end[ns]].GetOrig() << endl;
        }
		
        if (col[two_end[ns]]!=0) { 
            if (pIn[two_end[ns]].GetId()<0) pIn[two_end[ns]].SetId(-pIn[two_end[ns]].GetId());
        }
        else {
            if (pIn[two_end[ns]].GetId()>0) pIn[two_end[ns]].SetId(-pIn[two_end[ns]].GetId());
        }

        // Feed into PYTHIA
        ide=pIn[two_end[ns]].GetId();
        px=pIn[two_end[ns]].vGetP()[0];
        py=pIn[two_end[ns]].vGetP()[1];
        pz=pIn[two_end[ns]].vGetP()[2];
        ee=pIn[two_end[ns]].vGetP()[3];
        mm=pdt.m0(int(ide));
        ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
        if (col[two_end[ns]]==0 && acol[two_end[ns]]==0 && (ide==21 || abs(ide)<=6)) {
          cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
          exit(0);
        }
        event.append(int(ide),23,col[two_end[ns]],acol[two_end[ns]],px,py,pz,ee,mm);
        
	for (unsigned a=0; a<4; a++) {
          string_pos[a]/=sum_pt;
          //cout << " String Pos= " << string_pos[a] << endl;  
        }
    
        if (!hpythia.next()) continue;
        for (unsigned int ipart=0; ipart < event.size(); ++ipart)
        {
          if (event[ipart].isFinal())
          {
              int ide=event[ipart].id();
              vector<double> p {event[ipart].px(),event[ipart].py(),event[ipart].pz(),event[ipart].e()};
              hadron_list.push_back ( Parton ( p, 0., event[ipart].m(), 0, -1, -1, ide, "hadron", 0, 0, true ) );
    
              double had_pos[4];
	      had_pos[0] = string_pos[0] + event[ipart].xProd()*MM2FM;	    
	      had_pos[1] = string_pos[1] + event[ipart].yProd()*MM2FM;	    
	      had_pos[2] = string_pos[2] + event[ipart].zProd()*MM2FM;	    
	      had_pos[3] = string_pos[3] + event[ipart].tProd()*MM2FM;
	      //cout << " xprod= " << event[ipart].xProd()*MM2FM
		//   << " yprod= " << event[ipart].yProd()*MM2FM
		//   << " zprod= " << event[ipart].zProd()*MM2FM
		//   << " tprod= " << event[ipart].tProd()*MM2FM << endl;

              hadron_list[hadron_list.size()-1].SetPos(had_pos[0],had_pos[1],had_pos[2],had_pos[3]);

	      //Print on output file
              hadfile << event[ipart].px() << " " << event[ipart].py() << " " << event[ipart].pz() << " " << event[ipart].e() << " " 
		    << had_pos[0] << " " << had_pos[1] << " " << had_pos[2] << " " << had_pos[3] << " "
		    << event[ipart].id() << endl;
          }  
        }

    }

}

