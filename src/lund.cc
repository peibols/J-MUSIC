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

    // Don't let pi0 decay
    //hpythia.readString("111:mayDecay = off");
    // Don't let any hadron decay
    hpythia.readString("HadronLevel:Decay = off");

    // And initialize
    hpythia.init();

    hadron_list.clear();
}

void Jets::HadronizeTherm() {
    
    Event& event      = hpythia.event;
    ParticleData& pdt = hpythia.particleData;

    event.reset();

    vector<Parton> pIn = fin_and_therm_parton_list;
    
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
    
    double fin_pos[4];
    fin_pos[0] = hyper_point[0];
    fin_pos[1] = hyper_point[1];
    fin_pos[2] = hyper_point[3]*sinh(hyper_point[2]);
    fin_pos[3] = hyper_point[3]*cosh(hyper_point[2]);

    hpythia.next();
    hadfile <<"weight " << event_weight << " cross " << event_cross << endl;
    
    for (unsigned int ipart=0; ipart < event.size(); ++ipart)
    {
        if (event[ipart].isFinal())
        {
            int ide=hpythia.event[ipart].id();
            vector<double> p {hpythia.event[ipart].px(),hpythia.event[ipart].py(),hpythia.event[ipart].pz(),hpythia.event[ipart].e()};
            hadron_list.push_back ( Parton ( p, 0., hpythia.event[ipart].m(), 0, -1, -1, ide, "hadron", 0, 0, true ) );
    
            double had_pos[4];
	    had_pos[0] = fin_pos[0] + pythia.event[i].xProd()*MM2FM;	    
	    had_pos[1] = fin_pos[1] + pythia.event[i].yProd()*MM2FM;	    
	    had_pos[2] = fin_pos[2] + pythia.event[i].zProd()*MM2FM;	    
	    had_pos[3] = fin_pos[3] + pythia.event[i].tProd()*MM2FM;

            hadron_list[hadron_list.size()-1].SetPos(had_pos[0],had_pos[1],had_pos[2],had_pos[3]);

	    //Print on output file
            hadfile << hpythia.event[ipart].px() << " " << hpythia.event[ipart].py() << " " << hpythia.event[ipart].pz() << " " << hpythia.event[ipart].e() << " " 
		    << had_pos[0] << " " << had_pos[1] << " " << had_pos[2] << " " << had_pos[3] << " "
		    << hpythia.event[ipart].id() << endl;
        }
    }

    hadfile.close();   
}

void Jets::HadronizeJets() {
    
    Event& event      = hpythia.event;
    ParticleData& pdt = hpythia.particleData;

    event.reset();

    hadron_list.clear();

    double Lambda_QCD=0.2;
    double rempx=0.2;
    double rempy=0.2;
    double p_fake=1380.;
    double rempz=p_fake;
    double reme=std::sqrt(std::pow(rempx,2.)+std::pow(rempy,2.)+std::pow(rempz,2.));
    
    // Hadronize all showers together
    vector<Parton> pIn = final_parton_list;

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
        int tquark=one_end[ns];
        if (pIn[tquark].GetId()>0) col[tquark]=lab_col;
        else acol[tquark]=lab_col;
        lab_col+=1;
        int link=tquark;
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
                isdone[next_link]=1;
                if (col[link]==lab_col-1) col[next_link]=lab_col, acol[next_link]=lab_col-1;
                else col[next_link]=lab_col-1, acol[next_link]=lab_col;
                lab_col+=1;
                link=next_link;
            }
        } while (changes==1);
        // Attach second end
        if (col[link]==lab_col-1) col[two_end[ns]]=0, acol[two_end[ns]]=lab_col-1;
        else col[two_end[ns]]=lab_col-1, acol[two_end[ns]]=0;
    }
    // Changing identity of quarks to be consistent with color charge
    for( int iq=0; iq <  nquarks; ++iq)
    {
        if (col[isquark[iq]]!=0) { 
            if (pIn[isquark[iq]].GetId()<0) pIn[isquark[iq]].SetId(-pIn[isquark[iq]].GetId());
        }
        else {
            if (pIn[isquark[iq]].GetId()>0) pIn[isquark[iq]].SetId(-pIn[isquark[iq]].GetId());
        }
    }

    // Introduce partons into PYTHIA
    /*
    for (unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
      JSDEBUG << "Parton #" << ipart << " is a " << pIn[ipart]->pid() << "with energy = " << pIn[ipart]->e() << " with phi= " << pIn[ipart]->phi() << " and has col= " << col[ipart] << " and acol= " << acol[ipart];
    }
    */
      
/***************************************************************************************************************/
//
// Making collinear partons not collinear
//
// Could have dangerous effects, yet to be tested....
//
            
       for (unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
       {
           double px=pIn[ipart].vGetP()[0];
           double py=pIn[ipart].vGetP()[1];
           double pz=pIn[ipart].vGetP()[2];
           double ee=pIn[ipart].vGetP()[3];
           for (unsigned int j = ipart+1; j< pIn.size(); j++)
           {
               double p2x = pIn[j].vGetP()[0];
               double p2y = pIn[j].vGetP()[1];
               double p2z = pIn[j].vGetP()[2];
               double e2e = pIn[j].vGetP()[3];
               
               double diff = sqrt( pow(px-p2x,2) + pow(py-p2y,2) + pow(pz-p2z,2) );
               double f = 4.0;
               if (diff < f*Lambda_QCD)
               {
                   if ((pz>=0)&&(p2z>=0))
                   {
                     if (pz>=p2z)
                     {
                         pIn[ipart].SetP(px,py,pz+f*Lambda_QCD,sqrt(ee*ee + 2*pz*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD) );
                     }
                     else pIn[j].SetP(p2x,p2y,p2z+f*Lambda_QCD,sqrt(e2e*e2e + 2*p2z*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD));
                   }
                   else if ((pz>=0)&&(p2z<0))
                   {
                       if (abs(pz)>=abs(p2z))
                       {
                           pIn[ipart].SetP(px,py,pz+f*Lambda_QCD,sqrt(ee*ee + 2*pz*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD) );
                       }
                       else pIn[j].SetP(p2x,p2y,p2z-f*Lambda_QCD,sqrt(e2e*e2e - 2*p2z*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD));
                   }
                   else if ((pz<0)&&(p2z>=0))
                   {
                       if (abs(pz)>=abs(p2z))
                       {
                           pIn[ipart].SetP(px,py,pz-f*Lambda_QCD,sqrt(ee*ee - 2*pz*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD));
                       }
                       else pIn[j].SetP(p2x,p2y,p2z+f*Lambda_QCD,sqrt(e2e*e2e + 2*p2z*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD));
                   }
                   else
                   {
                       if (abs(pz)>=abs(p2z))
                       {
                           pIn[ipart].SetP(px,py,pz-f*Lambda_QCD,sqrt(ee*ee - 2*pz*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD));
                       }
                       else pIn[j].SetP(p2x,p2y,p2z-f*Lambda_QCD,sqrt(e2e*e2e - 2*p2z*f*Lambda_QCD + f*Lambda_QCD*f*Lambda_QCD));
                   }
               }
           }
       }
      
      
/**************************************************************************************************************************/

    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {  
        int ide=pIn[ipart].GetId();
        double px=pIn[ipart].vGetP()[0];
        double py=pIn[ipart].vGetP()[1];
        double pz=pIn[ipart].vGetP()[2];
        double ee=pIn[ipart].vGetP()[3];
        double mm=pdt.m0(int(ide));
        ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
        if (col[ipart]==0 && acol[ipart]==0 && (ide==21 || abs(ide)<=6)) {
            cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
            exit(0);
        }
        event.append(int(ide),23,col[ipart],acol[ipart],px,py,pz,ee,mm);
    }
    
    hpythia.next();
    hadfile <<"weight " << event_weight << " cross " << event_cross << endl;
    for (unsigned int ipart=0; ipart < event.size(); ++ipart)
    {
        if (event[ipart].isFinal())
        {
            int ide=hpythia.event[ipart].id();
            vector<double> p {hpythia.event[ipart].px(),hpythia.event[ipart].py(),hpythia.event[ipart].pz(),hpythia.event[ipart].e()};
            hadron_list.push_back ( Parton ( p, 0., hpythia.event[ipart].m(), 0, -1, -1, ide, "hadron", 0, 0, true ) );
            //Print on output file
            hadfile << hpythia.event[ipart].px() << " " << hpythia.event[ipart].py() << " " << hpythia.event[ipart].pz() << " " << hpythia.event[ipart].e() << " " << hpythia.event[ipart].id() << endl;
        }
    }

    hadfile.close();   
}

