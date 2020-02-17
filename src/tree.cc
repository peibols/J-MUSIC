#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Parton.h"
#include "jets.h"
#include <vector>

using std::vector;
using namespace std;
using namespace Pythia8;

void Jets::TreeInit()
{
	//Read cmnd file
	pythia.readFile("setup_pythia.cmnd");

	//Initialize random number generator
	randi.init(DATA.Jets_seed);
	
 	//Set Random Seed
        ostringstream seedstring;
        seedstring << "Random:seed = " << DATA.Jets_seed; 
        pythia.readString(seedstring.str().c_str());	
	
	//Initialize PYTHIA
	cout << " Initialize PYTHIA " << endl;
	pythia.init();
}

bool Jets::TreeDoer(vector<double> pos_coll, double total_cross, int icoll)
{
	if (!pythia.next()) return 0;
	//cout << " pythia.event.size= " << pythia.event.size() << endl;
        //cout << "cross section= " << scientific << setprecision(8) << pythia.info.sigmaGen() << endl;
        if (DATA.one_hard_collision==0) {
	  double cross_ratio=pythia.info.sigmaGen()/total_cross;
	  if (cross_ratio>1.) cout << " WTF CROSS RATIO > 1 " << endl;
	  if (randi.flat()>cross_ratio) return 0;
	}
	double weight=pythia.info.weight();
	//cout << "Hard Weight = " << setprecision(6) << weight << endl;
	//cout << "Hard Cross = " << setprecision(6) << pythia.info.sigmaGen() << endl;
	event_weight=weight;
	event_cross=pythia.info.sigmaGen();

	vector<Parton> partons=parton_list;

	//Add wrt current number of partons
	int current_npartons = partons.size();

	//Find Final Particles, excluding remnants
	for (int i = 0; i < pythia.event.size(); i++)
	{
		if (pythia.event[i].isFinal())
		{
			//Simply store remnants
			if (pythia.event[i].status() == 63)
			{
				partons.push_back ( Parton ( pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e(), 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "rem", pythia.event[i].col(), pythia.event[i].acol(), true ) );
				partons[partons.size()-1].vSetPos(pos_coll);
				partons[partons.size()-1].vSetPosIn(pos_coll);
				continue;
			}
			
			//Find first non-trivial mother
			int use = i;
			int m1 = 0;
			int m2 = 0;
			do
			{
				m1 = pythia.event[use].mother1();
				m2 = pythia.event[use].mother2();
				if (m1==m2) use = m1;
			} while (m1==m2);
			
			//Compute virtuality
			double virt=abs(sqrt(pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2()));
			
			//Add it to partons array
			partons.push_back ( Parton ( pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e(), virt, pythia.event[i].m(), m1, -1, -1, pythia.event[i].id(), "", pythia.event[i].col(), pythia.event[i].acol(), false ) );	

			//If mother is ISR initiator, say mother is -1
			if (pythia.event[m1].status()==-41)
			{
				partons[partons.size()-1].SetMom(-1);
				partons[partons.size()-1].SetOrig("isr");
				partons[partons.size()-1].SetIsDone(true);
				partons[partons.size()-1].vSetPos(pos_coll);
				partons[partons.size()-1].vSetPosIn(pos_coll);
			}

			//If mother is hard scattering initiator, say mother is -1
			if (pythia.event[m1].status()==-21)
			{
				//cout << " ahora " << endl;
				partons[partons.size()-1].SetMom(-1);
				partons[partons.size()-1].SetOrig("hs");
				partons[partons.size()-1].SetIsDone(true);
				partons[partons.size()-1].vSetPos(pos_coll);
				partons[partons.size()-1].vSetPosIn(pos_coll);
			}

		}
	}
		
	/*
	cout << "#Partons in array = " << partons.size() << endl;
	for (unsigned int i = 0; i < partons.size(); i++)
	{
		cout << " Parton " << i << " has mom= " << partons[i].GetMom() << " is done or not= " << partons[i].GetIsDone() << endl;
	}
	*/
	
	//Reconstruct tree, using momentum of daughters to get the one for mothers
	int changes=0;
	do
	{
		changes=1;
		unsigned int ps = partons.size();
		for (unsigned int i = current_npartons; i < ps; i++)
		{
			if (partons[i].GetIsDone()==false)
			{
				for (unsigned int j = current_npartons; j < ps; j++)
				{
					if (partons[i].GetMom()==partons[j].GetMom() && i!=j && partons[j].GetIsDone()==false)
					{
						int mom=partons[i].GetMom();

						//Set mother momentum
						double px=partons[i].vGetP()[0]+partons[j].vGetP()[0];
						double py=partons[i].vGetP()[1]+partons[j].vGetP()[1];
						double pz=partons[i].vGetP()[2]+partons[j].vGetP()[2];
						double en=partons[i].vGetP()[3]+partons[j].vGetP()[3];
						double virt=abs(sqrt(pow(en,2.)-pow(px,2.)-pow(py,2.)-pow(pz,2.)-pow(pythia.event[mom].m(),2.)));

						//Find first non-trivial mother of the mother
						int use = mom;
						int m1 = 0;
						int m2 = 0;
						do {
							m1 = pythia.event[use].mother1();
							m2 = pythia.event[use].mother2();
							if (m1==m2) use = m1;
						} while (m1==m2);

						//Fill it in partons array
						partons.push_back ( Parton ( px, py, pz, en, virt, pythia.event[mom].m(), m1, i, j, pythia.event[mom].id(), "", pythia.event[mom].col(), pythia.event[mom].acol(), false ) );
						//Update mother of daughters to point to position in partons array instead of pythia list, and declare as done
						partons[i].SetMom(partons.size()-1);
						partons[j].SetMom(partons.size()-1);
						partons[i].SetIsDone(true);
						partons[j].SetIsDone(true);

						//If mother is ISR initiator, say mother of mother is -1
						if (pythia.event[m1].status()==-41)
						{
							partons[partons.size()-1].SetMom(-1);
							partons[partons.size()-1].SetOrig("isr");
							partons[partons.size()-1].SetIsDone(true);
							partons[partons.size()-1].vSetPos(pos_coll);
						}

						//If mother is hard scattering initiator, say mother of mother is -1
						if (pythia.event[m1].status()==-21)
						{
							partons[partons.size()-1].SetMom(-1);
							partons[partons.size()-1].SetOrig("hs");
							partons[partons.size()-1].SetIsDone(true);
							partons[partons.size()-1].vSetPos(pos_coll);
						}
							//cout << " Hola " << partons.size() << " daugh1= " << i << " daugh2= " << j << " mom= " << mom << " array mom= " << partons[i].GetMom() << endl;
						changes=0;
						break;
					}
				}
			}
		}
	} while (changes==0);

	parton_list=partons;

/*	
	//Compute partons 4-momentum
	vector<double> pfour (4,0.);
	for (unsigned int iparton=current_npartons; iparton<partons.size(); iparton++) {
	    if (partons[iparton].GetD1()==-1 && partons[iparton].GetOrig()!="rem") {
                for (unsigned a=0; a<4; a++) {
	            pfour[a]+=partons[iparton].vGetP()[a];
	        }  
	    }
	}
	//Assign it to binary collision for hydro
	binary_list[icoll]->dpxdtau=-pfour[0];
	binary_list[icoll]->dpydtau=-pfour[1];
	binary_list[icoll]->dpzdtau=-pfour[2];
	binary_list[icoll]->dEdtau=-pfour[3];
*/

	//cout << " En= " << pfour[3] << " Px= " << pfour[0] << " Py= " << pfour[1] << " Pz= " << pfour[2] << endl;

	//cout << " I aki partons size = " << partons.size() << endl;
	/*
	for (unsigned int i = 0; i < partons.size(); i++)
	{
		cout << " Parton " << i << " has mom= " << partons[i].GetMom() << " is done or not " << partons[i].GetIsDone() << " orig= " << partons[i].GetOrig();
		cout << " and has daughters = " << partons[i].GetD1() << " and  " << partons[i].GetD2() << endl;
	}
	*/
	/*
	//Write to Output File
	for (unsigned int i = 0; i < partons.size(); i++)
	{
		geneal << partons[i].GetPx() << " " << partons[i].GetPy() << " " << partons[i].GetPz() << " " << partons[i].GetEn() << " " << partons[i].GetQ() << " ";
		geneal << partons[i].GetMom() << " " << partons[i].GetPx() << " "   	
	}
	*/	

	return 1;
}

