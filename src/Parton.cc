#include "Parton.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using std::vector;

Parton::Parton()
{

}

Parton::Parton(double px, double py, double pz, double en, double q, double mass, int mom, int d1, int d2, int id, std::string orig, int col, int acol, bool isdone)
{
	_p.push_back(px);
	_p.push_back(py);
	_p.push_back(pz);
	_p.push_back(en);

	_pi = _p;
	_pinh = _p;

	_f_dist=0.;

	for (unsigned a=0; a<4; a++) _x.push_back(-1.);
	for (unsigned a=0; a<4; a++) _xi.push_back(-1.);

	if (d1!=-1) {
	  if (q<0.01) q=0.01;
	  _tau_form=0.2*2.*en/pow(q,2.);
	}  
	else _tau_form=-1.;

	_lambda=1.;
	_splitTime=-1.;

	_q=q;
	_mass=mass;
	
	_mom=mom;
	_d1=d1;
	_d2=d2;

	_id=id;
	_orig=orig;

	_col=col;
	_acol=acol;

	_isdone=isdone;
}

Parton::Parton(vector<double> p, double q, double mass, int mom, int d1, int d2, int id, std::string orig, int col, int acol, bool isdone)
{
	_p=p;
	
	_pi = _p;
	_pinh = _p;

	_f_dist=0.;

	if (d1!=-1) {
	  if (q<0.01) q=0.01;
	  _tau_form=0.2*2.*p[3]/pow(q,2.);
	}  
	else _tau_form=-1.;

	for (unsigned a=0; a<4; a++) _x.push_back(-1.);
	for (unsigned a=0; a<4; a++) _xi.push_back(-1.);

	_lambda=1.;
	_splitTime=-1.;

	_q=q;
	_mass=mass;
	
	_mom=mom;
	_d1=d1;
	_d2=d2;

	_id=id;
	_orig=orig;

	_col=col;
	_acol=acol;

	_isdone=isdone;
}

Parton::~Parton()
{
	//std::cout << "Parton destructor called" << std::endl;
}

double Parton::delta_R(Parton other_parton)
{
    vector<double> op=other_parton.vGetP();
    double delta_eta=this->GetEta()-other_parton.GetEta();
    double delta_phi=acos((_p[0]*other_parton.vGetP()[0]+_p[1]*other_parton.vGetP()[1])/this->GetPt()/other_parton.GetPt());
    return sqrt(delta_phi*delta_phi+delta_eta*delta_eta);
}

double Parton::GetTau() const
{
    double tau=sqrt(_x[3]*_x[3]-_x[2]*_x[2]);
    if (tau!=tau) tau=0.;
    return tau;
}

void Parton::display() const
{
	std::cout << " Px = " << _p[0] << " Py= " << _p[1] << " Pz = " << _p[2] << " En= " << _p[3] << std::endl;
	std::cout << " Lambda= " << _lambda << " Px_in= " << _pi[0] << " Py_in= " << _pi[1] << " Pz_in= " << _pi[2] << " En_in= " << _pi[3] << std::endl;
	std::cout << " time= " << _x[3] << " mom= " << _mom << " d1= " << _d1 << " d2= " << _d2 << " id= " << _id << " orig= " << _orig << " IsDone?= " << _isdone << std::endl;
}

double Parton::GetTauForm()
{
    return _tau_form;
}

void Parton::SetsplitTime()	//Splitting time in tau coordinates
{
  double gamma=cosh(this->GetsEta());
  _splitTime=(_xi[3]+_tau_form)/gamma;
}

void Parton::reset_momentum(double px, double py, double pz)
{
    double en=sqrt(px*px+py*py+pz*pz+_mass*_mass);
    _p[0]=px;
    _p[1]=py;
    _p[2]=pz;
    _p[3]=en;
}

void Parton::vSetPos(vector<double> x)
{
      _x=x;
}

void Parton::vSetPosIn(vector<double> x)
{
      _xi=x;
}

vector<double> Parton::vGetPos() const
{
    return _x;
}

vector<double> Parton::vGetPosIn() const
{
    return _xi;
}

void Parton::vSetP(vector<double> p)
{
      _p=p;
}
void Parton::SetP(double px, double py, double pz, double en)
{
	_p[0]=px;
	_p[1]=py;
	_p[2]=pz;
	_p[3]=en;
}
vector<double> Parton::vGetP() const
{
	return _p;
}

void Parton::vSetPinh(vector<double> p)
{
      _pinh=p;
}
vector<double> Parton::vGetPinh() const
{
	return _pinh;
}

void Parton::vSetPIn(vector<double> p)
{
      _pi=p;
}
vector<double> Parton::vGetPIn() const
{
	return _pi;
}

double Parton::GetPt() const
{
	double pt=sqrt(pow(_p[0],2.)+pow(_p[1],2.));
	return pt;
}

double Parton::GetEta() const  //Pseudorapidity
{
	double pt=sqrt(pow(_p[0],2.)+pow(_p[1],2.));
	double eta=1./2.*log((sqrt(pow(pt,2.)+pow(_p[2],2.))+_p[2])/(sqrt(pow(pt,2.)+pow(_p[2],2.))-_p[2]));
	return eta;
}

double Parton::GetsEta() const  //Spatial Rapidity
{
  	double seta;
  	if (_x[2]!=0.) {	//Use actual spatial rapidity
		seta=1./2.*std::log((_x[3]+_x[2])/(_x[3]-_x[2]));
	}
	if (_x[2]==0.) {	//Use momentum rapidity for first step
		double modp=sqrt(_p[0]*_p[0]+_p[1]*_p[1]+_p[2]*_p[2]);
  		seta=1./2.*log((modp+_p[2])/(modp-_p[2]));
	}
	if (_x[3]<=fabs(_x[2]) && _x[3]!=0.) {	//Return large values for ill defined rapidities	
	    	std::cout << " PROBLEM, VERY HIGH RAPIDITY _x[2]= " << _x[2] << " _x[3]= " << _x[3] << std::endl;
  	    if (_x[2]>0.) seta=9999999.;
	    if (_x[2]<0.) seta=-9999999.;
        }
	
	return seta;
}

void Parton::SetQ(double q)
{
        _q=q;
}
double Parton::GetQ() const
{
        return _q;
}

void Parton::SetMass(double mass)
{
        _mass=mass;
}
double Parton::GetMass() const
{
        return _mass;
}

void Parton::SetMom(int mom)
{
        _mom=mom;
}
int Parton::GetMom() const
{
        return _mom;
}

void Parton::SetD1(int d1)
{
        _d1=d1;
}
int Parton::GetD1() const
{
        return _d1;
}

void Parton::SetD2(int d2)
{
        _d2=d2;
}
int Parton::GetD2() const
{
        return _d2;
}

void Parton::SetId(int id)
{
        _id=id;
}
int Parton::GetId() const
{
        return _id;
}

void Parton::SetOrig(std::string orig)
{
        _orig=orig;
}
std::string Parton::GetOrig() const
{
        return _orig;
}

void Parton::SetCol(int col)
{
        _col=col;
}
int Parton::GetCol() const
{
        return _col;
}

void Parton::SetAcol(int acol)
{
        _acol=acol;
}
int Parton::GetAcol() const
{
        return _acol;
}

void Parton::SetIsDone(bool isdone)
{
        _isdone=isdone;
}
bool Parton::GetIsDone() const
{
        return _isdone;
}
