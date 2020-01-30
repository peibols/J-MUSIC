#ifndef SRC_PARTON_H
#define SRC_PARTON_H

#include <vector>
#include <string>
#include <cmath>

using std::vector;

class Parton
{
	protected:
		vector<double> _p;	//Current momentum
		vector<double> _x;	//Current position

		vector<double> _pi;	//Initial, unquenched momentum
		vector<double> _pinh;	//Inherited momentum at creation
		vector<double> _xi;	//Position at creation

		double _lambda;

		double _f_dist;

		double _tau_form;
		double _splitTime;

		double _q;
		double _mass;
	
		int _mom;
		int _d1;
		int _d2;

		int _id;
		std::string _orig;

		int _col;
		int _acol;

		bool _isdone;

	public:
		Parton();
		Parton(vector<double> p, double q, double mass, int mom, int d1, int d2, int id, std::string orig, int col, int acol, bool isdone);
		Parton(double px, double py, double pz, double en, double q, double mass, int mom, int d1, int d2, int id, std::string orig, int col, int acol, bool isdone);
		~Parton();

		virtual void display() const;
		
		double delta_R(Parton other_parton);

		void reset_momentum(double px, double py, double pz);

		double GetFluidDist() const {return _f_dist;}
		void AddFluidDist(double delta_f_dist) {_f_dist+=delta_f_dist;}

		double GetTau() const;
		double GetsEta() const;

		void SetsplitTime();
		double splitTime() {return _splitTime;}

		void SetLambda(double lambda) { _lambda = lambda; }
		double GetLambda() {return _lambda;}

		double GetTauForm();

		void vSetPosIn(vector<double> x);
		vector<double> vGetPosIn() const;
		
		void vSetPosHy(vector<double> x);
		vector<double> vGetPosHy() const;

		void vSetPos(vector<double> x);
		void SetPos(double x, double y, double z, double t);
		vector<double> vGetPos() const;
	
		void vSetP(vector<double> p);
		void SetP(double px, double py, double pz, double en);
		vector<double> vGetP() const;

		void vSetPIn(vector<double> p);
		vector<double> vGetPIn() const;
                
		void vSetPinh(vector<double> p);
		vector<double> vGetPinh() const;
		
		double GetPt() const;

		double GetEta() const;

		void SetQ(double q);
                double GetQ() const;

		void SetMass(double mass);
                double GetMass() const;

		void SetMom(int mom);
                int GetMom() const;

		void SetD1(int d1);
                int GetD1() const;

		void SetD2(int d2);
                int GetD2() const;

		void SetId(int id);
                int GetId() const;

		void SetOrig(std::string orig);
                std::string GetOrig() const;

		void SetCol(int col);
                int GetCol() const;

		void SetAcol(int acol);
                int GetAcol() const;

		void SetIsDone(bool isdone);
                bool GetIsDone() const;
};

#endif // SRC_PARTON_H
