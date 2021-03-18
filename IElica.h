#ifndef _IELICA
#define _IELICA

#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"
#include <vector>
#include <iostream>

using namespace std;

class IElica{
  
 private:
  double m_x0;
  double m_y0;
  double m_z0;
  double m_R;
  double m_phi;
  int m_h;
  double m_lambda;
  
 public:
  IElica();
  IElica(double x0, double y0, double z0, double R, double phi, int h, double lambda);

  void SetParameters(double x0, double y0, double z0, double R, double phi, int h, double lambda);
  
  void SetX0(double x0);
  void SetY0(double y0);
  void SetZ0(double z0);
  double GetX0();
  double GetY0();
  double GetZ0();
  
  void SetRadius(double R);
  double GetRadius();

  void SetPhi(double phi);
  double GetPhi();

  void SetH(int h);
  int GetH();

  void SetLambda(double lambda);
  double GetLambda();

  TVector3 GetTangentVector(double s);
  void RandomizeFromTangent(TVector3 tan, double deviation, TRandom3 &rnd);

  vector<double> Eval(double s);

  void PrintParameters();

  void GetXIntersection(double x, vector<double> &pos);
  
};

#endif
