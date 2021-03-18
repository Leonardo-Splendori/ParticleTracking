#ifndef _RIEMANNFITTER
#define _RIEMANNFITTER

#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDEigen.h"
#include <iostream>

using namespace std;

class IRiemannFitter{
 private:
  vector<TVector2> m_points;
  vector<TVector3> m_points_sphere;
  double m_scaling_factor;
  void MapOnSphere();
  double GetPointP(unsigned int i);
  double GetCovariance(int i, int j);
  TVector3 GetPointsAverages();
 public:
  void SetPoints(vector<TVector2> &pos);
  void PrintPointsCartesian();
  void PrintPointsPolar();
  void PrintPointsSphere();
  TMatrixD GetCovarianceMatrix();
  void ScaleAndCenter(double b);
  void Reset();
  TVector3 Fit();
};





#endif
