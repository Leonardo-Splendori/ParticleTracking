#ifndef _IRILEVATORE
#define _IRILEVATORE

#include <iostream>
#include <vector>
#include "IElica.h"

#endif

using namespace std;

class IRilevatore{

 private:
  double m_position;
  double m_width;
  double m_height;
  double m_depth;
  double m_resolution;
  vector <vector<double>> m_hits;
  
 public:
  IRilevatore();
  IRilevatore(double position, double width, double height, double depth, double resolution);

  void SetHits(vector<vector<double>> &hits);
  void SetHit(double y, double z);
  void SetPosition(double position);
  void SetWidth(double width);
  void SetHeight(double height);
  void SetDepth(double depth);
  void SetResolution(double resolution);
  
  void GetHits(vector<vector<double>> &hits);
  double GetHitY(unsigned int i);
  double GetHitZ(unsigned int i);
  double GetPosition();
  double GetWidth();
  double GetHeight();
  double GetDepth();
  double GetResolution();
  int GetNumberOfHits();
  void Reset();
  
  void GenerateHits(double smin, double smax, IElica &Elica, double step);
};
