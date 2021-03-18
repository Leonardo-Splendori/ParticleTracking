#ifndef _IPATTERNFINDER
#define _IPATTERNFINDER

#include "TH1D.h"
#include "TVector2.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <vector>
#include <iostream>

using namespace std;

class IPatternFinder : public TH1D {
 private:
  int m_ntracks;
  int m_hits_per_track;
  vector<vector<int>> m_hit_register;
 public:
  IPatternFinder(const char *name, const char *title, int nbins, double xlow, double xup);

  void SetHitsPerTrack(int i);
  
  void AddHits(TGraphErrors &Gr);
  void AddHits(vector<TVector2> &hits);
  void GetPattern(vector<vector<int>> &result);

};

#endif
