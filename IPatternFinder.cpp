#include "IPatternFinder.h"

IPatternFinder::IPatternFinder(const char *name, const char *title, int nbins, double xlow, double xup) //consigliato xlow 1.5 xup 5
  : TH1D(name, title, nbins, xlow, xup){
  m_ntracks = 0;
  m_hits_per_track = 5;
  m_hit_register.resize(nbins + 10);
}

void IPatternFinder::SetHitsPerTrack(int i){
  if(i>0){
    m_hits_per_track = i;
  }else{
    std::cout << "Number of expected hits per track must be greater than 0" << std::endl;
  }
}

void IPatternFinder::AddHits(vector<TVector2> &hits){
  double arg, x, y;
  for(unsigned int i=0; i<hits.size(); i++){
    x = hits[i].X();
    y = hits[i].Y();
    if(x!=0){
      arg = atan(y/x);
    }else if(y>0){
      arg = TMath::PiOver2();
    }else if(y<0){
      arg = -1*TMath::PiOver2();
    }else if(y==0){
      arg = 0;
    }
    this->Fill(arg);
    //cout << arg << endl;
    //cout << "Adding to bin: " << this->FindBin(arg) <<endl;
    if((unsigned int)this->FindBin(arg)>m_hit_register.size()){
      m_hit_register.resize(this->FindBin(arg));
    }
    m_hit_register[this->FindBin(arg)].push_back(i);
  }
}

void IPatternFinder::AddHits(TGraphErrors &Gr){
  double X,Y,arg;
  for(int i=0; i<Gr.GetN(); i++){
    Gr.GetPoint(i, X, Y);
    arg = TMath::Pi() + TMath::ATan2(-1*Y, -1*X);
    this->Fill(arg);
    m_hit_register[this->FindBin(arg)].push_back(i);
  }
}

void IPatternFinder::GetPattern(vector<vector<int>> &result){
  int c = 0;
  for(int i=0; i<this->GetNbinsX(); i++){
    if(this->GetBinContent(i)>0){
      for(int j=0; j<this->GetBinContent(i); j++){
	result.resize(c+1);
	result[c].push_back(m_hit_register[i][j]);
      }
      c++;
      //cout << "HitGroupSize: " << this->GetBinContent(i) << endl;
    }
  }
}
