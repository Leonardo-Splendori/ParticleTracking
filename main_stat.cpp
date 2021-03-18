#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TView.h"
#include "TF1.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TPolyLine3D.h"
#include "TApplication.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TEllipse.h"
#include <vector>
#include "IElica.h"
#include "IRilevatore.h"
#include "IRiemannFitter.h"

using namespace std;

double RndmRange(double a, double b, TRandom3 &rnd){
  return rnd.Rndm()*(b-a)+a;
}

double RndmBinary(TRandom3 &rnd){
  double r = rnd.Rndm();
  if(r>0.5){
    return -1;
  }else{
    return 1;
  }
}

int main(){
  TApplication app("app",0,0);
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  
  TRandom3 rnd;
  rnd.SetSeed(123);
  
  TH3D *h3 = new TH3D("h3", "h3", 20, 0, 2, 20, -1.5, 1.5, 20, -1.5, 1.5);
  h3->GetXaxis()->SetTitle("X");
  h3->GetYaxis()->SetTitle("Y");
  h3->GetZaxis()->SetTitle("Z");

  TH1D *Initial_Momentum = new TH1D("h1_1", "Momenti", 100, 0.5, 100.5);
  TH1D *Final_Momentum = new TH1D("h1_2", "Momenti", 100, 0.5, 100.5);
  TH1D *Deviation = new TH1D("h1_3", "Scarto in P", 100, 0.5, 100.5);
  //TH1D *Error = new TH1D("h1_4", "Errore sullo scarto", 100, 0.5, 100.5);
  TGraph *Error = new TGraph();

  TEllipse *Ellipse = new TEllipse();

  TGraphErrors *Bending_Plane = new TGraphErrors();
  TGraphErrors *YZ = new TGraphErrors();
  TGraphErrors *XZ = new TGraphErrors();

  vector<vector<double>> momentums;
  momentums.resize(100);
  
  const int n = 10000;
  const int curve_length = 2;
  const double step = 0.001;
  const bool draw = false;
  
  TPolyLine3D line[n];
  int colour[n];

  double H = RndmBinary(rnd);
  
  IElica Helix;
  
  vector<double> pos;
  int i=0,i1=0;
  double x,y,z;

  IRilevatore Detector[5];
  Detector[0].SetPosition(0.2);
  Detector[1].SetPosition(0.4);
  Detector[2].SetPosition(0.6);
  Detector[3].SetPosition(0.8);
  Detector[4].SetPosition(1); 

  vector<double> intersect;

  double initial_p;
  double final_p;
  double R;
  double B = 3;

  vector<TVector2> positions;
  TVector2 temp;

  IRiemannFitter Riemann;

  TVector3 FitResult;

  for(int iter=0; iter<n; iter++){
    
    colour[iter] = RndmRange(1, 9, rnd);
    x = RndmRange(-0.01, 0.01, rnd);
    y = RndmRange(-0.01, 0.01, rnd);
    z = RndmRange(-0.01, 0.01, rnd);
    H = RndmBinary(rnd);
    initial_p = RndmRange(1, 100, rnd);
    if(initial_p>100){
      break;
    }
    R = initial_p/(B*0.299792458);

    Initial_Momentum->Fill(initial_p);
    
    Helix.SetParameters(x, y, z, R, RndmRange(H*(TMath::PiOver2()-0.2), H*(TMath::PiOver2()+0.2), rnd), H, RndmRange(TMath::Pi()-0.2, TMath::Pi()+0.2, rnd));
    //Helix.PrintParameters();

    i=0;
    if(draw){
      for(double s=0; s<curve_length; s += step){
	pos = Helix.Eval(s);
	line[iter].SetPoint(i,pos[0],pos[1],pos[2]);
	i++;
      }
    }
    
    for(int d=0; d<5; d++){
      Detector[d].GenerateHits(0, curve_length, Helix, Detector[d].GetDepth());
    }
    
    i=0;
    for(int d=0; d<5; d++){
      for(int e=0; e<Detector[d].GetNumberOfHits(); e++){
	Bending_Plane->SetPoint(i1, Detector[d].GetPosition(), Detector[d].GetHitY(e));
	Bending_Plane->SetPointError(i1, Detector[d].GetDepth(), Detector[d].GetResolution());
	YZ->SetPoint(i1, Detector[d].GetHitY(e), Detector[d].GetHitZ(e));
	YZ->SetPointError(i1, Detector[d].GetResolution(), Detector[d].GetResolution());
	XZ->SetPoint(i1, Detector[d].GetPosition(), Detector[d].GetHitZ(e));
	XZ->SetPointError(i1, Detector[d].GetDepth(), Detector[d].GetResolution());
	i1++;
      }
    }
    
    for(int d=0; d<5; d++){
      temp.SetX(Detector[d].GetPosition());
      temp.SetY(Detector[d].GetHitY(0));
      positions.push_back(temp);
    }
    
    Riemann.SetPoints(positions);
    
    FitResult = Riemann.Fit();
    
    final_p = 0.299792458*B*FitResult.Z();

    // cout << floor(final_p-1) << endl;
    momentums[floor(final_p-1)].push_back(final_p-initial_p);

    Final_Momentum->Fill(final_p);

    //RESET
    Riemann.Reset();
    for(int d=0; d<5; d++){
      Detector[d].Reset();
    }
    positions.clear();
    positions.shrink_to_fit();
    cout << "Event: " << iter << " done" <<endl;
  }
  
  cout << "Completed Monte Carlo" << endl;

  vector<double> avg;
  vector<double> deviations;
  avg.resize(100);
  deviations.resize(100);

  for(unsigned int iter=0; iter<avg.size(); iter++){
    for(unsigned int internal=0; internal<momentums[iter].size(); internal++){
      avg[iter] += momentums[iter][internal];
    }
    avg[iter] /= momentums[iter].size();
  }

  for(unsigned int iter=0; iter<deviations.size(); iter++){
    for(unsigned int internal=0; internal<momentums[iter].size(); internal++){
      deviations[iter] += pow(momentums[iter][internal]-(avg[iter]), 2);
    }
    deviations[iter] /= (momentums[iter].size());
    deviations[iter] = sqrt(deviations[iter]);
    //deviations[iter] /= sqrt(momentums[iter].size());
  }

  int n_bins = Initial_Momentum->GetNbinsX();
  for(int k=0; k<n_bins; k++){
    Deviation->AddBinContent(k+1, (abs(Initial_Momentum->GetBinContent(k) - Final_Momentum->GetBinContent(k))));
    //Error->AddBinContent(k+1, deviations[k]/(k+1));
    Error->SetPoint(k, k+1, deviations[k]/(k+1));
  }

  cout << "Chi2: " << Initial_Momentum->Chi2Test(Final_Momentum) << endl;

  TF1 *f = new TF1("f", "[0]*x+[1]", 0, 100);
  
  c1->Divide(2,2);
  c1->cd(1);

  h3->Draw();

  for(int c=0; c<n; c++){
    line[c].SetLineColor(colour[c]);
    line[c].Draw();
  }

  c1->cd(2);

  Ellipse->SetX1(FitResult[0]);
  Ellipse->SetY1(FitResult[1]);
  Ellipse->SetR1(FitResult[2]);
  Ellipse->SetR2(FitResult[2]);
  Ellipse->SetLineColor(colour[0]);
  Ellipse->SetFillStyle(0);

  Bending_Plane->SetLineWidth(0);
  Bending_Plane->SetMarkerStyle(2);
  Bending_Plane->GetXaxis()->SetTitle("X");
  Bending_Plane->GetYaxis()->SetTitle("Y");
  Bending_Plane->GetXaxis()->SetRangeUser(0, 2);
  Bending_Plane->GetYaxis()->SetRangeUser(-1, 1);
  Bending_Plane->Draw();
  Ellipse->Draw("SAME");

  c1->cd(3);

  XZ->SetLineWidth(0);
  XZ->SetMarkerStyle(2);
  XZ->GetXaxis()->SetTitle("X");
  XZ->GetYaxis()->SetTitle("Z");
  XZ->Draw();

  c1->cd(4);

  YZ->SetLineWidth(0);
  YZ->SetMarkerStyle(2);
  YZ->GetXaxis()->SetTitle("Y");
  YZ->GetYaxis()->SetTitle("Z");
  YZ->Draw();
  
  c2->Divide(2,1);
  c2->cd(1);

  Initial_Momentum->SetLineColor(2);
  Initial_Momentum->GetXaxis()->SetTitle("p[GeV]");
  Initial_Momentum->GetYaxis()->SetTitle("N");
  Initial_Momentum->Draw();
  Final_Momentum->Draw("SAME");

  c2->cd(2);

  Error->GetXaxis()->SetTitle("p[GeV]");
  Error->GetYaxis()->SetTitle("s(p)/p");
  Error->GetXaxis()->SetRangeUser(0, 100);
  Error->SetLineWidth(0);
  Error->SetMarkerStyle(20);
  Error->SetNameTitle("Errore su p", "Errore su p");
  Error->Fit(f);
  Error->Draw();
  f->Draw("SAME");

  app.Run(true);
  return 0;
}
