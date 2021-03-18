#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TView.h"
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
#include "IPatternFinder.h"

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
  
  TH3D *h3 = new TH3D("h3", "3D View", 20, 0, 2, 20, -1.5, 1.5, 20, -1.5, 1.5);
  h3->GetXaxis()->SetTitle("X");
  h3->GetYaxis()->SetTitle("Y");
  h3->GetZaxis()->SetTitle("Z");

  TH1D *Initial_Momentum = new TH1D("h1_1", "h1_1", 100, 0.5, 100.5);
  TH1D *Final_Momentum = new TH1D("h1_2", "h1_2", 100, 0.5, 100.5);
  TH1D *Deviation = new TH1D("h1_3", "P Difference", 100, 0.5, 100.5);

  TGraphErrors *Bending_Plane = new TGraphErrors();
  TGraphErrors *YZ = new TGraphErrors();
  TGraphErrors *XZ = new TGraphErrors();
  
  const int n = 3;
  const int curve_length = 2;
  const double step = 0.001;
  const bool draw = true;
  
  TPolyLine3D line[n];
  TEllipse Ellipse[n];
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

  for(int iter=0; iter<n; iter++){
    
    colour[iter] = RndmRange(1, 9, rnd);
    x = RndmRange(-0.01, 0.01, rnd);
    y = RndmRange(-0.01, 0.01, rnd);
    z = RndmRange(-0.01, 0.01, rnd);
    H = RndmBinary(rnd);
    initial_p = RndmRange(10, 50, rnd);
    R = initial_p/(B*0.299792458);

    Initial_Momentum->Fill(initial_p);
    
    Helix.SetParameters(x, y, z, R, RndmRange(H*(TMath::PiOver2()-0.2), H*(TMath::PiOver2()+0.2), rnd), H, RndmRange(TMath::Pi()-0.4, TMath::Pi()+0.4, rnd));
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

    //RESET
    for(int d=0; d<5; d++){
      Detector[d].Reset();
    }
    //cout << "Event: " << iter << " done" <<endl;
  }
  
  //cout << "Completed Monte Carlo" << endl;

  //Deviation = Initial_Momentum-(Final_Momentum);

  vector<vector<int>> HitGroups;
  vector<TVector2> positions;
  TVector2 temp;
  TVector3 FitResult;

  double X,Y;

  vector<TVector2> HitVector;

  for(int hit=0; hit<XZ->GetN(); hit++){
    XZ->GetPoint(hit, X, Y);
    temp.SetX(X);
    temp.SetY(Y);
    HitVector.push_back(temp);
  }

  IRiemannFitter Riemann;
  
  IPatternFinder Pattern("p", "angle", 50, -2, 2);
  Pattern.SetHitsPerTrack(5);
  Pattern.AddHits(HitVector);
  Pattern.GetPattern(HitGroups);

  //cout << HitGroups.size() << endl;
  
  for(unsigned int iter=0; iter<HitGroups.size(); iter++){
    for(unsigned int hit_id=0; hit_id<HitGroups[iter].size(); hit_id++){
      Bending_Plane->GetPoint(HitGroups[iter][hit_id], X, Y);
      temp.SetX(X);
      temp.SetY(Y);
      positions.push_back(temp);
      //cout << "Pushing: " << temp.X() << " " << temp.Y() << endl;
    }
    //cout << "Done Pushing" << endl;
    Riemann.SetPoints(positions);
    FitResult = Riemann.Fit();
    cout << "Fitted: " << iter << endl;
    Ellipse[iter].SetX1(FitResult.X());
    Ellipse[iter].SetY1(FitResult.Y());
    Ellipse[iter].SetR1(FitResult.Z());
    Ellipse[iter].SetR2(FitResult.Z());

    cout << "X0: " << FitResult.X() << " Y0: " << FitResult.Y() << " R: " << FitResult.Z() << endl;
    
    final_p = 0.299792458*B*FitResult.Z();

    Final_Momentum->Fill(final_p);
    
    positions.clear();
    positions.shrink_to_fit();
    Riemann.Reset();
  }

  cout << "Chi2: "  <<Initial_Momentum->Chi2Test(Final_Momentum) << endl;
  
  c1->Divide(2,2);
  c1->cd(1);

  h3->Draw();

  for(int c=0; c<n; c++){
    line[c].SetLineColor(colour[c]);
    line[c].Draw();
  }

  c1->cd(2);

  /* Ellipse->SetX1(FitResult[0]);
  Ellipse->SetY1(FitResult[1]);
  Ellipse->SetR1(FitResult[2]);
  Ellipse->SetR2(FitResult[2]);
  Ellipse->SetLineColor(colour[0]);
  Ellipse->SetFillStyle(0);
  */
  Bending_Plane->SetLineWidth(0);
  Bending_Plane->SetMarkerStyle(2);
  Bending_Plane->GetXaxis()->SetTitle("X");
  Bending_Plane->GetYaxis()->SetTitle("Y");
  Bending_Plane->GetXaxis()->SetRangeUser(0, 2);
  Bending_Plane->GetYaxis()->SetRangeUser(-1, 1);
  Bending_Plane->Draw();

  for(unsigned int iter=0; iter<n; iter++){
    Ellipse[iter].SetLineColor(colour[iter]);
    Ellipse[iter].SetFillStyle(0);
    Ellipse[iter].Draw("SAME");
  }
  
  c1->cd(3);

  XZ->SetLineWidth(0);
  XZ->SetMarkerStyle(2);
  XZ->GetXaxis()->SetTitle("X");
  XZ->GetYaxis()->SetTitle("Z");
  XZ->GetXaxis()->SetRangeUser(0, 2);
  XZ->GetYaxis()->SetRangeUser(-1, 1);
  XZ->Draw();

  c1->cd(4);

  YZ->SetLineWidth(0);
  YZ->SetMarkerStyle(2);
  YZ->GetXaxis()->SetTitle("Y");
  YZ->GetYaxis()->SetTitle("Z");
  YZ->GetXaxis()->SetRangeUser(0, 2);
  YZ->GetYaxis()->SetRangeUser(-1, 1);
  YZ->Draw();
  
  c2->Divide(2,1);
  c2->cd(1);

  //Initial_Momentum->SetLineColor(2);
  //Initial_Momentum->Draw();
  //Final_Momentum->Draw("SAME");
  Deviation->Rebin(10);
  Deviation->Draw();

  c2->cd(2);
  Pattern.Draw();

  app.Run(true);
  return 0;
}
