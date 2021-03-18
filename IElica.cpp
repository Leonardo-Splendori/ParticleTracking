#include "IElica.h"

IElica::IElica(){
  m_x0 = 0;
  m_y0 = 0;
  m_z0 = 0;
  m_R = 1;
  m_phi = 0;
  m_h = 1;
  m_lambda = TMath::PiOver4();
}

IElica::IElica(double x0, double y0, double z0, double R, double phi, int h, double lambda){
  m_x0 = x0;
  m_y0 = y0;
  m_z0 = z0;
  m_R = R;
  m_phi = phi;
  m_h = h;
  m_lambda = lambda;
}

void IElica::SetParameters(double x0, double y0, double z0, double R, double phi, int h, double lambda){
  m_x0 = x0;
  m_y0 = y0;
  m_z0 = z0;
  m_R = R;
  m_phi = phi;
  m_h = h;
  m_lambda = lambda;
}

void IElica::SetX0(double x0){
  m_x0 = x0;
}

void IElica::SetY0(double y0){
  m_y0 = y0;
}

void IElica::SetZ0(double z0){
  m_z0 = z0;
}

double IElica::GetX0(){
  return m_x0;
}

double IElica::GetY0(){
  return m_y0;
}

double IElica::GetZ0(){
  return m_z0;
}

void IElica::SetRadius(double R){
  m_R = R;
}

double IElica::GetRadius(){
  return m_R;
}

void IElica::SetPhi(double phi){
  m_phi = phi;
}

double IElica::GetPhi(){
  return m_phi;
}

void IElica::SetH(int h){
  m_h = h;
}

int IElica::GetH(){
  return m_h;
}

void IElica::SetLambda(double lambda){
  m_lambda = lambda;
}

double IElica::GetLambda(){
  return m_lambda;
}

TVector3 IElica::GetTangentVector(double s){
  double c = TMath::Cos(m_lambda);
  TVector3 vect;
  vect.SetX(-1*m_h*c*TMath::Sin(m_phi + m_h*s*c/m_R));
  vect.SetY(m_h*c*TMath::Cos(m_phi + m_h*s*c/m_R));
  vect.SetZ(TMath::Sin(m_lambda));
  return vect;
}

vector<double> IElica::Eval(double s){
  vector<double> pos;
  double c = TMath::Cos(m_lambda);
  
  pos.push_back(m_x0 + m_R*(TMath::Cos(m_phi + s*m_h*c/m_R) - TMath::Cos(m_phi)));
  pos.push_back(m_y0 + m_R*(TMath::Sin(m_phi + s*m_h*c/m_R) - TMath::Sin(m_phi)));
  pos.push_back(m_z0 + s*TMath::Sin(m_lambda));

  return pos;
}

void IElica::PrintParameters(){
  std::cout << "Initial position: (x, y, z) = (" << m_x0 << ", " <<  m_y0 << ", " << m_z0 << ")" << std::endl;
  std::cout << "Radius: " << m_R << std::endl;
  std::cout << "Phi: " << m_phi << std::endl;
  std::cout << "h: " << m_h << std::endl;
  std::cout << "Lambda: " << m_lambda << std::endl;
}

void IElica::GetXIntersection(double x, vector<double> &pos){
  
  pos.clear();
  pos.shrink_to_fit();

  double c = TMath::Cos(m_lambda);
  double s = m_R*(TMath::ACos((x-m_x0)/m_R + TMath::Cos(m_phi))-m_phi)/(m_h*c);
  cout << "s: " << s << endl;

  pos.push_back(m_y0 + m_R*(TMath::Sin(m_phi + s*m_h*c/m_R) - TMath::Sin(m_phi)));
  pos.push_back(m_z0 + s*TMath::Sin(m_lambda));
  
}
