#include "IRiemannFitter.h"

void IRiemannFitter::MapOnSphere(){
  m_points_sphere.clear();
  m_points_sphere.shrink_to_fit();
  double R;
  double Phi;
  TVector3 p;
  for(unsigned int i=0; i<m_points.size(); i++){
    R = m_points[i].Mod();
    Phi = m_points[i].Phi();
    p.SetX(R*TMath::Cos(Phi)/(1+pow(R, 2)));
    p.SetY(R*TMath::Sin(Phi)/(1+pow(R, 2)));
    p.SetZ(pow(R, 2)/(1+pow(R, 2)));
    m_points_sphere.push_back(p);
  }
}

double IRiemannFitter::GetPointP(unsigned int i){
  return pow(1+m_points[i].Mod2(), 2);
}

double IRiemannFitter::GetCovariance(int i, int j){
  TVector3 r0 = this->GetPointsAverages();
  double sum = 0;
  double sumP = 0;
  for(unsigned int a=0; a<m_points_sphere.size(); a++){
    sum += ((m_points_sphere[a](i) - r0(i))*(m_points_sphere[a](j) - r0(j)))*this->GetPointP(a);
    sumP += this->GetPointP(a);
  }
  return sum/sumP;
}

TVector3 IRiemannFitter::GetPointsAverages(){
  TVector3 Avg(0,0,0);
  double sumP = 0;
  for(unsigned int i=0; i<m_points_sphere.size(); i++){
    Avg = Avg + m_points_sphere[i]*this->GetPointP(i);
    sumP += this->GetPointP(i);
  }
  return Avg*(1/sumP);
}

void IRiemannFitter::SetPoints(vector<TVector2> &pos){
  m_points.clear();
  m_points.shrink_to_fit();
  for(unsigned int i=0; i<pos.size(); i++){
    m_points.push_back(pos[i]);
  }
}

void IRiemannFitter::PrintPointsCartesian(){
  if(m_points.size()>0){
    std::cout << "Printing points in cartesian coordinates..." << std::endl;
    for(unsigned int i=0; i<m_points.size(); i++){
      std::cout << "Point: " << i << " X: " << m_points[i].X() << " Y: " << m_points[i].Y() << std::endl;
    }
  }else{
    std::cout << "No points found" << std::endl;
  }
}

void IRiemannFitter::PrintPointsPolar(){
 if(m_points.size()>0){
   std::cout << "Printing points in polar coordinates..." << std::endl;
    for(unsigned int i=0; i<m_points.size(); i++){
      std::cout << "Point: " << i << " R: " << m_points[i].Mod() << " Phi: " << m_points[i].Phi() << std::endl;
    }
  }else{
    std::cout << "No points found" << std::endl;
  }
}

void IRiemannFitter::PrintPointsSphere(){
  if(m_points_sphere.size()<m_points.size()){
    this->MapOnSphere();
  }
  if(m_points_sphere.size()>0){
    std::cout << "Printing points on Riemann sphere..." << std::endl;
    for(unsigned int i=0; i<m_points_sphere.size(); i++){
      std::cout << "Point: " << i << " X: " << m_points_sphere[i].X() << " Y: " << m_points_sphere[i].Y() << " Z: " << m_points_sphere[i].Z() << std::endl;
    }
  }else{
    std::cout << "No points found" << std::endl;
  }
}

TMatrixD IRiemannFitter::GetCovarianceMatrix(){
  TMatrixD Cov(3, 3);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      Cov(i,j) = this->GetCovariance(i, j);
    }
  }
  return Cov;
}

void IRiemannFitter::ScaleAndCenter(double b){
  double sumX = 0, sumY = 0, modX = 0, modY = 0;
  vector<double> Xc;
  vector<double> Yc;
  for(unsigned int i=0; i<m_points.size(); i++){
    sumX += m_points[i].X();
    sumY += m_points[i].Y();
    Xc.push_back(m_points[i].X());
    Yc.push_back(m_points[i].Y());
  }
  for(unsigned int i=0; i<m_points.size(); i++){
    Xc[i] -= (sumX/m_points.size());
    Yc[i] -= (sumY/m_points.size());
    modX += pow(Xc[i], 2);
    modY += pow(Yc[i], 2);
  }
  m_scaling_factor = b/TMath::Sqrt((modX + modY)/((double)m_points.size()));
  cout << "Scaling with scaling factor: " << m_scaling_factor << endl;
  for(unsigned int i=0; i<m_points.size(); i++){
    m_points[i].SetX(m_scaling_factor*Xc[i]);
    m_points[i].SetY(m_scaling_factor*Yc[i]);
  }
  this->MapOnSphere();
}

void IRiemannFitter::Reset(){
  m_points.clear();
  m_points_sphere.clear();
  m_points.shrink_to_fit();
  m_points_sphere.shrink_to_fit();
}

TVector3 IRiemannFitter::Fit(){

  double c;
  TVector3 n, r0, parameters;
  
  if(m_points_sphere.size()<m_points.size()){
    this->MapOnSphere();
  }

  r0 = this->GetPointsAverages();

  const TMatrixD& CovarianceMatrix = this->GetCovarianceMatrix();
  const TMatrixDEigen eigen(CovarianceMatrix);
  const TVectorD EigenValues = eigen.GetEigenValuesRe();
  TMatrixD EigenVectors = eigen.GetEigenVectors();
  
  for(int i=0; i<EigenValues.GetNoElements(); i++){
    if(EigenValues[i]==EigenValues.Min()){
      n.SetX(EigenVectors(0, i));
      n.SetY(EigenVectors(1, i));
      n.SetZ(EigenVectors(2, i));
    }
  }

  c = -1*(n*r0);
											      
  parameters[0] = -1*n.X()/(2*(c+n.Z()));
  parameters[1] = -1*n.Y()/(2*(c+n.Z()));
  parameters[2] = TMath::Sqrt((pow(n.X(), 2) + pow(n.Y(), 2) - 4*c*(c+n.Z()))/(4*pow(c + n.Z(), 2)));
  
  return parameters;
}
