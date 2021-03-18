#include "IRilevatore.h"

IRilevatore::IRilevatore(){
  m_position = 1;
  m_width = 1;
  m_height = 1;
  m_depth = 150*1E-6;
  m_resolution = 50*1E-6;
}

IRilevatore::IRilevatore(double position, double width, double height, double depth, double resolution){
  m_position = position;
  m_width = width;
  m_height = height;
  m_depth = depth;
  m_resolution = resolution;
}

void IRilevatore::SetHits(vector<vector<double>> &hits){
  m_hits.resize(hits.size());
  for(unsigned int i=0; i<hits.size(); i++){
    m_hits[i].resize(hits[i].size());
    for(unsigned int j=0; j<hits[i].size(); j++){
      m_hits[i][j] = hits[i][j];
    }
  }
}

void IRilevatore::SetHit(double y, double z){
  m_hits.resize(m_hits.size()+1);
  m_hits[m_hits.size()-1].resize(2);
  m_hits[m_hits.size()-1][0] = y;
  m_hits[m_hits.size()-1][1] = z;
}

void IRilevatore::SetPosition(double position){
  m_position = position;
}

void IRilevatore::SetWidth(double width){
  m_width = width;
}

void IRilevatore::SetHeight(double height){
  m_height = height;
}

void IRilevatore::SetDepth(double depth){
  m_depth = depth;
}

void IRilevatore::SetResolution(double resolution){
  m_resolution = resolution;
}

void IRilevatore::GetHits(vector<vector<double>> &hits){
  hits.resize(m_hits.size());
  for(unsigned int i=0; i<m_hits.size(); i++){
    hits[i].resize(m_hits[i].size());
    for(unsigned int j=0; j<m_hits[i].size(); j++){
      hits[i][j] = m_hits[i][j];
    }
  }
}

double IRilevatore::GetHitY(unsigned int i){
  return m_hits[i][0];
}

double IRilevatore::GetHitZ(unsigned int i){
  return m_hits[i][1];
}

double IRilevatore::GetPosition(){
  return m_position;
}

double IRilevatore::GetWidth(){
  return m_width;
}

double IRilevatore::GetHeight(){
  return m_height;
}

double IRilevatore::GetDepth(){
  return m_depth;
}

double IRilevatore::GetResolution(){
  return m_resolution;
}

int IRilevatore::GetNumberOfHits(){
  return m_hits.size();
}

void IRilevatore::GenerateHits(double smin, double smax, IElica &Elica, double step){
  vector<double> point;
  for(double s=smin; s<smax; s+=step){
    point = Elica.Eval(s);
    if(point[0]>=m_position){
      this->SetHit(point[1], point[2]);
      break;
    }
  }
}

void IRilevatore::Reset(){
  for(unsigned int i=0; i<m_hits.size(); i++){
    m_hits[i].clear();
    m_hits[i].shrink_to_fit();
  }
  m_hits.clear();
  m_hits.shrink_to_fit();
}
