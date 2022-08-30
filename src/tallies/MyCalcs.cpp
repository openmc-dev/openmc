#include <iostream>
#include <cmath>
#include "/home/open_mc/openmc/include/openmc/position.h"
using namespace std;

void boostf( double A[4], double B[4], double X[4]);
double Vdot(double A[4],double B[4]);
void Vcros(double A[4],double B[4],double C[4]);
void Vunit(double A[4] ,double B[4]);

int getMu_lab(double x_det , double y_det , double z_det ,Position p_col , double awr , double incoming_mass )
{
  cout<<"  p col "<<p_col<<endl;
  double r1[4]= {0, x_det, y_det, z_det};  // detector position lab {ignor, x, y, z}
  double r2[4]= {0, p_col.r().x, p_col.r().y, p_col.r().z}; // collision position lab {ignor, x, y, z} 
  double r3[4]; // r1-r2 vector from collision to detector
  double m1= incoming_mass; // mass of incoming particle
  double m2= incoming_mass*awr; // mass of target matirial
  double m3= m1; // mass of outgoing particle to detector
  double m4= m2; // mass of recoil target  system
  double p1[3]={p_col.u().x, p_col.u().y, p_col.u().z}; // 3 momentum of incoming particle
  double p2[3]={0, 0, 0}; //3 momentum of target in lab // need seed 

    // calculate
  double Fp1[4]; //four momentum of incoming particle  
  double Fp2[4]; //four momentum of target
  double Fp3[4]; //four momentum of particle going to the detector
  double UFp3[4]; //unit 3  momentum of particle going to the detector
  double CM[4]; //four momentum of center of mass frame in lab
  double LAB[4]; //four momentum of lab in center of mass frame 
  double mCM; // mass of center of mass system
  double pCM; // momentum of center of mass system
  double p3LAB; // momentum of out particle
  
  double CMFp1[4]; //four momentum of incoming particle in CM  
  double CMFp3[4]; //four momentum of out particle in CM  
  double CMFp4[4]; //four momentum of out target in CM  
  double Fp4[4]; //four momentum of out target in LAB  
  double CMUp1[4]; //unit three vector of incoming particle in CM  
  double CMUp3[4]; //unit three vector of out particle in CM  
  double cosCM; // cosine of out going particle to detector in CM frame
  double cosLAB; // cosine of out going particle to detector in LAB frame
  double CME3; // Energy of out particle in CM
  double CMp3; // momentum of out particle in CM
  double Ur3[4], UCM[4];
  double aa, bb, cc;


  Fp1[0]=0;
  Fp2[0]=0;
  CM[0]=0;
  LAB[0]=0;
  r3[0]=0;


  for(int i=0; i<3; i++){
    Fp1[i+1]=p1[i];
    Fp2[i+1]=p2[i];
    CM[i+1]=Fp1[i+1]+Fp2[i+1];
    LAB[i+1]=-CM[i+1];
    r3[i+1]=r1[i+1]-r2[i+1];
  }
 
  Fp1[0]=sqrt(Fp1[1]*Fp1[1]+Fp1[2]*Fp1[2]+Fp1[3]*Fp1[3]+m1*m1);
  Fp2[0]=sqrt(Fp2[1]*Fp2[1]+Fp2[2]*Fp2[2]+Fp2[3]*Fp2[3]+m2*m2);
  CM[0]=Fp1[0]+Fp2[0];
  LAB[0]=CM[0];
  r3[0]=0;

  mCM=sqrt(CM[0]*CM[0]-CM[1]*CM[1]-CM[2]*CM[2]-CM[3]*CM[3]);
  pCM=sqrt(CM[1]*CM[1]+CM[2]*CM[2]+CM[3]*CM[3]);
  CME3=(mCM*mCM-m4*m4+m3*m3)/(2*mCM); // energy of out going particle  in CM
  CMp3=sqrt(CME3*CME3-m3*m3);

  cout<<mCM<<"  "<<pCM<<"  "<<CME3<<endl;

  boostf(CM,Fp1,CMFp1);

  Vunit(CM,UCM); 
  Vunit(r3,Ur3); 
  cosLAB=Vdot(UCM,Ur3);

  cout<<"cosLAB=  "<<cosLAB<<endl;
  Vunit(CMFp1,CMUp1);


  aa=pCM*pCM*cosLAB*cosLAB-CM[0]*CM[0];
  bb=2*pCM*cosLAB*CME3*mCM;
  cc=CME3*mCM*CME3*mCM-m3*m3*CM[0]*CM[0];

  p3LAB=0;
  p3LAB=(-bb+sqrt(bb*bb-4*aa*cc))/2.0/aa;
  if(p3LAB<=0) p3LAB=(-bb-sqrt(bb*bb-4*aa*cc))/2/aa;
  if(p3LAB<=0) {
    cout<<" detector out of range" <<endl;
    return -1;
  }


  cout<<"p3LAB= "<<p3LAB<<endl;

  Fp3[0]=sqrt(p3LAB*p3LAB+m3*m3);
  for(int i=0; i<3; i++){
    Fp3[i+1]=p3LAB*Ur3[i+1];
  }
  
  boostf(CM,Fp3,CMFp3);
  Vunit(CMFp3,CMUp3); 
  cosCM=Vdot(UCM,CMUp3);

  cout<<"cosCM= "<<cosCM<<endl;


}


void boostf( double A[4], double B[4], double X[4])
{
  //
  //     boosts B(labfram) to A rest frame and gives output in X
  //
  double W;
  int j;
  
  if ((A[0]*A[0]-A[1]*A[1]-A[2]*A[2]-A[3]*A[3])<=0) { 
    cout <<"negative sqrt in boostf"<<A[0]<<A[1]<<A[2]<<A[3]<<endl;}
      
  W=sqrt(A[0]*A[0]-A[1]*A[1]-A[2]*A[2]-A[3]*A[3]);

  if(W==0 || W==(-A[0])) cout <<"divid by 0 in boostf"<<endl;

  X[0]=(B[0]*A[0]-B[3]*A[3]-B[2]*A[2]-B[1]*A[1])/W;
    for(j=1; j<=3; j++) {
      X[j]=B[j]-A[j]*(B[0]+X[0])/(A[0]+W);
    } 
  return;
}




double Vdot(double A[4],double B[4])
{
  int j;

  double dot = 0;

  for(j=1; j<=3; j++) {
    dot = dot + A[j]*B[j];
  }

       
  return dot;   
}


void Vunit(double A[4] ,double B[4])
{
  double fff;
  int j;

  fff = 0;
  
  for(j=1; j<=3; j++) {
    fff = fff + A[j]*A[j];
  }
  
  if (fff==0) {
    cout <<"in vunit divid by zero" << endl;
    return;
  }
  
  for(j=1; j<=3; j++) {
    B[j] = A[j]/sqrt(fff);
  }
  B[0] = 0;
  
  return;   
}


void Vcros(double A[4],double B[4],double C[4])
{
  C[1] = A[2]*B[3]-A[3]*B[2];
  C[2] = A[3]*B[1]-A[1]*B[3];
  C[3] = A[1]*B[2]-A[2]*B[1];
  C[0] = 0;

  if (C[1]==0 && C[2]==0 && C[3]==0.) { 
    cout << "vcross zero" << endl;
  }
  
  return;   
}

