#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"Test1.h"

#define N_PTCL        5000        //系内の全粒粒子数
#define N_Satb_PTCL   0         //考えている系の前後にそれぞれ配置する粒子数
#define N_ALL   N_PTCL + N_Satb_PTCL
#define N1  N_PTCL * Rho1 / (Rho1 + Rho2)
#define N2  N_PTCL * Rho2 / (Rho1 + Rho2)
#define dl_N1  2 / (sqrt(9 + 8 * N1) - 5)
#define dl_N2  2 / (sqrt(9 + 8 * N2) + 1)

#define Gamma   1.4              //比熱比
#define StepN   5000             //ステップ数
#define TOL     0.000001          //p_starの収束条件
#define DIM     2                 //二次元
#define Mesh    200               //網の目の数

void Current_time();

double KernelFunc  ( int i,int j, double r[][DIM], double h[] );
double DifferKernelFunc ( int i, int j, int k, double r[][DIM], double h[] ); //kはベクトル成分
double DifferKernelFunc_hj ( int i, int j, int k, double r[][DIM], double h[] ); //kはベクトル成分
double Grad_h_term (int i, int k,double m[], double r[][DIM], double d, double h[]);

void Mesh_Condi (double mesh_r[Mesh+1][Mesh+1][DIM]);
void Mesh_Density (double m[], double r[][DIM], double mesh_d[Mesh+1][Mesh+1], double mesh_r[Mesh+1][Mesh+1][DIM]);
double Mesh_KernelFunc (int k, double x, double y, double r[][DIM], double h );
void Mesh_PrintData    (FILE *file, double r[Mesh+1][Mesh+1][DIM], double d[Mesh+1][Mesh+1]);
void Mesh_Periodic_Boundary_Conditions(double mesh_x, double mesh_y, double r[][DIM]);
void Periodic_Boundary_Conditions(int i, double r[][DIM]);
void Return_Periodic_Boundary_Conditions(double r[][DIM]);

double ViscosityTerm ( int i, int j, double r[][DIM], double v[][DIM], double d[], double u[], double h[]);
double Smoothing_Length ( int i, double m[], double r[][DIM],double d[]);
double Time_Step   ( double r[][DIM], double v[][DIM], double p[], double d[], double u[], double h[]);

void InicialCondi  ( double m[], double r[][DIM], double v[][DIM], double a[][DIM], double p[], double d[], double u[], double du[], double h[] );
void RungeKutta    ( double dt, double m[], double r[][DIM], double v[][DIM], double a[][DIM], double p[], double d[], double u[], double du[], double h[] );
void PrintData     ( FILE *file, double m[], double r[][DIM], double v[][DIM], double a[][DIM], double p[], double d[], double u[], double du[], double h[] );

//double用比較関数
int asc(const void *a, const void *b) ;

//配列内の最大値・最小値を返す関数
double MaxArray( double n[], int l );
double MinArray( double n[], int l );
double Compare_Max( double x, double y);


int main (void)
{
  Current_time();
  clock_t start,end;
  start = clock();
  //質量、位置、速度、加速度、圧力、密度、内部エネルギー、内部エネルギーの時間変化、影響半径
  double mass[N_ALL + 1], pos[N_ALL+1][DIM], vel[N_ALL+1][DIM], acc[N_ALL+1][DIM], press[N_ALL+1],
                dens[N_ALL+1], energy[N_ALL+1], difenergy[N_ALL+1], len[N_ALL+1];
  //画像表示用の配列準備
  static double mesh_pos[Mesh+1][Mesh+1][DIM], mesh_dens[Mesh+1][Mesh+1];
  double timestep, totaltime = 0.0;
  FILE *fp;
  FILE *fl;
  int i = 0;
  char filename[100], filename_2[100];
  for(i = 0; i < StepN; i++){
    sprintf(filename, "data_%d.ptcl",i);
    sprintf(filename_2,"Mesh_data_%d.ptcl",i);
    fp = fopen(filename, "w");
    fl = fopen(filename_2,"w");
    if( fp == NULL || fl == NULL) {
      perror("ファイルの読み込みに失敗！\n");
      return 1;
    }
    if(i==0){
      InicialCondi(mass, pos, vel, acc, press, dens, energy, difenergy, len);
      PrintData(fp, mass, pos, vel, acc, press, dens, energy, difenergy, len);
      Mesh_Condi(mesh_pos);
      Mesh_Density(mass, pos, mesh_dens, mesh_pos);
      Mesh_PrintData(fl, mesh_pos, mesh_dens);
    }else{
      //timestepの計算
      timestep = Time_Step(pos, vel, press, dens, energy, len);
      totaltime = totaltime + timestep;
      if (totaltime > EndTime){
        totaltime = totaltime - timestep;
        timestep = EndTime - totaltime;
        RungeKutta(timestep, mass, pos, vel, acc, press, dens, energy, difenergy, len);
        PrintData(fp, mass, pos, vel, acc, press, dens, energy, difenergy, len);
        fclose(fp);

        Mesh_Density(mass, pos,mesh_dens,mesh_pos);
        Mesh_PrintData(fl, mesh_pos, mesh_dens);
        fclose(fl);
        end = clock();
        printf("%.2f秒\n",(double)(end-start)/CLOCKS_PER_SEC);
        exit(0);
      }
      RungeKutta(timestep, mass, pos, vel, acc, press, dens, energy, difenergy, len);
      PrintData(fp, mass, pos, vel, acc, press, dens, energy, difenergy, len);
      Mesh_Density(mass, pos, mesh_dens, mesh_pos);
      Mesh_PrintData(fl, mesh_pos, mesh_dens);
    }
    fclose(fp);
    fclose(fl);
  }
  return 0;
}

void Mesh_Condi (double mesh_r[Mesh+1][Mesh+1][DIM])
{
  int i,j;
  double a = 0.0, b = 0.0;
  double dx = XMax / Mesh;
  double dy = YMax / Mesh;
  for(i = 0; i < Mesh+1; i++){
    b = 0.0;
    for(j = 0; j < Mesh+1; j++){
      mesh_r[i][j][1] = a*dx;
      mesh_r[i][j][2] = b*dy;
      b =  b +1.0 ;
    }
    a = a + 1.0;
  }
}
//周期的境界条件を満たすように、i粒子の位置によって周りの粒子の位置を変更する
void Periodic_Boundary_Conditions(int i, double r[][DIM])
{
  int k;
  if (r[i][1] > 0.75){  //i粒子が0.75よりも大きい位置にいる場合、0<ｘ<0.25位置にある粒子をx=1よりも大きい領域になければならない
    if (r[i][2] > 0.75){  //上右端
      for (k=0;k<N_ALL;k++){
        if (r[k][2] < 0.25 && r[k][1] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
          r[k][1] = r[k][1] + XMax + dl_N1; //0<x<0.25 の粒子を1<xに置く →左端の粒子を右端に置く
        }else if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1] < 0.25 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] + XMax + dl_N2;  //Rho2領域
        }else if (r[k][1] < 0.25 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] + XMax + dl_N1 ;  //Rho1領域
        }
      }
    }else if (r[i][2] < 0.25){    //下右端
      for (k=0;k<N_ALL;k++){
        if (r[k][2] > 0.75 && r[k][1] < 0.25){
          r[k][2] = r[k][2] - YMax - dl_N1; //0.75<y<1 粒子をy<0に置く
          r[k][1] = r[k][1] + XMax + dl_N1; //0<x<0.25 の粒子を1<xに置く →左上端の粒子を右下端に置く
        }else if (r[k][2] > 0.75){
          r[k][2] = r[k][2] - YMax - dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1] < 0.25 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] + XMax + dl_N2;
        }else if (r[k][1] < 0.25 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] + XMax + dl_N1 ;
        }
      }
    }else {
      for (k=0;k<N_ALL;k++){
        if (r[k][1] < 0.25 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] + XMax + dl_N2;
        }else if (r[k][1] < 0.25 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] + XMax + dl_N1 ;
        }
      }
    }

  }else if (r[i][1] < 0.25){
    if (r[i][2] > 0.75){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] < 0.25 && r[k][1] > 0.75){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
          r[k][1] = r[k][1] - XMax - dl_N1; //0.75<x<1 の粒子をx<0に置く →右下端の粒子を左上端に
        }else if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1]> 0.75 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] - XMax - dl_N2;
        }else if (r[k][1] > 0.75 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] - XMax - dl_N1;
        }
      }
    }else if(r[i][2] < 0.25){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] > 0.75 && r[k][1] > 0.75){
          r[k][2] = r[k][2] - YMax - dl_N1; //0.75<y<1 粒子をy<0に置く
          r[k][1] = r[k][1] - XMax - dl_N1; //0.75<x<1 の粒子をx<0に置く →右上端の粒子を左下端に置く
        }else if (r[k][2] > 0.75){
          r[k][2] = r[k][2] - YMax - dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1]> 0.75 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] - XMax - dl_N2;
        }else if (r[k][1] > 0.75 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] - XMax - dl_N1;
        }
      }
    }else{
      for(k = 0; k < N_ALL; k++){
        if (r[k][1]> 0.75 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] - XMax - dl_N2;
        }else if (r[k][1] > 0.75 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] - XMax - dl_N1;
        }
      }
    }
  }else {
    if(r[i][2] < 0.25){
      for(k = 0; k < N_ALL; k++){
        if ( r[k][2] > 0.75 ){
          r[k][2] = r[k][2] - YMax - dl_N1  ;
          printf("*\n" );
        }
      }
    }else if (r[i][2] > 0.75){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1;
          printf("**\n" );
        }
      }
    }
  }
}

void Mesh_Periodic_Boundary_Conditions(double mesh_x, double mesh_y, double r[][DIM])
{
  int k;
  if (mesh_x > 0.75){  //i粒子が0.75よりも大きい位置にいる場合、0<ｘ<0.25位置にある粒子をx=1よりも大きい領域になければならない
    if (mesh_y > 0.75){  //上右端
      for (k=0;k<N_ALL;k++){
        if (r[k][2] < 0.25 && r[k][1] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
          r[k][1] = r[k][1] + XMax + dl_N1; //0<x<0.25 の粒子を1<xに置く →左端の粒子を右端に置く
        }else if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1] < 0.25 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] + XMax + dl_N2;  //Rho2領域
        }else if (r[k][1] < 0.25 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] + XMax + dl_N1 ;  //Rho1領域
        }
      }
    }else if (mesh_y < 0.25){    //下右端
      for (k=0;k<N_ALL;k++){
        if (r[k][2] > 0.75 && r[k][1] < 0.25){
          r[k][2] = r[k][2] - YMax - dl_N1; //0.75<y<1 粒子をy<0に置く
          r[k][1] = r[k][1] + XMax + dl_N1; //0<x<0.25 の粒子を1<xに置く →左上端の粒子を右下端に置く
        }else if (r[k][2] > 0.75){
          r[k][2] = r[k][2] - YMax - dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1] < 0.25 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] + XMax + dl_N2;
        }else if (r[k][1] < 0.25 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] + XMax + dl_N1 ;
        }
      }
    }else {
      for (k=0;k<N_ALL;k++){
        if (r[k][1] < 0.25 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] + XMax + dl_N2;
        }else if (r[k][1] < 0.25 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] + XMax + dl_N1 ;
        }
      }
    }

  }else if (mesh_x < 0.25){
    if (mesh_y > 0.75){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] < 0.25 && r[k][1] > 0.75){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
          r[k][1] = r[k][1] - XMax - dl_N1; //0.75<x<1 の粒子をx<0に置く →右下端の粒子を左上端に
        }else if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1]> 0.75 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] - XMax - dl_N2;
        }else if (r[k][1] > 0.75 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] - XMax - dl_N1;
        }
      }
    }else if(mesh_y < 0.25){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] > 0.75 && r[k][1] > 0.75){
          r[k][2] = r[k][2] - YMax - dl_N1; //0.75<y<1 粒子をy<0に置く
          r[k][1] = r[k][1] - XMax - dl_N1; //0.75<x<1 の粒子をx<0に置く →右上端の粒子を左下端に置く
        }else if (r[k][2] > 0.75){
          r[k][2] = r[k][2] - YMax - dl_N1; //0<y<0.25の粒子を1<yに置く
        }
        if (r[k][1]> 0.75 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] - XMax - dl_N2;
        }else if (r[k][1] > 0.75 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] - XMax - dl_N1;
        }
      }
    }else{
      for(k = 0; k < N_ALL; k++){
        if (r[k][1]> 0.75 && r[k][2] > 0.25 && r[k][2] < 0.75){
          r[k][1] = r[k][1] - XMax - dl_N2;
        }else if (r[k][1] > 0.75 && (r[k][2] < 0.25 || r[k][2] > 0.75) ){
          r[k][1] = r[k][1] - XMax - dl_N1;
        }
      }
    }
  }else {
    if(mesh_y < 0.25){
      for(k = 0; k < N_ALL; k++){
        if ( r[k][2] > 0.75 ){
          r[k][2] = r[k][2] - YMax - dl_N1  ;
        }
      }
    }else if (mesh_y > 0.75){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax + dl_N1;
        }
      }
    }
  }
}

void Mesh_Density (double m[], double r[][DIM], double mesh_d[Mesh+1][Mesh+1], double mesh_r[Mesh+1][Mesh+1][DIM])
{
  int i, j, k;
  double h = 0.03; //smoothing lengthの目安　→hが収束するまで、密度計算、h計算を繰り返す

  for (i = 0; i < Mesh+1; i++){
    for(j = 0; j < Mesh+1; j++){
      Mesh_Periodic_Boundary_Conditions(mesh_r[i][j][1] , mesh_r[i][j][2], r);
      mesh_d[i][j] = 0.0;
      for(k = 0; k<N_ALL;k++)
        mesh_d[i][j] = mesh_d[i][j] + m[k] * Mesh_KernelFunc(k, mesh_r[i][j][1] , mesh_r[i][j][2], r, h);
    }
    Return_Periodic_Boundary_Conditions(r);
  }
}


void Return_Periodic_Boundary_Conditions(double r[][DIM])
{
  int k;
  for(k = 0; k < N_ALL; k++){
    if (r[k][1] > 1.0 && (r[k][2] < 0.25 || r[k][2] > 0.75)){
      r[k][1] = r[k][1] - XMax - dl_N1;
    }else if (r[k][1] > 1.0 && r[k][2] > 0.25 && r[k][2] < 0.75){
      r[k][1] = r[k][1] - XMax - dl_N2;
    }else if(r[k][1] < 0.0 && (r[k][2] < 0.25 || r[k][2] > 0.75)){
      r[k][1] = r[k][1] + XMax + dl_N1;
    }else if (r[k][1] < 0.0 && r[k][2] > 0.25 && r[k][2] < 0.75){
      r[k][1] = r[k][1] + XMax + dl_N2;
    }else if (r[k][2] > 1.0){
      r[k][2] = r[k][2] - YMax - dl_N1;
    }else if (r[k][2] < 0.0){
      r[k][2] = r[k][2] + YMax + dl_N1;
    }
  }
}
//カーネル関数Wijの定義　粒子の位置とメッシュの位置が違うVer.
double Mesh_KernelFunc (int k, double x, double y, double r[][DIM], double h )
{
  double dx, dy;
  dx = x - r[k][1];
  dy = y - r[k][2];
  return exp(-1.0*(dx*dx + dy*dy)/h/h) /h/h/M_PI;
}

//カーネル関数Wijの定義
double KernelFunc ( int i, int j, double r[][DIM], double h[] )
{
  double dx, dy;
  dx = r[i][1] - r[j][1];
  dy = r[i][2] - r[j][2];
  return exp(-1.0*(dx*dx + dy*dy)/h[i]/h[i]) /h[i]/h[i]/M_PI;
}
//カーネル関数の微分∇Wijの定義
double DifferKernelFunc ( int i, int j, int k, double r[][DIM], double h[] ) //kはベクトル成分
{
  double dx,dy;
  dx = r[i][1] - r[j][1];
  dy = r[i][2] - r[j][2];
  if (k==1){
    return -2.0*dx*exp(-1.0*(dx*dx+dy*dy)/h[i]/h[i])/h[i]/h[i]/h[i]/h[i]/M_PI;
  }else if (k == 2){
    return -2.0*dy*exp(-1.0*(dx*dx+dy*dy)/h[i]/h[i])/h[i]/h[i]/h[i]/h[i]/M_PI;
  }
}
//カーネル関数の微分∇Wij(hj)の定義
double DifferKernelFunc_hj ( int i, int j, int k, double r[][DIM], double h[]) //kはベクトル成分
{
  double dx,dy;
  dx = r[i][1] - r[j][1];
  dy = r[i][2] - r[j][2];
  if (k==1){
    return -2.0*dx*exp(-1.0*(dx*dx+dy*dy)/h[j]/h[j])/h[j]/h[j]/h[j]/h[j]/M_PI;
  }else if (k == 2){
    return -2.0*dy*exp(-1.0*(dx*dx+dy*dy)/h[j]/h[j])/h[j]/h[j]/h[j]/h[j]/M_PI;

  }
}
//人工粘性項の定義
double ViscosityTerm ( int i, int j, double r[][DIM], double v[][DIM], double d[], double u[], double h[])
{
  double   e = 0.01;
  double dx = 0.0, dy = 0.0, dvx = 0.0, dvy = 0.0, c = 0.0, mu = 0.0, w = 0.0, v_sig = 0.0;
  dx = r[i][1] - r[j][1];
  dy = r[i][2] - r[j][2];
  dvx = v[i][1] - v[j][1];
  dvy = v[i][2] - v[j][2];
  //音速の計算
  w = (dx * dvx + dy * dvy) / sqrt(dx*dx+dy*dy);
  v_sig =  sqrt((Gamma - 1.0) * u[i]) + sqrt((Gamma-1.0) * u[j]) - 3* w;
  //c = (sqrt((Gamma - 1.0) * u[i]) + sqrt((Gamma-1.0) * u[j])) / 2.0;
  //mu = h[i]*dv*dx /(dx * dx + Epsilon * h[i] * h[i]);

  if ( w < 0){
    //return 2.0*(-1.0 * Alpha * c * mu + Beta * mu * mu) / (d[i] + d[j]);
    return - 1.0 * Alpha * v_sig * w /(d[i]+d[j]);
  }else{
    return 0.0;
  }
}
//Smoothing Length の計算　i粒子に対する近傍粒子数に応じて長さを変える
double Smoothing_Length ( int i, double m[], double r[][DIM], double d[])
{
  return sqrt(Neighbor_PTCL * m[i] / d[i]);
  /*
  int j = 0;
  double dx[N_PTCL], dy[N_PTCL];
  for( j = 0; j < N_PTCL; j++){
    dx[j] = fabs(r[i][1] - r[j][1]);
    dy[j] = fabs(r[i][2] - r[j][2]);
  }
  qsort(dx, N_PTCL, sizeof(double), asc);
  return dx[Neighbor_PTCL];
  */
}
//最も小さいdtの計算
double Time_Step ( double r[][DIM], double v[][DIM], double p[], double d[], double u[], double h[])
{
  double t[N_PTCL], mu[N_PTCL];
  int i = 0, j = 0;
  double dx = 0.0, dy = 0.0, dvx = 0.0, dvy = 0.0,c = 0.0;

  for ( i = 0; i < N_PTCL; i++){
    c = sqrt((Gamma-1) * u[i]);
    for( j = 0; j < N_PTCL; j++){
      dx = r[i][1] - r[j][1];
      dy = r[i][2] - r[j][2];
      dvx = v[i][1] - v[j][1];
      dvy = v[i][2] - v[j][2];
      mu[j] = h[i] * (dvx * dx  + dvy * dy) / (dx * dx + dy * dy + Epsilon * h[i] * h[i]);
    }
    t[i] = h[i] / ( c + 0.6 * (Alpha * c + Beta * MaxArray(mu, (int) sizeof(mu) / sizeof(mu[0])) ));
  }
  return  0.25 * MinArray(t, (int) sizeof(t) / sizeof(t[0]));
}
//初期条件(密度、位置、速さ、圧力、加速度)の定義：密度を1,0.25と仮定→hを固定
void InicialCondi  ( double m[], double r[][DIM], double v[][DIM], double a[][DIM], double p[], double d[], double u[], double du[], double h[] )
{
  int i = 0, j = 0;

  double difker_i_x = 0.0 , difker_j_x = 0.0, difker_ij_x = 0.0 , difker_i_y = 0.0 , difker_j_y = 0.0, difker_ij_y = 0.0 , vis = 0.0;

  for (i = 0; i < N_ALL; i++ ){
    m[i] = (Rho1 + Rho2) /2.0/N_PTCL;
    //初期化
    a[i][1] = 0.0;
    a[i][2] = 0.0;
    du[i] = 0.0;
  }
  for(i=0; i<N_ALL; i++){
    if (i<= N1/2){
      if (i==0){
        r[0][1] = XMin;
        r[0][2] = YMin;
      }else{
        r[i][1] = r[i-1][1] + dl_N1;
        r[i][2] = r[i-1][2];
        if (r[i][1] > XMax){
          r[i][1] = 0.0;
          r[i][2] = r[i][2] + dl_N1;
        }
      }
      p[i] = Press1;
      d[i] = Rho1;
      v[i][1] = Vel1;
      v[i][2] = W0 *sin(2*M_PI*r[i][1]/Lambda)*(exp(-(r[i][2]-0.25)*(r[i][2]-0.25)/2/Sigma/Sigma)+exp(-(r[i][2]-0.75)*(r[i][2]-0.75)/2/Sigma/Sigma));
    }else if (N1 / 2 < i && i <= N2 + N1 / 2){
      r[i][1] = r[i-1][1] + dl_N2;
      r[i][2] = r[i-1][2];
      if (r[i][1] > XMax){
        r[i][1] = 0.0;
        r[i][2] = r[i][2] + dl_N2;
      }
      p[i] = Press2;
      d[i] = Rho2;
      v[i][1] = Vel2;
      v[i][2] = W0 *sin(2*M_PI*r[i][1]/Lambda)*(exp(-(r[i][2]-0.25)*(r[i][2]-0.25)/2/Sigma/Sigma)+exp(-(r[i][2]-0.75)*(r[i][2]-0.75)/2/Sigma/Sigma));
    }else{
      r[i][1] = r[i-1][1] + dl_N1;
      r[i][2] = r[i-1][2];
      if(r[i][1] > XMax){
        r[i][1] = 0.0;
        r[i][2] = r[i][2] + dl_N1;
      }
      p[i] = Press1;
      d[i] = Rho1;
      v[i][1] = Vel1;
      v[i][2] = W0 *sin(2*M_PI*r[i][1]/Lambda)*(exp(-(r[i][2]-0.25)*(r[i][2]-0.25)/2/Sigma/Sigma)+exp(-(r[i][2]-0.75)*(r[i][2]-0.75)/2/Sigma/Sigma));
    }
  }
  //smoothing length の計算
  for(i = 0; i < N_ALL; i++){
    h[i] = Smoothing_Length(i, m, r, d);
  }
  //内部エネルギーの初期値
  for(i = 0; i < N_ALL; i++){
    u[i] = p[i] / (Gamma -1.0) / d[i];
  }

  //加速度とエネルギーの時間微分の計算
  for(i = 0; i < N_ALL; i++){
    for(j = 0; j < N_ALL; j++){
      difker_i_x =  DifferKernelFunc(i, j, 1, r, h);
      difker_j_x =  DifferKernelFunc_hj(i, j, 1, r, h);
      difker_ij_x = (difker_i_x+difker_j_x)/2.0;
      difker_i_y =  DifferKernelFunc(i, j, 2, r, h);
      difker_j_y =  DifferKernelFunc_hj(i, j, 2, r, h);
      difker_ij_y = (difker_i_y+difker_j_y)/2.0;
      a[i][1] = a[i][1] + (-1.0) * m[j] * (p[i] / d[i] / d[i] + p[j] / d[j] / d[j] ) * difker_ij_x;
      a[i][2] = a[i][2] + (-1.0) * m[j] * (p[i] / d[i] / d[i] + p[j] / d[j] / d[j] ) * difker_ij_y;
      //エネルギー積分値は初速度0のとき初期値0
      du[i]= du[i] + (p[i] / d[i] / d[i] )* m[j] * ((v[i][1] - v[j][1] ) * difker_ij_x + (v[i][2] - v[j][2] ) * difker_ij_y);

    }
  }
}
double Grad_h_term (int i, int k, double m[], double r[][DIM], double d, double h[])//kはベクトル成分
{
  double dens_h = 0;
  double dr = 0.0;
  int j = 0;
  for (j = 0; j < N_PTCL; j++){
    dr = r[i][k] - r[j][k];
    dens_h =dens_h + 2*m[j]*dr*dr*KernelFunc(i,j,r,h);
  }
  return h[i]*h[i]*d/dens_h;
}
void RungeKutta ( double dt, double m[], double r[][DIM], double v[][DIM], double a[][DIM], double p[], double d[], double u[], double du[], double h[] )
{
  double vp[N_ALL][DIM], up[N_ALL], a1[N_ALL][DIM], du1[N_ALL];
  int i = 0, j = 0, k = 0;
  double difker_i_x = 0.0 , difker_j_x = 0.0, difker_ij_x = 0.0, difker_i_y = 0.0 , difker_j_y = 0.0, difker_ij_y = 0.0 ,vis = 0.0;
  double fi = 0.0, fj =0.0;
  for(i=0; i < N_ALL; i++){
    //ステップ後の位置を計算 系内の粒子のみ
    r[i][1] = r[i][1] + v[i][1] * dt + a[i][1] * dt *dt / 2.0;    //ステップ後の位置はこの値となる
    if (r[i][1] > XMax){
      r[i][1] = r[i][1] - XMax;
    }else if (r[i][1] < XMin){
      r[i][1] = r[i][1] + XMax;
    }
    r[i][2] = r[i][2] + v[i][2] * dt + a[i][2] * dt *dt / 2.0;    //ステップ後の位置はこの値となる
    if (r[i][2] > YMax){
      r[i][2] = r[i][2] - YMax;
    }else if (r[i][2] < YMin){
      r[i][2] = r[i][2] + YMax;
    }
  }

  for( i = 0; i < N_ALL; i++){
    vp[i][1] = v[i][1] + a[i][1] *dt;
    vp[i][2] = v[i][2] + a[i][2] *dt;
    up[i] = u[i] + du[i] *dt;
  }

  //新しい位置とuで密度と圧力を計算　系内の粒子のみ、全粒子に対して
  for(i=0; i < N_ALL; i++){
    d[i] = 0;
    //smoothing lengthと密度の計算のため、境界近くにある粒子にも、周期的境界条件が適用されるように、各粒子を移動させる。
    Periodic_Boundary_Conditions(i,r);
    for(j=0; j < N_ALL; j++){
      d[i] = d[i] + m[j] * KernelFunc(i, j, r, h); //ステップ後の密度はこの値となる
    }
    p[i] = (Gamma - 1.0)* d[i] * up[i];
    h[i] = Smoothing_Length(i, m, r, d);
    //周期的境界条件のために移動させた位置をもとに戻す
    Return_Periodic_Boundary_Conditions(r);

  }
  //ステップ後の加速度とエネルギーの時間微分を計算　系内の粒子のみ、全粒子に対して
  for(i = 0; i < N_ALL; i++){
    if (r[i][1] > 0.75){  //i粒子が0.75よりも大きい位置にいる場合、ｘの負の位置にある粒子をx=1よりも大きい領域になければならない
      for( k = 0; k < N_ALL; k++){
        if (r[k][1] < 0.25){
          r[k][1] = r[k][1] + XMax;
        }
      }
    }else if (r[i][1] < 0.25){
      for(k = 0; k < N_ALL; k++){
        if ( r[k][1] > 0.75){
          r[k][1] = r[k][1] - XMax;
        }
      }
    }else if (r[i][2] > 0.75){
      for(k = 0; k < N_ALL; k++){
        if (r[k][2] < 0.25){
          r[k][2] = r[k][2] + YMax;
        }
      }
    }else if (r[i][2] < 0.25){
      for(k = 0; k < N_ALL; k++){
        if ( r[k][2] > 0.75){
          r[k][2] = r[k][2] - YMax;
        }
      }
    }

    a1[i][1] = 0.0;
    a1[i][2] = 0.0;
    du1[i] = 0.0;
    //fi = Grad_h_term(i, m, r, d[i], h);
    for(j = 0; j < N_ALL; j++){
      difker_i_x =  DifferKernelFunc(i, j, 1, r, h);
      difker_j_x =  DifferKernelFunc_hj(i, j, 1, r, h);
      difker_i_y =  DifferKernelFunc(i, j, 2, r, h);
      difker_j_y =  DifferKernelFunc_hj(i, j, 2, r, h);
      difker_ij_x = (difker_i_x + difker_j_x)/2.0;
      difker_ij_y = (difker_i_y + difker_j_y)/2.0;
      vis = ViscosityTerm(i, j, r, vp, d, up ,h);
      //fj = Grad_h_term(j, m, r, d[j], h);
      a1[i][1] = a1[i][1] + (-1.0) * m[j] * (p[i] / d[i] / d[i]  +  p[j] / d[j] / d[j])*difker_ij_x
              -1.0 * m[j] * vis * difker_ij_x;
      a1[i][2] = a1[i][2] + (-1.0) * m[j] * (p[i] / d[i] / d[i]  +  p[j] / d[j] / d[j])*difker_ij_y
              -1.0 * m[j] * vis * difker_ij_y;
      du1[i] = du1[i] + (p[i] / d[i] / d[i] ) * m[j] * ((vp[i][1] - vp[j][1] ) * difker_ij_x + (vp[i][2] - vp[j][2])*difker_ij_y)
              +1.0 * m[j] * vis *  ((vp[i][1] - vp[j][1] ) * difker_ij_x + (vp[i][2] - vp[j][2])*difker_ij_y)/ 2.0 ;
    }
    //周期的境界条件のために移動させた位置をもとに戻す
    for(k = 0; k < N_ALL; k++){
      if (r[k][1] > 1){
        r[k][1] = r[k][1] - XMax;
      }else if(r[k][1] < 0){
        r[k][1] = r[k][1] + XMax;
      }else if (r[k][2] > 1){
        r[k][2] = r[k][2] - YMax;
      }else if (r[k][2] < 0){
        r[k][2] = r[k][2] + YMax;
      }
    }
  }
  //ステップ後の速度と内部エネルギーの計算 系内の粒子のみ
  for(i = 0; i < N_PTCL; i++){
    v[i][1] = v[i][1] + (a[i][1] + a1[i][1]) * dt / 2.0;      //ステップ後の速度
    v[i][2] = v[i][2] + (a[i][2] + a1[i][2]) * dt / 2.0;      //ステップ後の速度
    u[i] = u[i] + (du[i]+ du1[i]) * dt /2.0;      //ステップ後のエネルギー
  }
  //加速度とエネルギーの時間微分の1ステップ後を本来の配列に代入
  for(i = 0; i < N_PTCL; i++){
    a[i][1] = a1[i][1];
    a[i][2] = a1[i][2];
    du[i] = du1[i];
    p[i] = (Gamma - 1.0)* d[i] * u[i];
  }

}
void PrintData    (FILE *file, double m[], double r[][DIM], double v[][DIM], double a[][DIM], double p[], double d[], double u[], double du[], double h[] )
{
  int i = 0;
  for(i = 0; i < N_PTCL; i++){
    if (i < N_PTCL - 1){
      fprintf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", m[i], r[i][1], r[i][2], v[i][1], v[i][2], a[i][1],a[i][2], p[i], d[i], u[i], du[i], h[i]);
    }else{
      fprintf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n\n", m[i], r[i][1], r[i][2], v[i][1], v[i][2], a[i][1],a[i][2], p[i], d[i], u[i], du[i], h[i]);
    }
  }
}
void Mesh_PrintData    (FILE *file, double r[Mesh+1][Mesh+1][DIM], double d[Mesh+1][Mesh+1])
{
  int i = 0, j = 0;
  for(i = 0; i < Mesh+1; i++){
    for(j = 0; j<Mesh+1;j++){
      if (j < Mesh){
        fprintf(file, "%lf, %lf, %lf \n", r[i][j][1], r[i][j][2], d[i][j]);
      }else{
        fprintf(file, "%lf, %lf, %lf \n\n", r[i][j][1], r[i][j][2], d[i][j]);
      }
    }
  }
}
//double用比較関数
int asc(const void *a, const void *b) {
  double *A = (double *)a;
  double *B = (double *)b;
  if (*A > *B) return 1;
  if (*A < *B) return -1;
  return 0;
}
//配列内の最大値を返す関数
double MaxArray( double n[], int l )
{
  int i = 0;
  double max;
  max = n[0];
  for(i = 0; i< l; i++){
    if(max < n[i]) max = n[i];
  }
  return max;
}
//配列内の最小値を返す関数
double MinArray( double n[], int l )
{
  int i = 0;
  double min;
  min = n[0];
  for(i = 0; i< l; i++){
    if(min > n[i]) min = n[i];
  }
  return min;
}
//2数を比較し、大きい方の数字を返す
double Compare_Max( double x, double y)
{
  if(x < y){
    return y;
  }else {
    return x;
  }
}

void Current_time()
{
  time_t now = time(NULL);
  struct tm *pnow = localtime(&now);
  char buff[128]="";
  sprintf(buff,"%d:%d:%d",pnow->tm_hour,pnow->tm_min,pnow->tm_sec);
  printf("%s\n",buff);
}
