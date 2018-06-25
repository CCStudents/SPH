#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N_PTCL        500       //系内の全粒粒子数
#define N_Satb_PTCL   50        //考えている系の前後にそれぞれ配置する粒子数
#define Neighbor_PTCL 4        //近傍粒子数
#define N_ALL   N_PTCL + N_Satb_PTCL
#define Gamma   1.4         //比熱比
#define StepN   7000         //ステップ数
#define Rho1    5.99924         //初期密度1
#define Rho2    5.99242        //初期密度2
#define Vel1    19.5975         //初期速度1
#define Vel2    -6.19633         //初期速度2
#define Press1  460.894         //初期圧力1
#define Press2  46.0950      //初期圧力2
#define Alpha   1.0         //粘性項のα
#define Beta    2.0         //粘性項のβ
#define Epsilon 0.01        //粘性項のε


double KernelFunc  ( int i,int j, double x[], double h[] );
double DifferKernelFunc  ( int i,int j, double x[], double h[] );
double ViscosityTerm ( int i, int j, double x[], double v[], double d[], double u[], double h[]);
double Smoothing_Length ( int i, double x[]);
double Time_Step   ( double x[], double v[], double p[], double d[],double h[]);

void InicialCondi  ( double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] );
void RungeKutta    ( double dt, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] );
void PrintData     ( FILE *file, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] );

//double用比較関数
int asc(const void *a, const void *b) ;

//配列内の最大値・最小値を返す関数
double MaxArray( double n[], int l );
double MinArray( double n[], int l );

int main (void)
{
  //質量、位置、速度、加速度、圧力、密度、内部エネルギー、内部エネルギーの時間変化、影響半径
  static double mass[N_ALL], pos[N_ALL], vel[N_ALL], acc[N_ALL], press[N_ALL],
                dens[N_ALL], energy[N_ALL], difenergy[N_ALL], len[N_ALL];
  double timestep, totaltime = 0.0;
  ///*
  FILE *fp;
  int i = 0;
  char filename[100];
  for(i = 0; i < StepN; i++){
    sprintf(filename, "data_%d.ptcl",i);
    fp = fopen(filename, "w");
    if( fp == NULL ) {
      perror("ファイルの読み込みに失敗！\n");
      return 1;
    }

    if(i==0){
      InicialCondi(mass, pos, vel, acc, press, dens, energy, difenergy, len);
      PrintData(fp, mass, pos, vel, acc, press, dens, energy, difenergy, len);
    }else{
      //timestepの計算
      timestep = Time_Step(pos, vel, press, dens, len);
      totaltime = totaltime + timestep;
      if (totaltime > 0.25){
        timestep = totaltime - 0.25;
        RungeKutta(timestep, mass, pos, vel, acc, press, dens, energy, difenergy, len);
        PrintData(fp, mass, pos, vel, acc, press, dens, energy, difenergy, len);
        printf("%f\n",timestep);
        fclose(fp);
        exit(0);
      }
      RungeKutta(timestep, mass, pos, vel, acc, press, dens, energy, difenergy, len);
      PrintData(fp, mass, pos, vel, acc, press, dens, energy, difenergy, len);
      printf(" %f\n",totaltime);
    }
    fclose(fp);
  }
  return 0;
}

//カーネル関数Wijの定義
double KernelFunc ( int i, int j, double x[], double h[] )
{
  double dx = 0.0;
  dx = x[i] - x[j];
  return exp(-1*dx*dx/h[i]/h[i]) /h[i]/sqrt(M_PI);
}

//カーネル関数の微分∇Wijの定義
double DifferKernelFunc ( int i, int j, double x[], double h[] )
{
  double dx = 0.0;
  dx = x[i] - x[j];
  return -2*dx*exp(-1*dx*dx/h[i]/h[i]) /h[i]/h[i]/h[i]/sqrt(M_PI);
}

//人工粘性項の定義
double ViscosityTerm ( int i, int j, double x[], double v[], double d[], double u[], double h[])
{
  double   e = 0.01;
  double dx = 0.0, dv = 0.0, c = 0.0, mu = 0.0;
  dx = x[i] - x[j];
  dv = v[i] - v[j];
  //音速の計算
  c = (sqrt((Gamma - 1) * u[i]) + sqrt((Gamma-1) * u[j])) / 2;
  mu = h[i]*dv*dx /(dx * dx + Epsilon * h[i] * h[i]);
  if ( dx * dv < 0){
    return 2*(-1 * Alpha * c * mu + Beta * mu * mu) / (d[i] + d[j]);
  }else{
    return 0.0;
  }
}
//Smoothing Length の計算　i粒子に対する近傍粒子数に応じて長さを変える
double Smoothing_Length ( int i, double x[] )
{
  int j = 0;
  double dx[N_ALL];
  for( j = 0; j < N_ALL; j++){
    dx[j] = fabs(x[i] - x[j]);
  }
  qsort(dx, N_ALL, sizeof(double), asc);
  return dx[Neighbor_PTCL];
}


//最も小さいdtの計算
double Time_Step ( double x[], double v[], double p[], double d[],double h[])
{
  double t[N_PTCL], mu[N_PTCL];
  int i = 0, j = 0;
  double dx = 0.0, dv = 0.0, c = 0.0;

  for ( i = 0; i < N_PTCL; i++){
    c = sqrt(Gamma * p[i] / d[i]);
    for( j = 0; j < N_PTCL; j++){
      dx = x[i] - x[j];
      dv = v[i] - v[j];
      mu[j] = h[i] * dv * dx / (dx * dx + Epsilon * h[i] * h[i]);
    }
    t[i] = h[i] / ( c + 0.6 * (Alpha * c + Beta * MaxArray(mu, (int) sizeof(mu) / sizeof(mu[0])) ));
  }
  return  0.5 * MinArray(t, (int) sizeof(t) / sizeof(t[0]));
}


//初期条件(密度、位置、速さ、圧力、加速度)の定義：密度を1,0.25と仮定→hを固定
void InicialCondi  ( double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] )
{
  int i = 0, j = 0;
  double N1 = N_PTCL * Rho1 / (Rho1 + Rho2);
  double N2 = N_PTCL * Rho2 / (Rho1 + Rho2);
  for(i = 0; i < N_ALL; i++){
    m[i] = (Rho1 + Rho2) / 2 / N_PTCL;
    //密度を1と0.25にするために位置を調節
    //ifの中は0-0.5の範囲（粒子数多め）
    if( i <= N1 ){
      if( i == 0 ){
        x[i] = 0;
      }else{
        x[i] = x[i-1] + 0.5 / N1;
      }
      p[i] = Press1;
      d[i] = Rho1;
      v[i] = Vel1;
      //h[i] = 2*m[i] / Rho1;
    }else if ( N1 < i && i <= N_PTCL + N_Satb_PTCL / 2 ){
      x[i] = x[i-1] + 0.5 / N2;
      p[i] = Press2;
      d[i] = Rho2;
      v[i] = Vel2;
      //h[i] = 2*m[i] / Rho2;
    }else{
      if ( i == N_PTCL + N_Satb_PTCL / 2 + 1 ){
        x[i] = x[0] - 0.5 / N1;
      }else{
        x[i] = x[i-1] - 0.5 / N1;
      }
      p[i] = Press1;
      d[i] = Rho1;
      v[i] = Vel1;
      //h[i] = 2* m[i] / Rho1;
    }
    //初期化
    a[i] = 0.0;
    du[i] = 0.0;
  }
  //smoothing length の計算
  for(i = 0; i < N_ALL; i++){
    h[i] = Smoothing_Length(i, x);
  }
  //現在の粒子位置における密度を計算
  for(i = 0; i < N_PTCL; i++){
    d[i] = 0.0;
    for(j = 0; j < N_ALL; j++){
      d[i] = d[i] + m[j] * KernelFunc(i, j, x, h);
    }
  }



  //内部エネルギーの初期値
  for(i = 0; i < N_ALL; i++){
    u[i] = p[i] / (Gamma -1) / d[i];
  }

  //加速度とエネルギーの時間微分の計算
  for(i = 0; i < N_ALL; i++){
    for(j = 0; j < N_ALL; j++){
      a[i] = a[i] + (-1) * m[j] * (p[i] / d[i] / d[i] + p[j] / d[j] / d[j] ) * DifferKernelFunc(i, j, x, h);
      //エネルギー積分値は初速度0なので初期値0
      du[i]= du[i] + (p[i] / d[i] / d[j] )* m[j] * (v[i] - v[j] ) * DifferKernelFunc(i, j, x, h);
    }
  }
}

void RungeKutta ( double dt, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] )
{
  double vp[N_ALL], up[N_ALL], a1[N_ALL], du1[N_ALL];
  int i = 0, j = 0, k = 0, l = 0, n = 0, o = 0, q = 0;

  for(i=0; i < N_PTCL; i++){
    //ステップ後の位置を計算 系内の粒子のみ
    x[i] = x[i] + v[i] * dt + a[i] * dt *dt / 2;    //ステップ後の位置はこの値となる
    //境界に行った場合、速度を逆向きにする
    if( x[i] < 0 || x[i] > 1){
      v[i] = -1 * v[i];
    }
  }

  for( i = 0; i < N_ALL; i++){
    vp[i] = v[i] + a[i] *dt;
    up[i] = u[i] + du[i] *dt;
  }

  //新しい位置とuで密度と圧力を計算　系内の粒子のみ、全粒子に対して
  for(i=0; i < N_PTCL; i++){
    d[i] = 0;
    for(j=0; j < N_ALL; j++){
      d[i] = d[i] + m[j] * KernelFunc(i, j, x, h); //ステップ後の密度はこの値となる
    }
    p[i] = (Gamma - 1)* d[i] * up[i];
  }
  //ステップ後の加速度とエネルギーの時間微分を計算　系内の粒子のみ、全粒子に対して
  for(i = 0; i < N_PTCL; i++){
    a1[i] = 0;
    du1[i] = 0;
    for(j = 0; j < N_ALL; j++){
      a1[i] = a1[i] + (-1) * m[j] * (p[i] / d[i] / d[i] + p[j] / d[j] / d[j] ) * DifferKernelFunc(i, j, x, h)
              -1 * m[j] * ViscosityTerm(i, j, x, vp, d, up ,h) * DifferKernelFunc(i, j, x, h) ;
      du1[i] = du1[i] + (p[i] / d[i] / d[i] ) * m[j] * (vp[i] - vp[j] ) * DifferKernelFunc(i, j, x, h)
              +1 * m[j] * ViscosityTerm(i, j, x, vp, d, up, h) * (vp[i] - vp[j] ) * DifferKernelFunc(i, j, x, h) / 2 ;
    }
  }
  //ステップ後の速度と内部エネルギーの計算 系内の粒子のみ
  for(i = 0; i < N_PTCL; i++){
    v[i] = v[i] + (a[i] + a1[i]) * dt / 2;      //ステップ後の速度
    u[i] = u[i] + (du[i]+ du1[i]) * dt /2;      //ステップ後のエネルギー
  }
  //加速度とエネルギーの時間微分の1ステップ後を本来の配列に代入 影響半径の計算
  for(i = 0; i < N_PTCL; i++){
    a[i] = a1[i];
    du[i] = du1[i];
    h[i] = Smoothing_Length(i, x);
  }
}

void PrintData    (FILE *file, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] )
{
  int i = 0;
  for(i = 0; i < N_ALL; i++){
    fprintf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", m[i], x[i], v[i], a[i], p[i], d[i], u[i], du[i], h[i]);
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
