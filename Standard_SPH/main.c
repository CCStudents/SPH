#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"Test4.h"

#define N_PTCL        500        //系内の全粒粒子数
#define N_Satb_PTCL   100         //考えている系の前後にそれぞれ配置する粒子数
#define N_ALL   N_PTCL + N_Satb_PTCL
#define Gamma   1.4              //比熱比
#define StepN   5000             //ステップ数
#define TOL     0.000001          //p_starの収束条件

void Current_time();

double KernelFunc  ( int i,int j, double x[], double h[] );
double DifferKernelFunc  ( int i,int j, double x[], double h[] );
double DifferKernelFunc_hj ( int i, int j, double x[], double h[] );
double Grad_h_term (int i, double m[], double x[], double d, double h[]);


double ViscosityTerm ( int i, int j, double x[], double v[], double d[], double u[], double h[]);
double Smoothing_Length ( int i, double m[], double x[],double d[]);
double Time_Step   ( double x[], double v[], double p[], double d[], double u[], double h[]);

void InicialCondi  ( double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] );
void RungeKutta    ( double dt, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] );
void PrintData     ( FILE *file, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] );

double Star_Press (void);
double Star_Velocity(double p_star);
double Star_Dens  (char i, double p_star);
double Speed (char i, char j, double p_star, double u_star);
void Analytical_solution(double x[], double v[], double p[], double d[], double u[]);
double State_fan(int i, int j, char str, double x[]);    //j=1 密度、j=2 速度、 j=3 圧力
double Press_Two_Shock_Approximation(void);

double Func_shock ( double p_star, double press_k, double rho_k);
double Func_rarefac (double p_star, double press_k, double rho_k );
double Dif_Func_shock ( double p_star, double press_k, double rho_k);
double Dif_Func_rarefac (double p_star, double press_k, double rho_k );

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
  static double mass[N_ALL], pos[N_ALL], vel[N_ALL], acc[N_ALL], press[N_ALL],
                dens[N_ALL], energy[N_ALL], difenergy[N_ALL], len[N_ALL];
  static double vel_solution[N_ALL], press_solution[N_ALL], dens_solution[N_ALL], energy_solution[N_ALL];
  double timestep, totaltime = 0.0;
  printf("%f, %f, %f, %f\n",Star_Press(),Star_Velocity(Star_Press()),Star_Dens('L',Star_Press()), Star_Dens('R',Star_Press()));
  printf("%f, %f, %f, %f\n",Speed('L','H',Star_Press(),Star_Velocity(Star_Press())),
                            Speed('L','T',Star_Press(),Star_Velocity(Star_Press())),
                            Speed('R','H',Star_Press(),Star_Velocity(Star_Press())),
                            Speed('R','T',Star_Press(),Star_Velocity(Star_Press()))
                            );

  FILE *fp;
  int i = 0;
  for(i = 0; i < StepN; i++){

    if(i==0){
      InicialCondi(mass, pos, vel, acc, press, dens, energy, difenergy, len);
    }else{
      //timestepの計算
      timestep = Time_Step(pos, vel, press, dens, energy, len);
      totaltime = totaltime + timestep;
      if (totaltime > EndTime){
        totaltime = totaltime - timestep;
        timestep = EndTime - totaltime;
        RungeKutta(timestep, mass, pos, vel, acc, press, dens, energy, difenergy, len);
        //解析解の出力
        fp = fopen("analytical_data.ptcl", "w");
        Analytical_solution(pos, vel_solution, press_solution, dens_solution, energy_solution);
        PrintData(fp, pos, vel, vel_solution, press, press_solution, dens, dens_solution, energy, energy_solution);
        fclose(fp);
        end = clock();
        printf("%.2f秒\n",(double)(end-start)/CLOCKS_PER_SEC);
        return 0;
        exit(0);
      }
      RungeKutta(timestep, mass, pos, vel, acc, press, dens, energy, difenergy, len);
    }
  }
}

double Star_Press (void) //Star領域での圧力
{
  double rho_min, rho_max, p_min, p_max;
  double p_0 = 0.001;//(Press1 + Press2) / 2 ;
  double p_1 = 0;
  double function, dif_function;
  double du = Vel2 - Vel1;
  double CHA = 1; //収束条件の判定
  if( Press1 <= Press2 ){
    p_min = Press1;
    p_max = Press2;
    rho_min = Rho1;
    rho_max = Rho2;
  }else{
    p_min = Press2;
    rho_min = Rho2;
    p_max = Press1;
    rho_max = Rho1;
  }
  while(CHA > TOL){
    if(p_0 <= p_min){
      function = Func_rarefac(p_0, p_min,rho_min) + Func_rarefac(p_0, p_max, rho_max) +du;
      dif_function = Dif_Func_rarefac(p_0, p_min,rho_min) + Dif_Func_rarefac(p_0, p_max, rho_max);
    }else if( p_min < p_0 && p_0 < p_max){
      function = Func_shock(p_0, p_min,rho_min) + Func_rarefac(p_0, p_max, rho_max) +du;
      dif_function = Dif_Func_shock(p_0, p_min,rho_min) + Dif_Func_rarefac(p_0, p_max, rho_max);
    }else{
      function = Func_shock(p_0, p_min,rho_min) + Func_shock(p_0, p_max, rho_max) +du;
      dif_function = Dif_Func_shock(p_0, p_min,rho_min) + Dif_Func_shock(p_0, p_max, rho_max) ;
    }
    p_1 = p_0 - function / dif_function;
    CHA = 2 * fabs(p_1 - p_0)/ (p_1 + p_0);
    p_0 = p_1;
  }
  return p_1;
}
double Press_Two_Shock_Approximation(void)
{
  double a_L = sqrt(Gamma*Press1/Rho1);     //音速
  double a_R = sqrt(Gamma*Press2/Rho2);     //音速
  double p_0 = Compare_Max( TOL, (Press1+Press2)/2-(Vel2-Vel1)*(Rho1+Rho2)*(a_L+a_R)/8);
  double A_L = 2 / (Gamma + 1) / Rho1;
  double B_L = (Gamma - 1) / (Gamma + 1)* Press1;
  double A_R = 2 / (Gamma + 1) / Rho2;
  double B_R = (Gamma - 1) / (Gamma + 1)* Press2;
  double g_L = sqrt(A_L / (p_0 + B_L));
  double g_R = sqrt(A_R / (p_0 + B_R));
  double p_TS = (g_L * Press1 + g_R * Press2 -Vel2 + Vel1) / (g_L + g_R);
  return Compare_Max(TOL, p_TS);
}

double Star_Velocity(double p_star)  //u_star接触不連続面の特性線の計算
{
  double func;
  if( p_star > Press1 && p_star > Press2){
    func = Func_shock(p_star, Press2, Rho2) - Func_shock(p_star, Press1, Rho1);
  }else if( p_star > Press1 && p_star <= Press2 ){
    func = Func_rarefac(p_star, Press2, Rho2) - Func_shock(p_star, Press1, Rho1);
  }else if( p_star <= Press1 && p_star > Press2){
    func = Func_shock(p_star, Press2, Rho2) - Func_rarefac(p_star, Press1, Rho1);
  }else{
    func = Func_rarefac(p_star, Press2, Rho2) - Func_rarefac(p_star, Press1, Rho1);
  }
  return (Vel1 + Vel2 + func) / 2;
}

double Star_Dens (char i, double p_star)
{
  double g = (Gamma - 1) / (Gamma + 1);
  if(i == 'L' || i== 'l'){
    double p = Press1;
    if(p < p_star){                        //shock
      return Rho1 * ((g + (p_star / p) )/( (g * (p_star / p) + 1)));
    }else{                                //rarefaction
      return Rho1 * pow(p_star / p, 1 / Gamma);
    }
  }else if (i =='R'|| i == 'r'){
    double p = Press2;
    if(p < p_star){                        //shock
      return Rho2 * ((g + (p_star / p) )/ ((g * (p_star / p) + 1)));
    }else{                                //rarefaction
      return Rho2 * pow(p_star / p, 1 / Gamma);
    }
  }else{
    return 0.0;
  }
}

double Speed (char i, char j, double p_star, double u_star)//もし引数jにSが入っていたら、それはshock
{
  if(i == 'L' || i== 'l'){
    double p = Press1, u = Vel1;
    double a = sqrt(Gamma * p / Rho1);
    double a_star = sqrt(Gamma * p_star / Star_Dens('L', Star_Press()));
    if(p < p_star){                        //shock
      return u - a * sqrt((p_star*(Gamma+1)/p+(Gamma-1))/2/Gamma);
    }else{                                 //rarefaction
      if(j == 'H'){
        return u - a;
      }else if(j =='T'){
        return u_star - a_star;
      }else{
        return 0.0;
      }
    }
  }else if (i =='R'|| i == 'r'){
    double p = Press2, u = Vel2;
    double a = sqrt(Gamma * p / Rho2);
    double a_star = sqrt(Gamma * p_star /  Star_Dens('R', Star_Press()));
    if(p < p_star){                        //shock
      return u + a * sqrt((p_star*(Gamma+1)/p+(Gamma-1))/2/Gamma);
    }else{                                 //rarefaction
      if(j == 'H'){
        return u + a;
      }else if(j =='T'){
        return u_star + a_star;
      }else{
        return 0.0;
      }
    }
  }
}

double Func_shock ( double p_star, double press_k, double rho_k)
{
  double A_k = 2 / (Gamma +1) / rho_k;                //定数A
  double B_k = press_k * (Gamma - 1.0) / (Gamma + 1.0);   //定数B
  return (p_star - press_k) * sqrt(A_k / (p_star + B_k));
}

double Func_rarefac ( double p_star, double press_k, double rho_k)
{
  double a_k = sqrt(Gamma * press_k / rho_k);         //音速
  return 2 * a_k / (Gamma - 1.0) * (pow (p_star / press_k , (Gamma-1.0)/2/Gamma) - 1);
}
double Dif_Func_shock ( double p_star, double press_k, double rho_k)
{
  double A_k = 2 / (Gamma +1) / rho_k;                //定数A
  double B_k = press_k * (Gamma - 1.0) / (Gamma + 1.0);   //定数B
  return sqrt(A_k / (p_star + B_k))*(1 - (p_star - press_k)/ 2 / (B_k + p_star));
}

double Dif_Func_rarefac ( double p_star, double press_k, double rho_k)
{
  double a_k = sqrt(Gamma * press_k / rho_k);         //音速
  return (pow (p_star / press_k , -1*(Gamma+1)/2/Gamma))/ rho_k / a_k;
}

void Analytical_solution(double x[], double v[], double p[], double d[], double u[])
{
  double p_star = Star_Press();
  double u_star = Star_Velocity(Star_Press());
  int i = 0;
  for(i = 0; i < N_ALL; i++){
    if((x[i] - Xmid)<=u_star*EndTime){  //スター領域よりも左側にある粒子を計算
      //左側の領域
      if(Press1 < p_star){                      //shock
        if((x[i] - Xmid) < Speed('L','S',p_star,u_star)*EndTime){
          p[i] = Press1;
          v[i] = Vel1;
          d[i] = Rho1;
          u[i] = p[i] / (Gamma - 1) / d[i];
        }else if(Speed('L','S',p_star,u_star)*EndTime<=(x[i] - Xmid) && (x[i] - Xmid) <= u_star*EndTime){
          p[i] = p_star;
          v[i] = u_star;
          d[i] = Star_Dens('L',Star_Press());
          u[i] = p[i] / (Gamma - 1) / d[i];
        }
      }else{                                //rarefaction
        if((x[i] - Xmid) < Speed('L','H',p_star,u_star)*EndTime){
          p[i] = Press1;
          v[i] = Vel1;
          d[i] = Rho1;
          u[i] = p[i] / (Gamma - 1) / d[i];
        }else if(Speed('L','H',p_star,u_star)*EndTime <= (x[i] - Xmid) && (x[i] - Xmid) <= Speed('L','T',p_star,u_star)*EndTime){
          p[i] = State_fan(i,3,'L',x);
          v[i] = State_fan(i,2,'L',x);
          d[i] = State_fan(i,1,'L',x);
          u[i] = p[i] / (Gamma - 1) / d[i];
        }else if (Speed('L','T',p_star,u_star)*EndTime <=(x[i] - Xmid) ){
          p[i] = p_star;
          v[i] = u_star;
          d[i] = Star_Dens('L',Star_Press());
          u[i] = p[i] / (Gamma - 1) / d[i];
        }
      }
    }else if(u_star*EndTime <=(x[i] - Xmid)){ //左側の固定粒子をを除く
      //右側の領域
      if(Press2 < p_star){                      //shock
        if((x[i] - Xmid) >  Speed('R','S',p_star,u_star)*EndTime){
          p[i] = Press2;
          v[i] = Vel2;
          d[i] = Rho2;
          u[i] = p[i] / (Gamma - 1) / d[i];
        }else if(Speed('R','S',p_star,u_star)*EndTime>=(x[i] - Xmid) && (x[i] - Xmid) >= u_star*EndTime){
          p[i] = p_star;
          v[i] = u_star;
          d[i] = Star_Dens('R',Star_Press());
          u[i] = p[i] / (Gamma - 1) / d[i];
        }
      }else{                                //rarefaction
        if(u_star*EndTime <= (x[i] - Xmid)&&(x[i] - Xmid)<= Speed('R','T',p_star,u_star)*EndTime){
          p[i] = p_star;
          v[i] = u_star;
          d[i] = Rho2 * pow(p[i]/Press2 , 1/Gamma);
          u[i] = p[i] / (Gamma - 1) / d[i];
        }else if(Speed('R','H',p_star,u_star)*EndTime >= (x[i] - Xmid) && (x[i] - Xmid) >= Speed('R','T',p_star,u_star)*EndTime){
          p[i] = State_fan(i,3,'R',x);
          v[i] = State_fan(i,2,'R',x);
          d[i] = State_fan(i,1,'R',x);
          u[i] = p[i] / (Gamma - 1) / d[i];
        }else if((x[i] - Xmid) >= Speed('R','H',p_star,u_star)*EndTime){
          p[i] = Press2;
          v[i] = Vel2;
          d[i] = Rho2;
          u[i] = p[i] / (Gamma - 1) / d[i];
        }
      }
    }
  }
}
//rarefactionにおける位置・時間依存の物理状態
double State_fan(int i, int j, char str, double x[])    //j=1 密度、j=2 速度、 j=3 圧力
{
  double u, d, p, a, g;
  if(str == 'L' || str == 'l'){
    u = Vel1;
    d = Rho1;
    p = Press1;
    a = sqrt(Gamma * p / d);
    g = 2/(Gamma+1) +(Gamma-1)/(Gamma+1)/a*(u-(x[i]-Xmid)/EndTime);
    if(j==1){
      return d * pow(g,2 / (Gamma-1));
    }else if(j ==2){
      return 2 / (Gamma + 1)* (a + (Gamma -1) / 2 *u + (x[i]-Xmid)/EndTime);
    }else if (j==3){
      return p * pow(g, 2 * Gamma / (Gamma-1));
    }
  }else if(str =='R' || str == 'r'){
    u = Vel2;
    d = Rho2;
    p = Press2;
    a = sqrt(Gamma * p / d);
    g = 2 / (Gamma + 1) - (Gamma - 1) / (Gamma + 1) / a * (u- (x[i]-Xmid) / EndTime);
    if(j==1){
      return d * pow(g, 2 / (Gamma-1));
    }else if(j ==2){
      return 2 / (Gamma + 1) * (- 1 * a + (Gamma -1) / 2 *u + (x[i]-Xmid)/EndTime);
    }else if (j==3){
      return p * pow(g, 2 * Gamma / (Gamma-1));
    }
  }
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
//カーネル関数の微分∇Wijの定義
double DifferKernelFunc_hj ( int i, int j, double x[], double h[] )
{
  double dx = 0.0;
  dx = x[i] - x[j];
  return -2*dx*exp(-1*dx*dx/h[j]/h[j]) /h[j]/h[j]/h[j]/sqrt(M_PI);
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
double Smoothing_Length ( int i, double m[], double x[], double d[])
{
  /*
  int j = 0;
  double dx[N_ALL];
  for( j = 0; j < N_ALL; j++){
    dx[j] = fabs(x[i] - x[j]);
  }
  qsort(dx, N_ALL, sizeof(double), asc);
  return dx[Neighbor_PTCL];
  */
  return Neighbor_PTCL * m[i] / d[i];
}
//最も小さいdtの計算
double Time_Step ( double x[], double v[], double p[], double d[], double u[], double h[])
{
  double t[N_PTCL], mu[N_PTCL];
  int i = 0, j = 0;
  double dx = 0.0, dv = 0.0, c = 0.0;

  for ( i = 0; i < N_PTCL; i++){
    c = sqrt((Gamma-1) * u[i]);
    for( j = 0; j < N_PTCL; j++){
      dx = x[i] - x[j];
      dv = v[i] - v[j];
      mu[j] = h[i] * dv * dx / (dx * dx + Epsilon * h[i] * h[i]);
    }
    t[i] = h[i] / ( c + 0.6 * (Alpha * c + Beta * MaxArray(mu, (int) sizeof(mu) / sizeof(mu[0])) ));
  }
  return  0.25 * MinArray(t, (int) sizeof(t) / sizeof(t[0]));
}
//初期条件(密度、位置、速さ、圧力、加速度)の定義：密度を1,0.25と仮定→hを固定
void InicialCondi  ( double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] )
{
  int i = 0, j = 0;
  double N1 = N_PTCL * Rho1 / (Rho1 + Rho2);
  double N2 = N_PTCL * Rho2 / (Rho1 + Rho2);
  double difker_i = 0.0 , difker_j = 0.0, difker_ij = 0.0 ,vis = 0.0;
  for(i = 0; i < N_ALL; i++){
    m[i] = (Rho1 + Rho2) * DeltaX / N_PTCL;
    //密度を1と0.25にするために位置を調節
    //ifの中は0-0.5の範囲（粒子数多め）
    if( i <= N1 ){
      if( i == 0 ){
        x[i] = XMin;
      }else{
        x[i] = x[i-1] + DeltaX / N1;
      }
      p[i] = Press1;
      d[i] = Rho1;
      v[i] = Vel1;
    }else if ( N1 < i && i <= N_PTCL + N_Satb_PTCL / 2.0 ){
      x[i] = x[i-1] + DeltaX / N2;
      p[i] = Press2;
      d[i] = Rho2;
      v[i] = Vel2;
    }else{
      if ( i == N_PTCL + N_Satb_PTCL / 2.0 + 1.0 ){
        x[i] = x[0] - DeltaX / N1;
      }else{
        x[i] = x[i-1] - DeltaX / N1;
      }
      p[i] = Press1;
      d[i] = Rho1;
      v[i] = Vel1;
    }
    //初期化
    a[i] = 0.0;
    du[i] = 0.0;
  }
  //smoothing length の計算
  for(i = 0; i < N_ALL; i++){
    h[i] = Smoothing_Length(i, m, x, d);
  }
  //内部エネルギーの初期値
  for(i = 0; i < N_ALL; i++){
    u[i] = p[i] / (Gamma -1.0) / d[i];
  }

  //加速度とエネルギーの時間微分の計算
  for(i = 0; i < N_ALL; i++){
    for(j = 0; j < N_ALL; j++){
      difker_i =  DifferKernelFunc(i, j, x, h);
      difker_j =  DifferKernelFunc_hj(i, j, x, h);
      difker_ij = (difker_i+difker_j)/2.0;
      a[i] = a[i] + (-1.0) * m[j] * (p[i] / d[i] / d[i] + p[j] / d[j] / d[j] ) * difker_ij;
      //エネルギー積分値は初速度0のとき初期値0
      du[i]= du[i] + (p[i] / d[i] / d[i] )* m[j] * (v[i] - v[j] ) * difker_ij;
    }
  }
}

double Grad_h_term (int i, double m[], double x[], double d, double h[])
{
  double dens_h = 0;
  double dx = 0.0;
  int j = 0;
  for (j = 0; j < N_ALL; j++){
    dx = x[i]-x[j];
    dens_h =dens_h + 2*m[j]*dx*dx*KernelFunc(i,j,x,h);
  }
  return h[i]*h[i]*d/dens_h;
}



void RungeKutta ( double dt, double m[], double x[], double v[], double a[], double p[], double d[], double u[], double du[], double h[] )
{
  double vp[N_ALL], up[N_ALL], a1[N_ALL], du1[N_ALL];
  int i = 0, j = 0;
  double difker_i = 0.0 , difker_j = 0.0, difker_ij = 0.0, vis = 0.0;
  double fi = 0.0, fj =0.0;
  for(i=0; i < N_ALL; i++){
    //ステップ後の位置を計算 系内の粒子のみ
    x[i] = x[i] + v[i] * dt + a[i] * dt *dt / 2.0;    //ステップ後の位置はこの値となる
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
    p[i] = (Gamma - 1.0)* d[i] * up[i];
    h[i] = Smoothing_Length(i, m, x, d);
  }
  //ステップ後の加速度とエネルギーの時間微分を計算　系内の粒子のみ、全粒子に対して
  for(i = 0; i < N_PTCL; i++){
    a1[i] = 0.0;
    du1[i] = 0.0;
    fi = Grad_h_term(i, m, x, d[i], h);
    for(j = 0; j < N_ALL; j++){
      difker_i =  DifferKernelFunc(i, j, x, h);
      difker_j =  DifferKernelFunc_hj(i, j, x, h);
      difker_ij = (difker_i+difker_j)/2.0;
      vis = ViscosityTerm(i, j, x, vp, d, up ,h);
      fj = Grad_h_term(j, m, x, d[j], h);
      a1[i] = a1[i] + (-1.0) * m[j] * (fi*p[i] / d[i] / d[i] *difker_i + fj* p[j] / d[j] / d[j]*difker_j)
              -1.0 * m[j] * vis * difker_ij;
      du1[i] = du1[i] + fi * (p[i] / d[i] / d[i] ) * m[j] * (vp[i] - vp[j] ) * difker_i
              +1.0 * m[j] * vis * (vp[i] - vp[j] ) * difker_ij / 2.0 ;
    }
  }
  //ステップ後の速度と内部エネルギーの計算 系内の粒子のみ
  for(i = 0; i < N_PTCL; i++){
    v[i] = v[i] + (a[i] + a1[i]) * dt / 2.0;      //ステップ後の速度
    u[i] = u[i] + (du[i]+ du1[i]) * dt /2.0;      //ステップ後のエネルギー
  }
  //加速度とエネルギーの時間微分の1ステップ後を本来の配列に代入
  for(i = 0; i < N_PTCL; i++){
    a[i] = a1[i];
    du[i] = du1[i];
    p[i] = (Gamma - 1.0)* d[i] * u[i];
  }
  printf("%f\n",p[200] );
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
