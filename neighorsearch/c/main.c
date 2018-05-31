// 練習問題　その１　一様球の近傍粒子探査(C base)
// ############################################################
// ライブラリ（Python の import from ... と似た作業）
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// ############################################################
// マクロ定義　この文字列は後ろの数字で置き換えられる　型は自動で指定される
#define N_PTCL    1024 // 全粒子数
#define DIM       3    // 次元数
#define M_SPHERE  1.0  // 一様球の質量
#define R_SPHERE  1.0  // 一様球の半径
#define r_VIRIAL  0.5  // ビリアル比
// ############################################################
// 関数定義　
// このプログラムのように main() の後に具体的な内容を書く場合もあるが、
// 基本的にはヘッダーファイルにまとめてしまえば良い
void gaussian   ( void );
void makeSphere ( double m[], double x[][DIM], double v[][DIM] );
void printData  ( double m[], double x[][DIM], double v[][DIM] );
// ############################################################
// メイン関数の中身が実行ファイルで実行される
int main ( void ) // int型の返り値を持つ void mainでも良い
{
  // 粒子分布生成のために必要な配列を準備(struct を使うこともできるが)
  // static は静的メモリ確保。　メモリを保持し続けることができる。
  static double mass[N_PTCL], pos[N_PTCL][DIM], vel[N_PTCL][DIM], acc[N_PTCL][DIM];

  // 一様球の作成
  makeSphere(mass, pos, vel);

  // *************************************
  // 課題： ここで各粒子の近くにある粒子を探査すること
  // *************************************


  // 粒子データの出力
  printData(mass, pos, vel);
  return 0; // int の返り値として 0 を return する
}
// ############################################################
// 関数の具体的な内容記述
double gaussian ( void )
{
  double x, y, r2;
  do{
    x  = 2.0 * drand48() - 1.0;
    y  = 2.0 * drand48() - 1.0;
    r2 = (x * x) + (y * y);
  }while((r2 >= 1.0) || (r2 == 0.0));
  // ガウシアンに従う確率分布で[0:1]の値を返す
  return sqrt(-2.0 * log(r2) / r2) * x;
}
void makeSphere ( double m[], double x[][DIM], double v[][DIM] )
{
  // C言語では変数の初期化を意識する
  int nCount = 0, i = 0, j = 0, k = 0;
  double r = 0.0; // 粒子の半径

  // 粒子の位置計算　モンテカルロ法
  while(nCount < N_PTCL){
    // 実行文が１文だけの場合は {} を省略できる
    // for(k = 0; k < DIM; k++) x[nCount][k] = 2.0 * drand48() - 1.0;

    // ループアンローリングという手法　for文の中身をそのまま書くほうが早かったりする
    x[nCount][0] = 2.0 * drand48() - 1.0;
    x[nCount][1] = 2.0 * drand48() - 1.0;
    x[nCount][2] = 2.0 * drand48() - 1.0;

    // r = sqrt(pow(x[nCount][0], 2.0) + pow(x[nCount][1], 2.0) + pow(x[nCount][2], 2.0));
    // 半径の計算　pow() は遅いので 単純な掛け算で書くことが多い
    r = sqrt( x[nCount][0] * x[nCount][0]
            + x[nCount][1] * x[nCount][1]
            + x[nCount][2] * x[nCount][2]);
    if(r < R_SPHERE) nCount++;
  }

  double potEnergy = 0.0;
  double dr[DIM];
  for(i = 0; i < N_PTCL - 1; i++){
    // 各粒子の質量計算　N_PTCL はおそらく int 型なので(double)でキャストする
    m[i] = (double)(M_SPHERE  / N_PTCL);
    for(j = i + 1; j < N_PTCL; j++){
      m[j] = (double)(M_SPHERE  / N_PTCL);
      // 相対距離の計算
      dr[0] = x[j][0] - x[i][0];
      dr[1] = x[j][1] - x[i][1];
      dr[2] = x[j][2] - x[i][2];
      // 相対距離の絶対値の計算
      r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
      potEnergy += -(m[i] * m[j]) / r;
    }
  }

  // 速度分散の計算
  double vSigma = sqrt((2.0 * r_VIRIAL * fabs(potEnergy)) / (3.0 * M_SPHERE));
  // 速度分散 vSigma のガウシアンの確率分布で速度を与える
  for(i = 0; i < N_PTCL; i++){
    v[i][0] = vSigma * gaussian();
    v[i][1] = vSigma * gaussian();
    v[i][2] = vSigma * gaussian();
  }
}
void printData  ( double m[], double x[][DIM], double v[][DIM] )
{
  // 粒子データの出力
  int i = 0;
  for(i = 0; i < N_PTCL; i++){
    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
            m[i], x[i][0], x[i][1], x[i][2], v[i][0], v[i][1], v[i][2]);
  }
}
