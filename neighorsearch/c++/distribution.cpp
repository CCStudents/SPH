#include "header.hpp"

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
void makeSphere ( Particle * pFirst )
{
  Particle * pTmp = pFirst, * pTmpI, * pTmpJ;
  // C言語では変数の初期化を意識する
  int nCount = 0;
  double r = 0.0; // 粒子の半径

  // 粒子の位置計算　モンテカルロ法
  vector x = 0.0;
  while(nCount < N_PTCL){
    // 実行文が１文だけの場合は {} を省略できる
    // for(k = 0; k < DIM; k++) x[k] = 2.0 * drand48() - 1.0;

    // ループアンローリングという手法　for文の中身をそのまま書くほうが早かったりする
    x[0] = 2.0 * drand48() - 1.0;
    x[1] = 2.0 * drand48() - 1.0;
    x[2] = 2.0 * drand48() - 1.0;

    // 半径の計算
    r = sqrt(x * x);
    if(r < R_SPHERE){
      pTmpI = (pTmp + nCount);
      pTmpI->setPos(x);
      nCount++;
    }
  }

  double potEnergy = 0.0;
  vector dr = 0.0;
  for(int i = 0; i < N_PTCL - 1; i++){
    // 各粒子の質量計算　N_PTCL は const int 型なので(double)でキャストする
    pTmpI = (pTmp + i);
    pTmpI->setMass((double)(M_SPHERE  / N_PTCL));
    for(int j = i + 1; j < N_PTCL; j++){
      pTmpJ = (pTmp + j);
      pTmpJ->setMass((double)(M_SPHERE  / N_PTCL));
      // 相対距離の計算
      dr = pTmpJ->getPos() - pTmpI->getPos();
      // 相対距離の絶対値の計算
      r = sqrt(dr * dr);
      potEnergy += -(pTmpJ->getMass() * pTmpI->getMass()) / r;
    }
  }

  // 速度分散の計算
  double vSigma = sqrt((2.0 * r_VIRIAL * fabs(potEnergy)) / (3.0 * M_SPHERE));
  vector v = 0.0;
  // 速度分散 vSigma のガウシアンの確率分布で速度を与える
  for(int i = 0; i < N_PTCL; i++){
    v[0] = vSigma * gaussian();
    v[1] = vSigma * gaussian();
    v[2] = vSigma * gaussian();
    pTmpI = (pTmp + i);
    pTmpI->setVel(v);
  }
}
