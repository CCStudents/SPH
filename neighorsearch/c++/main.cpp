// 練習問題　その１　一様球の近傍粒子探査(C++ base)

// ライブラリ（Python の import from ... と似た作業）
#include "header.hpp"  // headerファイルをまとめた

// メイン関数の中身が実行ファイルで実行される
int main ( void ) // int型の返り値を持つ void mainでも良い
{
  // 粒子分布生成のために必要なParticle classを準備
  // // まずParticle systemのポインタを作成
  // // Heap領域にnewで粒子数分のメモリを確保する
  Particle * ptclSys;
  ptclSys = new Particle[N_PTCL];

  // 一様球の作成
  makeSphere(ptclSys);

  // *************************************
  // 課題： ここで各粒子の近くにある粒子を探査すること
  // *************************************


  // 粒子データの出力
  printData(ptclSys);
  return 0; // int の返り値として 0 を return する
}
