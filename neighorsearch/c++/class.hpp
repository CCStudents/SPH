#pragma once

class Particle{
private: // Particle class内でしかアクセスできない
  vector pos;   // 位置
  vector vel;   // 速度
  double mass;  // 質量
public: // どこからでもアクセス（読み書き・呼び出し）できる
  Particle () // コンストラクタ
  {
    // Particleを作った時に自動で呼び出される初期化関数と思えば良い
    pos = vel = 0.0;
    mass = 0.0;
  }
  // Particleの要素（メンバ変数）を得るための関数（メンバ関数）
  vector getPos  ( void ) {return pos;}
  vector getVel  ( void ) {return vel;}
  double getMass ( void ) {return mass;}

  // Particleの要素に値を入れるための関数
  void setPos  ( vector newVec )  {pos  = newVec;}
  void setVel  ( vector newVec )  {vel  = newVec;}
  void setMass ( double newMass ) {mass = newMass;}
};
