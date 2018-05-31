#pragma once

class Particle{
private: // Particle class内でしかアクセスできない
  vector pos;   // 位置
  vector vel;   // 速度
  vector acc;   // 加速度
  double mass;  // 質量
public: // どこからでもアクセス（読み書き・呼び出し）できる
  Particle () // コンストラクタ
  {
    // Particleを作った時に自動で呼び出される初期化関数と思えば良い
    pos = vel = acc = 0.0;
    mass = 0.0;
  }
  // Particleの要素（メンバ変数）を得るための関数（メンバ関数）
  void getPos  ( void ) {return pos;}
  void getVel  ( void ) {return vel;}
  void getAcc  ( void ) {return acc;}
  void getMass ( void ) {return mass;}
};
