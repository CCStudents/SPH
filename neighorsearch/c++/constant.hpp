// 一回しかincludeしないようにするコマンド　２回以上ヘッダーファイルを読み込んではいけない
#pragma once
// ############################################################
// C++ではdefineよりも，const typeで定義すると良い　型も自分で指定できる
const int  N_PTCL   = 1024; // 全粒子数
const int  DIM      = 3;    // 次元数
const real M_SPHERE = 1.0;  // 一様球の質量
const real R_SPHERE = 1.0;  // 一様球の半径
const real r_VIRIAL = 0.5;  // ビリアル比
