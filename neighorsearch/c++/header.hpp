// ライブラリ（Python の import from ... と似た作業）
//// C versionの課題１から増やした
//// いらないパッケージも増えてるが，これくらいいれておけばだいたいのことはできる
#include <cmath>        // C言語の数学関数ライブラリ
#include <cstdio>       // C言語の標準ライブラリ
#include <cstdlib>      // C言語の標準ライブラリ
#include <unistd.h>     // getopt()を使うためのライブラリ
#include <iostream>     // 基本的なstream入出力機能
#include <sstream>      // stringstream利用機能
#include <fstream>      // filestream ファイルの入出力関数
#include <vector>       // vector neighbor listの作成
#include <string>       // 文字列クラス
#include <limits>       // min,maxなどの算術関数
#include <sys/stat.h>   // ファイルの状態を得る（フォルダ作成時のchmod用）

#define real double         // myvector.hのためにreal型はdoubleで考える
#include "myvector.hpp"     // vectorを自分で定義したヘッダーファイル
#include "constant.hpp"     // 定数定義(const type)
#include "class.hpp"        // c++の特色クラス定義 Particle
#include "distribution.hpp" // 粒子分布の作成
#include "io.hpp"           // input/output
