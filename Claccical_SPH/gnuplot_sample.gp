reset;
set size square;
set grid;
set xr[-0.2:1.2];
set yr[0:1.2]
set term gif animate                    # 出力をgifアニメに設定
set output "animation.gif"              # 出力ファイル名の設定
n0 = 1                                  # ループ変数の初期値
n1 = 150                                # ループ変数の最大値
dn = 1                                  # ループ変数の増加間隔
load "make_animation.plt"

set term png
load "make_figure.plt"
