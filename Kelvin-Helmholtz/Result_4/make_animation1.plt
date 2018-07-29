set term gif animate                    # 出力をgifアニメに設定
set output "animation.gif"              # 出力ファイル名の設定
n0 = 1                                  # ループ変数の初期値
!find ./ -name "*.ptcl" | wc -l | awk '{print"n1="$1}'>tmp.gp  #.dataのファイル数を取得、その値をn1に代入するという命令をtmp.gpに書き込む
load "tmp.gp"                           # ループ変数の最大値
dn = 1                                  # ループ変数の増加間隔
set xr[0:1]
set yr[0:1.5]
load "../make_animation2.plt"
