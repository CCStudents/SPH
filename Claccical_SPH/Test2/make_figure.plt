reset;
set size square;
set grid;
!find ./ -name "*.ptcl" | wc -l | awk '{print"n1="$1}'>tmp.gp  #.dataのファイル数を取得、その値をn1に代入するという命令をtmp.gpに書き込む
load "tmp.gp"                           # ループ変数の最大値
set term png


#速度のグラフ書き出し
set xr[-0.5:1.5];
set yr[-2.1:2.1];
set output sprintf("data_%d_velocity.png" , n1-1)
set title "velocity"
unset key
filename = sprintf("data_%d.ptcl", n1-1)     # n番目のデータファイルの名前の生成
plot filename using 2:3 with points pt 7 # 1:質量 2:位置 3:速度 4:加速度 5:圧力 6:密度 7:内部エネルギー 8:内部エネルギー時間微分 9:smoothing length

#密度のグラフ書き出し
set xr[-0.5:1.5];
set yr[-0.2:1.2];
set output sprintf("data_%d_density.png" , n1-1)
set title "density"
unset key
filename = sprintf("data_%d.ptcl", n1-1)     # n番目のデータファイルの名前の生成
plot filename using 2:6 with points pt 7 # 1:質量 2:位置 3:速度 4:加速度 5:圧力 6:密度 7:内部エネルギー 8:内部エネルギー時間微分 9:smoothing length

#圧力のグラフ書き出し
set xr[-0.5:1.5];
set yr[0:0.6];
set output sprintf("data_%d_pressure.png" , n1-1)
set title "pressure"
unset key
filename = sprintf("data_%d.ptcl", n1-1)     # n番目のデータファイルの名前の生成
plot filename using 2:5 with points pt 7 # 1:質量 2:位置 3:速度 4:加速度 5:圧力 6:密度 7:内部エネルギー 8:内部エネルギー時間微分 9:smoothing length

#内部エネルギーのグラフ書き出し
set xr[-0.5:1.5];
set yr[0:1.1];
set output sprintf("data_%d_energy.png" , n1-1)
set title "energy"
unset key
filename = sprintf("data_%d.ptcl", n1-1)     # n番目のデータファイルの名前の生成
plot filename using 2:7 with points pt 7 # 1:質量 2:位置 3:速度 4:加速度 5:圧力 6:密度 7:内部エネルギー 8:内部エネルギー時間微分 9:smoothing length
