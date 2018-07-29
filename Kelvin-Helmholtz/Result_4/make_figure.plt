reset;
set size square;
set grid;
set term png
set autoscale
#set xr[0.3:1.2];
load "tmp.gp"                           # ループ変数の最大値
filename = sprintf("data_%d.ptcl", n1-1)     # n番目のデータファイルの名前の読み取り

# 1:質量 2:位置x 3:位置y 4:速度x 5:速度y 6:加速度x
#速度のグラフ書き出し
set output "analytical_data_velocity.png"
set title "velocity"
unset key
plot filename using 1:2 with points pt 7

#圧力のグラフ書き出し
set output "analytical_data_pressure.png"
set title "pressure"
unset key
plot filename using 1:4 with points pt 7


#密度のグラフ書き出し
set output "analytical_data_density.png"
set title "density"
unset key
plot filename using 1:6 with points pt 7

#内部エネルギーのグラフ書き出し
set output "analytical_data_energy.png"
set title "energy"
unset key
plot filename using 1:8 with points pt 7
