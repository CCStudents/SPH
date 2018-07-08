reset;
set size square;
set grid;
set term png
#set xr[0:1];
set autoscale
filename = "analytical_data.ptcl"

# 1:位置 2:速度 3:速度の解 4:圧力 5:圧力の解 6:密度 7:密度の解 8:内部エネルギー 9:内部エネルギーの解
#速度のグラフ書き出し
set output "analytical_data_velocity.png"
set title "velocity"
unset key
plot filename using 1:2 with points pt 7,filename using 1:3 with lines lw 2

#圧力のグラフ書き出し
set output "analytical_data_pressure.png"
set title "pressure"
unset key
plot filename using 1:4 with points pt 7,filename using 1:5 with lines lw 2


#密度のグラフ書き出し
set output "analytical_data_density.png"
set title "density"
unset key
plot filename using 1:6 with points pt 7,filename using 1:7 with lines lw 2

#内部エネルギーのグラフ書き出し
set output "analytical_data_energy.png"
set title "energy"
unset key
plot filename using 1:8 with points pt 7,filename using 1:9 with lines lw 2
