#!/bin/bash
# プログラムのコンパイルとデータ作成を行う shell script
# ** $ bash test.sh で実行できる

#make clean
rm data/*
rm *.png
make

./sph.out
#ファイル内をソートして、gnuplotにおける線描画をきれいにする
sort -t , -k 1 -n analytical_data.ptcl -o analytical_data.ptcl
mv -f *.ptcl data/

# gnuplot で位置情報を描画　q をコンソールに入力してグラフモード終了
gnuplot << EOF
  load "gnuplot_sample.gp"
EOF

paplay /usr/share/sounds/freedesktop/stereo/complete.oga
