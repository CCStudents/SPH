#!/bin/bash
# プログラムのコンパイルとデータ作成を行う shell script
# ** $ bash test.sh で実行できる

#make clean
rm data/*
rm *.png
make

./sph.out

mv -f *.ptcl data/

# gnuplot で位置情報を描画　q をコンソールに入力してグラフモード終了
gnuplot << EOF
  load "gnuplot_sample.gp"
EOF

#終了時に音を鳴らす
#paplay /usr/share/sounds/freedesktop/stereo/complete.oga
