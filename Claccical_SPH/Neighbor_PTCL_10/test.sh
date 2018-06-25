#!/bin/bash
# プログラムのコンパイルとデータ作成を行う shell script
# ** $ bash test.sh で実行できる

#make clean
rm data/*
make

./sph.out

# gnuplot で位置情報を描画　q をコンソールに入力してグラフモード終了
gnuplot << EOF
  load "gnuplot_sample.gp"
EOF

mv -f *.ptcl data/
