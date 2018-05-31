# プログラムのコンパイルとデータ作成を行う shell script
# ** $ bash test.sh で実行できる
make clean
make
# 標準出力をdata.ptclに上書き保存する
./search.out > data.ptcl
# gnuplot で位置情報を描画　q をコンソールに入力してグラフモード終了
gnuplot -e "
  reset;
  set size square;
  set view equal xyz;
  set grid;
  set ticslevel 0;
  set xr[-2:2];
  set yr[-2:2];
  set zr[-2:2];
  spl 'data.ptcl' u 2:3:4 w dot title 'position';
  pause -1
"
