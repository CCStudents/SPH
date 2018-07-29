if(exist("n")==0 || n<0) n = n0             # ループ変数の初期化
filename = sprintf("data_%d.ptcl", n-1)     # n番目のデータファイルの名前の読み取り
filetitle = sprintf("n = %d", n-1)
plot filename using 2:3 title filetitle with points pt 7 # 1:質量 2:位置x 3:位置y 4:速度x 5:速度y 6:加速度x
n = n + dn            # ループ変数の増加
if ( n < n1 ) reread  # ループの評価
undefine n            # ループ変数の削除
