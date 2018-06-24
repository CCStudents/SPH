if(exist("n")==0 || n<0) n = n0             # ループ変数の初期化
set output sprintf("data_%d.png" , n-1)     #output先の指定
filename = sprintf("data_%d.ptcl", n-1)     # n番目のデータファイルの名前の生成
filetitle = sprintf("n = %d", n-1)
plot filename using 2:6 title filetitle with points pt 7 # 1:質量 2:位置 3:速度 4:加速度 5:圧力 6:密度 7:内部エネルギー 8:内部エネルギー時間微分 9:smoothing length
n = n + dn            # ループ変数の増加
if ( n <= n1 ) reread  # ループの評価
undefine n            # ループ変数の削除
