if(exist("n")==0 || n<0) n = n0             # ループ変数の初期化
filetitle = sprintf("n = %d", n-1)
filename = sprintf("data_%d.ptcl", n-1)     # n番目のデータファイルの名前の読み取り
splot filename title filetitle
n = n + dn            # ループ変数の増加
if ( n < n1 ) reread  # ループの評価
undefine n            # ループ変数の削除
