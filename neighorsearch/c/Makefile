# コンパイラの指定
CC = gcc
# 最適化オプション
CFLAGS += -O3 -funroll-loops -ffast-math -march=native

# ライブラリリンク
LIBS = -lm

# main() の中でマクロ定義した #define はここで指定するようにしても良い
# ただしどちらか片方のみ
# The number of particles
# MACROS  = -DN_PTCL=1024
# Dimension
# MACROS += -DDIM=1
# The total mass of sphere
# MACROS += -DM_SPHERE=1.0
# The radius of sphere
# MACROS += -DR_SPHERE=1.0
# virial ratio
# MACROS += -Dr_VIRIAL=0.5

# ソースファイル（main.c）と実行ファイル名の指定
SOURCE  = main.c
PROGRAM = search.out

$(PROGRAM):
	$(CC) $(SOURCE) $(CFLAGS) $(MACROS) $(LIBS) -o $(PROGRAM)

.PHONY: clean
clean:
	rm -f *.out *.log
