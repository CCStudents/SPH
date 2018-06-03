# コンパイラの指定
CC = g++
# 最適化オプション
CFLAGS += -O3 -funroll-loops -ffast-math -march=native

# ライブラリリンク
LIBS = -lm

# ソースファイル（main.c）と実行ファイル名の指定
CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.hpp)
PROGRAM = search.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(CPPOBJS) $(LIBS) -o $(PROGRAM)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CC) -c $< $(CFLAGS) $(LIBS)
	@echo "[$< OK]"

clean:
	-rm *.out *.o
