# Number theoretic transform サンプル実装

## 概要

Number theoretic transform (NTT) のサンプル実装です．

動作環境は以下のとおりです．
* OS: Linux (Ubuntu 18.04.3 LTS)
* コンパイラ: g++ (Ubuntu 7.4.0-1ubuntu1~18.04.1) 7.4.0

## ディレクトリ構成

```
|- ntt/
   |- Doxyfile               - doxygen設定ファイル
   |- Makefile               - makeファイル
   |- README.md              - 本ファイル
   |- makeenv.sh             - 環境構築用スクリプト
   |- include/               - ヘッダファイル
   |  |- montgomery.hpp
   |  |- ntt.hpp
   |  |- util.hpp
   |
   |- main/                  - メインファイル
   |  |- main.cpp
   |
   |- src/                   - ソースファイル
   |  |- montgomery.cpp
   |  |- ntt.cpp
   |  |- util.cpp
   |
   |- test/                  - テストファイル
      |- gtest_montgomery.cpp
```

## 準備と使いかた

### 準備

環境構築のために以下を実行して下さい．

```
$ source makeenv.sh
```

スクリプト実行には `cmake` が必要です．    
エラーが発生したら，以下を実行して下さい．
```
$ sudo apt install cmake
```

## コンパイルと実行

メインプログラムのコンパイルとプログラムは，以下のコマンドで実行できます．

```
$ make
$ ./main.o
```

## テストの実行

テストプログラムのコンパイルと単体テストの実行は以下のコマンドで実行できます．    
`make` は `make test` でも可能です．

```
$ make 
$ ./test.o
```

## Doxygenの生成

以下のコマンドで詳細仕様が記述された html ファイルが生成できます．    
`doxygen/html/index.html` をブラウザで開いて読むことができます．

```
$ make docs
```

## 実行結果

`Z/19529729Z` における 131072 (=2^17) 次多項式の畳込みの計算を，    
以下の3つの方法で実行した場合の実行時間を記載します．
* 素朴な畳み込み (Naive)
* NTT (NTT)
* モンゴメリ乗算を利用したNTT (NTT+Montgomery)

実行環境は以下のとおりです．    
* CPU: Intel(R) Core(TM) it-6200U CPU @ 2.30GHz
* メモリ: 4GB
* OS: Linux (Ubuntu 18.04.3 LTS)

テスト問題は以下のとおりです．
* `a = [1, 2, 3, 4, 0, ..., 0]`
* `b = [1, 2, 3, 4, 0, ..., 0]`

したがって，積 `ab` は以下のとおりです．
* `ab = [1, 4, 10, 20, 25, 24, 16, 0, ..., 0]`

実行時間は以下のとおりです．    
高速フーリエ変換と，モンゴメリ乗算による高速化の効果が現れています．

|方法          |時間 [ms]|備考         |
|--------------|---------|-------------|
|Naive         |>28800000|8時間を超えた|
|NTT           |943.1    |10回の平均   |
|NTT+Montgomery|309.2    |10回の平均   |

作者環境でのプログラムの実行結果は以下のとおりです．

```
$ ./main.o
---- NTT (Basic)       ----
a    : 1	2	3	4	0	0	0	
b    : 1	2	3	4	0	0	0	
a * b: 1	4	10	20	25	24	16	

running......end

elapsed: 
- repeat: 10
- mean  : 943.1 [ms]
- std   : 1.86815

---- NTT (Montgomery)  ----
a    : 1	2	3	4	0	0	0	
b    : 1	2	3	4	0	0	0	
a * b: 1	4	10	20	25	24	16	

running......end

elapsed: 
- repeat: 10
- mean  : 309.2 [ms]
- std   : 4.0694
```
