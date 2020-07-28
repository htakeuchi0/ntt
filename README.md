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

