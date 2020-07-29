/**
 * @file main.cpp
 * @brief メインメソッドをもつソースファイル．
 */

#include "include/util.hpp"
#include "include/montgomery.hpp"
#include "include/ntt.hpp"
#include <array>
#include <iostream>
#include <chrono>
#include <cmath>

namespace {

/**
 * Number theoretic transform の実行サンプルを出力する．
 *
 * @param[in] ntt NTTオブジェクト
 * @param[in] is_show_mode 標準出力する場合true
 * @return double 実行時間 [ms]
 */
double NttSample(const ntt::Ntt& ntt, bool is_show_mode) {
    int size = ntt.N();

    ntt::ll *a = new ntt::ll[size];
    ntt::ll *b = new ntt::ll[size];
    ntt::ll *c = new ntt::ll[size];

    for (int i = 0; i < size; i++) {
        a[i] = 0;
        b[i] = 0;
        c[i] = 0;
    }

    a[0] = b[0] = 1;
    a[1] = b[1] = 2;
    a[2] = b[2] = 3;
    a[3] = b[3] = 4;
    
    if (is_show_mode) {
        std::cout << "a    : ";
        for (int i = 0; i < 7; i++) {
            std::cout << a[i] << "\t";
        }
        std::cout << std::endl;

        std::cout << "b    : ";
        for (int i = 0; i < 7; i++) {
            std::cout << b[i] << "\t";
        }
        std::cout << std::endl;
    }

    auto begin = std::chrono::system_clock::now();
    ntt.Mult(a, b, c);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    double elapsed_time = static_cast<double>(elapsed.count());

    if (is_show_mode) {
        std::cout << "a * b: ";
        for (int i = 0; i < 7; i++) {
            std::cout << c[i] << "\t";
        }
        std::cout << std::endl;
    }

    delete[] a;
    a = nullptr;

    delete[] b;
    b = nullptr;

    delete[] c;
    c = nullptr;

    return elapsed_time;
}

/**
 * 素朴な実装の畳み込み演算の実行サンプルを出力する．
 *
 * @param[in] is_show_mode 標準出力する場合true
 * @return double 実行時間 [ms]
 */
double NttNaiveSample(bool is_show_mode) {
    ntt::NttNaiveMod19529729Deg131072 ntt;
    double elapsed = NttSample(ntt, is_show_mode);
    return elapsed;
}

/**
 * 素朴な実装の Number theoretic transform の実行サンプルを出力する．
 *
 * @param[in] is_show_mode 標準出力する場合true
 * @return double 実行時間 [ms]
 */
double NttBasicSample(bool is_show_mode) {
    ntt::NttMod19529729Deg131072 ntt;
    double elapsed = NttSample(ntt, is_show_mode);
    return elapsed;
}

/**
 * モンゴメリ乗算を利用した Number theoretic transform の実行サンプルを出力する．
 *
 * @param[in] is_show_mode 標準出力する場合true
 * @return double 実行時間 [ms]
 */
double NttMontgomerySample(bool is_show_mode) {
    ntt::NttMod19529729Deg131072M ntt;
    double elapsed = NttSample(ntt, is_show_mode);
    return elapsed;
}

/**
 * 配列の要素の平均を計算して返す．
 *
 * @param[in] array 配列
 * @return double 平均
 */
template <size_t N>
double Mean(const std::array<double, N>& array) {
    double mean = 0.0;
    for (double entry : array) {
        mean += entry;
    }
    mean /= array.size();
    return mean;
}

/**
 * 配列の要素の標準偏差を計算して返す．
 *
 * @param[in] array 配列
 * @return double 標準偏差
 */
template <size_t N>
double Std(const std::array<double, N>& array) {
    double mean = Mean(array);

    double std = 0.0;
    for (double entry : array) {
        double diff = entry - mean;
        std += diff * diff;
    }
    std /= array.size();
    std = ::sqrt(std);
    return std;
}

/**
 * サンプルを表示する．
 *
 * @param[in] sample_name サンプル名
 * @param[in] sample サンプルメソッド
 * @param[in] times 実行回数
 */
void ShowSample(std::string sample_name, double (*sample)(bool), int times = 10) {
    std::cout << sample_name << std::endl;
    std::array<double, 10> elapsed_times;

    double elapsed_time = sample(true);
    elapsed_times[0] = elapsed_time;

    std::cout << "\n";
    std::cout << "running......" << std::flush;
    for (int i = 1; i < times; i++) {
        double elapsed_time = sample(false);
        elapsed_times[i] = elapsed_time;
    }
    std::cout << "end\n" << std::endl;

    double mean = Mean(elapsed_times);
    double std = Std(elapsed_times);

    std::cout << "elapsed: \n";
    std::cout << "- repeat: " << times << "\n";
    std::cout << "- mean  : " << mean << " [ms]\n";
    std::cout << "- std   : " << std << "\n";
    std::cout << std::endl;
}

} // namespace

/**
 * ヘルプを表示する．
 *
 * @param[in, out] out 出力先
 * @param[in] program_name プログラム名
 */
void ShowHelp(std::ostream& out, const std::string& program_name) {
    out << "Usage:\n";
    out << "$ " << program_name << " [--naive]\n";
    out << "$ " << program_name << " [--help]\n";
    out << "\n";
    out << "Arguments:\n";
    out << "--naive, -n : Also execute the naive NTT\n";
    out << "--help, -h  : Show the help message and exit" << std::endl;
}

/**
 * メインメソッド
 *
 * @param[in] argc コマンドライン引数の数
 * @param[in] argv コマンドライン引数
 * @return int 終了コード
 */
int main(int argc, char **argv) {
    bool do_naive = false;
    if (argc >= 2) {
        std::string arg(argv[1]);
        do_naive = (arg == std::string("--naive") || arg == std::string("-n"));

        if (arg == std::string("--help") || arg == std::string("-h")) {
            ShowHelp(std::cout, std::string(argv[0]));
            return 0;
        }
    }

    if (do_naive) {
        int times = 1;
        ShowSample("---- NTT (Naive)       ----", NttNaiveSample, times);
    }

    ShowSample("---- NTT (Basic)       ----", NttBasicSample);
    ShowSample("---- NTT (Montgomery)  ----", NttMontgomerySample);
    return 0;
}
