/**
 * @file main.cpp
 * @brief メインメソッドをもつソースファイル．
 */

#include "include/util.hpp"
#include "include/montgomery.hpp"
#include "include/ntt.hpp"
#include <iostream>

namespace {

/**
 * べき乗を計算して返す．
 *
 * @param [in] x 基数
 * @param [in] k 指数
 * @param [in] n モジュラス
 * @return ll x の k 乗
 */
ntt::ll Pow(ntt::ll x, ntt::ll k, ntt::ll n) {
    ntt::ll p = x;
    ntt::ll v = 1;
    if (k == 0) {
        return v;
    }

    while (k >= 1) {
        if ((k & 1) == 1) {
            v = (v * p) % n;
        }
        k >>= 1;
        p = (p * p) % n;
    }

    return v;
}

/**
 * モンゴメリ乗算の実行サンプルを出力する．
 */
void MontgomerySample() {
    ntt::NttNativeMod19529729Deg131072 ntt;
    ntt::MontgomeryMod19529729R25 montgomery;

    ntt::ll a = 12345678;
    ntt::ll b = 87654321;
    ntt::ll ab = montgomery.Mult(a, b);
    ntt::ll a_pow_b = montgomery.Pow(a, b);

    std::cout << "(montgomery) ab = " << ab << std::endl;
    std::cout << "(naive     ) ab = " << (a * b) % montgomery.N() << std::endl;

    std::cout << "(montgomery) a^b = " << a_pow_b << std::endl;
    std::cout << "(naive     ) a^b = " << Pow(a, b, montgomery.N()) << std::endl;
}

/**
 * Number theoretic transform の実行サンプルを出力する．
 */
void NttSample() {
    ntt::NttMod19529729Deg131072 ntt;
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

    ntt.Mult(a, b, c);
    std::cout << "a * b: ";
    for (int i = 0; i < 7; i++) {
        std::cout << c[i] << "\t";
    }
    std::cout << std::endl;

    delete[] a;
    a = nullptr;

    delete[] b;
    b = nullptr;

    delete[] c;
    c = nullptr;
}

/**
 * サンプルを表示する．
 *
 * @param [in] sample_name サンプル名
 * @param [in] sample サンプルメソッド
 */
void ShowSample(std::string sample_name, void (*sample)(void)) {
    std::cout << sample_name << std::endl;
    sample();
    std::cout << std::endl;
}

} // namespace

/**
 * メインメソッド
 *
 * @param[in] argc コマンドライン引数の数
 * @param[in] argv コマンドライン引数
 * @return int 終了コード
 */
int main(int argc, char **argv) {
    ShowSample("---- Montgomery ----", MontgomerySample);
    ShowSample("----    Ntt     ----", NttSample);
    return 0;
}
