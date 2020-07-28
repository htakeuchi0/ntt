/**
 * @file util.cpp
 * @brief ユーティリティクラスを定義するソースファイル．
 */

#include "include/util.hpp"

/*
 * Number theoretic transform 向け名前空間
 */
namespace ntt {

/*
 * 逆数を返す．
 *
 * @param [in] x 値
 * @param [in] n モジュラス
 * @return ll 逆数
 */
ll Utility::InvMod(ll x, ll n) {
    ll a = n;
    ll b = x;
    ll b_pre = 0;

    ll q = a / b;
    ll r = a - b * q;
    ll b_current = 1;

    ll tmp = b_current;

    ll qn = (q / n) + 1;
    b_current = (b_pre + (((qn * n - q)*b_current) % n)) % n;
    b_pre = tmp;
    a = b;
    b = r;

    while (r > 1) {
        q = a / b;
        r = a - b * q;

        ll tmp = b_current;
        ll qn = (q / n) + 1;
        b_current = (b_pre + (((qn * n - q)*b_current) % n)) % n;
        b_pre = tmp;
        a = b;
        b = r;
    }

    return b_current % n;
}

} // namespace ntt
