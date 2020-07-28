/**
 * @file gtest_montgomery.cpp
 * @brief モンゴメリ乗算のテストファイル．
 */

#include "gtest/gtest.h"
#include "include/montgomery.hpp"

namespace ntt {

/**
 * テストケースのサンプルクラス．
 */
class MontgomeryTest : public ::testing::Test {
protected:
    /**
     * べき乗を計算して返す．
     *
     * @param [in] x 基数
     * @param [in] k 指数
     * @param [in] n モジュラス
     * @return ll x の k 乗
     */
    ll Pow(ll x, ll k, ll n);
};

/*
 * モンゴメリ乗算による積が正しく計算できることを確認する．
 */
TEST_F(MontgomeryTest, Mult) {
    MontgomeryMod19529729R25 montgomery;

    ll a = 12345678;
    ll b = 87654321;
    ll actual = montgomery.Mult(a, b);

    ll expected = (a * b) % montgomery.N();

    ASSERT_EQ(expected, actual);
}

/*
 * モンゴメリ乗算によるべき乗が正しく計算できることを確認する．
 */
TEST_F(MontgomeryTest, Pow) {
    MontgomeryMod19529729R25 montgomery;

    ll a = 12345678;
    ll b = 87654321;
    ll actual = montgomery.Pow(a, b);

    ll expected = Pow(a, b, montgomery.N());

    ASSERT_EQ(expected, actual);
}

/*
 * べき乗を計算して返す．
 *
 * @param [in] x 基数
 * @param [in] k 指数
 * @param [in] n モジュラス
 * @return ll x の k 乗
 */
ll MontgomeryTest::Pow(ll x, ll k, ll n) {
    ll p = x;
    ll v = 1;
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

} // namespace ntt

