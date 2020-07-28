/**
 * @file util.hpp
 * @brief ユーティリティクラスを定義するヘッダファイル．
 */

#ifndef FFT_UTIL_HPP_
#define FFT_UTIL_HPP_

/*
 * Number theoretic transform 向け名前空間
 */
namespace ntt {

/** 64ビット整数型 */
using ll = long long int;

/**
 * ユーティリティクラス．
 */
class Utility {

public:
    /**
     * 逆数を返す．
     *
     * @param[in] x 値
     * @param[in] n モジュラス
     * @return ll 逆数
     */
    static ll InvMod(ll x, ll n);
};

} // namespace ntt

#endif // #ifndef FFT_UTIL_HPP_
