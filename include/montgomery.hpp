/**
 * @file montgomery.hpp
 * @brief モンゴメリ乗算を実装するためのヘッダファイル．
 */

#ifndef FFT_MONTGOMERY_HPP_
#define FFT_MONTGOMERY_HPP_

/*
 * Number theoretic transform 向け名前空間
 */
namespace ntt {

/* 64ビット整数型 */
using ll = long long int;

/**
 * モンゴメリ乗算を行うためのクラス
 */
class Montgomery {

public:
    /**
     * コンストラクタ．
     *
     * @param [in] n モジュラスN
     * @param [in] log2r R>N が 2 の何乗か
     */
    Montgomery(int n, int log2r);

    /**
     * mod R で NN' = -1 を満たす N を計算して返す．
     */
    ll ComputeNn(); 

    /**
     * モンゴメリリダクションを返す．
     *
     * @param [in] t リダクションを計算する値
     */
    ll Reduction(ll t);

    /**
     * mod N で積を計算して返す．
     *
     * @param [in] a 値
     * @param [in] b 値
     * @return ll a と b の積
     */
    ll Mult(ll a, ll b);

    /**
     * mod N でべき乗を計算して返す．
     *
     * @param [in] a 基数
     * @param [in] k 指数
     * @return ll a の b 乗
     */
    ll Pow(ll a, ll k);

    /**
     * モジュラスを返す．
     *
     * @return ll モジュラス
     */
    ll N() { return n_; }

private:
    /** モジュラス N */
    ll n_;

    /** R > N を満たす2べきの値 */
    ll r_;

    /** R が 2 の何乗か */
    ll log2r_;

    /** mod N における R の 2 乗 */
    ll r2_;

    /** mod R における NN' = -1 を満たす N' */
    ll nn_;
};

/**
 * N=19529729, R=2^25 としたときのモンゴメリ乗算を行うためのクラス
 */
class MontgomeryMod19529729R25 : public Montgomery {

public:
    /** コンストラクタ． */
    MontgomeryMod19529729R25() : Montgomery(19529729, 25) {}
};

} // namespace ntt

#endif // #ifndef FFT_MONTGOMERY_HPP_
