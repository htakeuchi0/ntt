/**
 * @file montgomery.cpp
 * @brief モンゴメリ乗算を実装するためのソースファイル．
 */

#include "include/montgomery.hpp"

/*
 * Number theoretic transform 向け名前空間
 */
namespace ntt {

/*
 * コンストラクタ．
 *
 * @param[in] n モジュラスN
 * @param[in] log2r R>N が 2 の何乗か
 */
Montgomery::Montgomery(int n, int log2r) : n_(n), r_(1 << log2r), log2r_(log2r) {
    nn_ = ComputeNn(); 
    r2_ = (r_ * r_) % n_;
}

/*
 * mod R で NN' = -1 を満たす N を計算して返す．
 *
 * @return ll mod R で NN' = -1 を満たす N
 */
ll Montgomery::ComputeNn() const {
    ll nn = 0;
    ll t = 0;
    for (int i = 0; i < log2r_; i++) {
        if ((t & 1) == 0) {
            t += n_;
            nn += (1 << i);
        }
        t >>= 1;
    }

    return nn;
}

/*
 * モンゴメリリダクションを返す．
 *
 * @param[in] t リダクションを計算する値
 * @return ll モンゴメリリダクション
 */
ll Montgomery::Reduction(ll t) const {
    t = (t + ((t * nn_) & (r_ - 1))*n_) >> log2r_;
    return (t >= n_) ? t - n_ : t;
}

/*
 * mod N で積を計算して返す．
 *
 * @param[in] a 値
 * @param[in] b 値
 * @return ll a と b の積
 */
ll Montgomery::Mult(ll a, ll b) const {
    ll c = Reduction(Reduction(a * b) * r2_);
    return c;
}

/*
 * mod N でべき乗を計算して返す．
 *
 * @param[in] a 基数
 * @param[in] k 指数
 * @return ll a の b 乗
 */
ll Montgomery::Pow(ll a, ll k) const {
    ll p = a;
    ll v = 1;
    if (k == 0) {
        return v;
    }

    ll mr_p = Reduction(p * r2_);
    ll mr_v = Reduction(v * r2_);

    while (k >= 1) {
        if ((k & 1) == 1) {
            mr_v = Reduction(mr_v * mr_p);
        }
        k >>= 1;
        mr_p = Reduction(mr_p * mr_p);
    }

    v = Reduction(mr_v);

    return v;
}

} // namespace ntt
