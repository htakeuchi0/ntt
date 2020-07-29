/**
 * @file ntt.cpp
 * @brief Number theoretic transform を実装するソースファイル．
 */

#include "include/ntt.hpp"
#include <iostream>

/*
 * Number theoretic transform 向け名前空間
 */
namespace ntt {

/*
 * 数列をビット反転で並び替えて返す．
 *
 * @param [in/out] a 数列．変換後の数列を上書きして返す．
 */
void Ntt::Reverse(ll *a) const {
    ll j = 0;
    ll n = N();
    for (ll i = 0; i < n; i++) {
        if (j > i) {
            ll tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }

        ll m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

/*
 * 数列の要素ごとの積を計算して返す．
 *
 * @param [in] a 数列．
 * @param [in] b 数列．
 * @param [out] c 数列 a と b の要素ごとの積．
 */
void Ntt::MultVec(ll *a, ll *b, ll *c) const {
    ll n = N();
    ll mod = Mod();

    for (ll i = 0; i < n; i++) {
        c[i] = (a[i] * b[i]) % mod;
    }
}

/*
 * 数列の畳み込みを計算して返す．
 *
 * @param [in] a 数列．
 * @param [in] b 数列．
 * @param [out] c 数列 a と b の畳み込み．
 */
void Ntt::Mult(ll *a, ll *b, ll *c) const {
    Dft(a);
    Dft(b);
    MultVec(a, b, c);
    Idft(c);
}

/*
 * 数列の離散フーリエ変換を計算して返す．
 *
 * @param [in/out] 数列．変換後の数列を上書きして返す．
 */
void NttNaive::Dft(ll *a) const {
    ll *c = new ll[n_];
    for (ll i = 0; i < n_; i++) {
        c[i] = 0;
        for (ll j = 0; j < n_; j++) {
            c[i] += Pow(omega_, j*i) * a[j];
        }
        c[i] %= mod_;
    }

    for (ll i = 0; i < n_; i++) {
        a[i] = c[i];
    }

    delete[] c; c = nullptr;
}

/*
 * 数列の逆離散フーリエ変換を計算して返す．
 *
 * @param [in/out] a 数列．変換後の数列を上書きして返す．
 */
void NttNaive::Idft(ll *a) const {
    ll *c = new ll[n_];
    for (ll i = 0; i < n_; i++) {
        c[i] = 0;
        for (ll j = 0; j < n_; j++) {
            c[i] += Pow(phi_, j*i) * a[j];
        }
        c[i] %= mod_;
    }

    for (ll i = 0; i < n_; i++) {
        c[i] = (c[i] * n_inv_) % mod_;
    }

    for (ll i = 0; i < n_; i++) {
        a[i] = c[i];
    }

    delete[] c; c = nullptr;
}

/*
 * べき乗を計算して返す．
 *
 * @param [in] x 基数
 * @param [in] k 指数
 * @return ll x の k 乗
 */
ll NttNaive::Pow(ll x, ll k) const {
    ll p = x;
    ll v = 1;
    if (k == 0) {
        return v;
    }

    while (k >= 1) {
        if ((k & 1) == 1) {
            v = (v * p) % mod_;
        }
        k >>= 1;
        p = (p * p) % mod_;
    }

    return v;
}

/*
 * 数列の離散フーリエ変換を計算して返す．
 *
 * @param [in/out] 数列．変換後の数列を上書きして返す．
 */
void NttBase::Dft(ll *a) const {
    Reverse(a);

    ll m = log_n_;

    for (ll l = 1; l <= m; l++) {
        ll max_q = (1 << (m - l));
        for (ll q = 0; q < max_q; q++) {
            ll max_r = (1 << (l - 1));
            for (ll r = 0; r < max_r; r++) {
                ll k = (q << l) + r;
                Butterfly(a[k], a[k + max_r], r << (m - l));
            }
        }
    }
}

/*
 * 数列の逆離散フーリエ変換を計算して返す．
 *
 * @param [in/out] a 数列．変換後の数列を上書きして返す．
 */
void NttBase::Idft(ll *a) const {
    Reverse(a);

    ll m = log_n_;

    for (ll l = 1; l <= m; l++) {
        ll max_q = (1 << (m - l));
        for (ll q = 0; q < max_q; q++) {
            ll max_r = (1 << (l - 1));
            for (ll r = 0; r < max_r; r++) {
                ll k = (q << l) + r;
                ButterflyInv(a[k], a[k + max_r], r << (m - l));
            }
        }
    }

    for (ll i = 0; i < n_; i++) {
        a[i] = (a[i] * n_inv_) % mod_;
    }
}

/*
 * バタフライ演算を実行して結果を返す．
 *
 * @param [in/out] a 要素
 * @param [in/out] b 要素
 * @param [in] k 指数
 */
void NttBase::Butterfly(ll& a, ll& b, ll k) const {
    ll tmp = (PowOmega(k) * b) % mod_;
    ll minus_tmp = (mod_ - tmp) % mod_;

    b = (a + minus_tmp) % mod_;
    a = (a + tmp) % mod_;
}

/*
 * 逆離散フーリエ変換でのバタフライ演算を実行して結果を返す．
 *
 * @param [in/out] a 要素
 * @param [in/out] b 要素
 * @param [in] k 指数
 */
void NttBase::ButterflyInv(ll& a, ll& b, ll k) const {
    ll tmp = (PowPhi(k) * b) % mod_;
    ll minus_tmp = (mod_ - tmp) % mod_;

    b = (a + minus_tmp) % mod_;
    a = (a + tmp) % mod_;
}

/*
 * べき乗を計算して返す．
 *
 * @param [in] x 基数
 * @param [in] k 指数
 * @return ll x の k 乗
 */
ll NttBase::Pow(ll x, ll k) const {
    ll p = x;
    ll v = 1;
    if (k == 0) {
        return v;
    }

    while (k >= 1) {
        if ((k & 1) == 1) {
            v = (v * p) % mod_;
        }
        k >>= 1;
        p = (p * p) % mod_;
    }

    return v;
}

/*
 * 1 の n 乗根のべき乗を計算して返す．
 *
 * @param [in] k 指数
 * @return ll 1 の n 乗根の k 乗
 */
ll NttMod337Deg8::PowOmega(ll k) const {
    return omega_pows_[k % n_];
}

/*
 * 1 の n 乗根の逆元のべき乗を計算して返す．
 *
 * @param [in] k 指数
 * @return ll 1 の n 乗根の逆数の k 乗
 */
ll NttMod337Deg8::PowPhi(ll k) const {
    return phi_pows_[k % n_];
}

/*
 * コンストラクタ．
 *
 * @param [in] mod モジュラス．
 * @param [in] omega 1 の n 乗根．
 * @param [in] phi 1 の n 乗根の逆元．
 * @param [in] n 次数．
 * @param [in] n_inv 次数の逆元．
 */
NttNaive::NttNaive(ll mod, ll omega, ll phi, ll n, ll n_inv) : 
        mod_(mod), omega_(omega), phi_(phi), n_(n), n_inv_(n_inv){}

/*
 * コンストラクタ．
 *
 * @param [in] mod モジュラス．
 * @param [in] omega 1 の n 乗根．
 * @param [in] phi 1 の n 乗根の逆元．
 * @param [in] n 次数．
 * @param [in] n_inv 次数の逆元．
 * @param [in] log_n 次数が 2 の何乗か
 */
NttBase::NttBase(ll mod, ll omega, ll phi, ll n, ll n_inv, ll log_n) :
        mod_(mod),
        omega_(omega),
        phi_(phi),
        n_(n),
        n_inv_(n_inv),
        log_n_(log_n) {}

/* コンストラクタ */
NttMod337Deg8::NttMod337Deg8() :
        NttBase(kMod, kOmega, kPhi, kN, kNInv, kLogN) {

    ll omega_pows[] { 85, 148, 111, 136, 252, 189, 226 };

    omega_pows_[0] = 1;
    phi_pows_[0] = 1;

    for (ll i = 1; i < n_; i++) {
        omega_pows_[i] = omega_pows[i - 1];
        phi_pows_[i] = omega_pows[n_ - i - 1];
    }
}

/* コンストラクタ */
NttNaiveMod337Deg8::NttNaiveMod337Deg8() :
        NttNaive(kMod, kOmega, kPhi, kN, kNInv) {}

/* コンストラクタ */
NttMod19529729Deg131072::NttMod19529729Deg131072() :
        NttBase(kMod, kOmega, kPhi, kN, kNInv, kLogN) {
}

/* コンストラクタ */
NttMod19529729Deg131072M::NttMod19529729Deg131072M() :
        NttBase(kMod, kOmega, kPhi, kN, kNInv, kLogN) {
}

/* コンストラクタ */
NttNaiveMod19529729Deg131072::NttNaiveMod19529729Deg131072() :
        NttNaive(kMod, kOmega, kPhi, kN, kNInv) {
}

/*
 * 1 の n 乗根のべき乗を計算して返す．
 *
 * @param [in] k 指数
 * @return ll 1 の n 乗根の k 乗
 */
ll NttMod19529729Deg131072::PowOmega(ll k) const {
    return Pow(omega_, k);
}

/*
 * 1 の n 乗根の逆元のべき乗を計算して返す．
 *
 * @param [in] k 指数
 * @return ll 1 の n 乗根の逆数の k 乗
 */
ll NttMod19529729Deg131072::PowPhi(ll k) const {
    return Pow(phi_, k);
}

/*
 * 数列の要素ごとの積を計算して返す．
 *
 * @param [in] a 数列．
 * @param [in] b 数列．
 * @param [out] c 数列 a と b の要素ごとの積．
 */
void NttMod19529729Deg131072M::MultVec(ll *a, ll *b, ll *c) const {
    for (ll i = 0; i < n_; i++) {
        c[i] = montgomery_.Mult(a[i], b[i]);
    }
}

/*
 * 数列の逆離散フーリエ変換を計算して返す．
 *
 * @param [in/out] a 数列．変換後の数列を上書きして返す．
 */
void NttMod19529729Deg131072M::Idft(ll *a) const {
    Reverse(a);

    ll m = log_n_;

    for (ll l = 1; l <= m; l++) {
        ll max_q = (1 << (m - l));
        for (ll q = 0; q < max_q; q++) {
            ll max_r = (1 << (l - 1));
            for (ll r = 0; r < max_r; r++) {
                ll k = (q << l) + r;
                ButterflyInv(a[k], a[k + max_r], r << (m - l));
            }
        }
    }

    for (ll i = 0; i < n_; i++) {
        a[i] = montgomery_.Mult(a[i], n_inv_);
    }
}

/*
 * バタフライ演算を実行して結果を返す．
 *
 * @param [in/out] a 要素
 * @param [in/out] b 要素
 * @param [in] k 指数
 */
void NttMod19529729Deg131072M::Butterfly(ll& a, ll& b, ll k) const {
    ll tmp = montgomery_.Mult(PowOmega(k), b);
    ll minus_tmp = mod_ - tmp;

    b = a + minus_tmp;
    b = (b >= mod_) ? b - mod_ : b;

    a = a + tmp;
    a = (a >= mod_) ? a - mod_ : a;
}

/*
 * 逆離散フーリエ変換でのバタフライ演算を実行して結果を返す．
 *
 * @param [in/out] a 要素
 * @param [in/out] b 要素
 * @param [in] k 指数
 */
void NttMod19529729Deg131072M::ButterflyInv(ll& a, ll& b, ll k) const {
    ll tmp = montgomery_.Mult(PowPhi(k), b);
    ll minus_tmp = mod_ - tmp;

    b = a + minus_tmp;
    b = (b >= mod_) ? b - mod_ : b;

    a = a + tmp;
    a = (a >= mod_) ? a - mod_ : a;
}

/*
 * 1 の n 乗根のべき乗を計算して返す．
 *
 * @param [in] k 指数
 * @return ll 1 の n 乗根の k 乗
 */
ll NttMod19529729Deg131072M::PowOmega(ll k) const {
    return montgomery_.Pow(omega_, k);
}

/*
 * 1 の n 乗根の逆元のべき乗を計算して返す．
 *
 * @param [in] k 指数
 * @return ll 1 の n 乗根の逆数の k 乗
 */
ll NttMod19529729Deg131072M::PowPhi(ll k) const {
    return montgomery_.Pow(phi_, k);
}

} // namespace ntt
