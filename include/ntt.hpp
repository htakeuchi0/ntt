/**
 * @file ntt.hpp
 * @brief Number theoretic transform を実装するヘッダファイル．
 */

#ifndef FFT_NTT_HPP_
#define FFT_NTT_HPP_

/**
 * Number theoretic transform 向け名前空間
 */
namespace ntt {

/** 64ビット整数型 */
using ll = long long int;

/**
 * Number theoretic transform のための基本クラス．
 */
class Ntt {

public:
    /**
     * 次数を返す．
     *
     * @return ll 次数
     */
    virtual ll N() const = 0;

    /**
     * モジュラスを返す．
     *
     * @return ll モジュラス
     */
    virtual ll Mod() const = 0;

    /**
     * 数列の離散フーリエ変換を計算して返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Dft(ll *a) const = 0;

    /**
     * 数列の逆離散フーリエ変換を計算して返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Idft(ll *a) const = 0;

    /**
     * 数列をビット反転で並び替えて返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Reverse(ll *a) const;

    /**
     * 数列の要素ごとの積を計算して返す．
     *
     * @param[in] a 数列．
     * @param[in] b 数列．
     * @param[out] c 数列 a と b の要素ごとの積．
     */
    virtual void MultVec(ll *a, ll *b, ll *c) const;

    /**
     * 数列の畳み込みを計算して返す．
     *
     * @param[in] a 数列．
     * @param[in] b 数列．
     * @param[out] c 数列 a と b の畳み込み．
     */
    virtual void Mult(ll *a, ll *b, ll *c) const;
};

/**
 * Number theoretic transform のための素朴な実装のためのクラス．
 */
class NttNative : public Ntt {

public:
    /**
     * コンストラクタ．
     *
     * @param[in] mod モジュラス．
     * @param[in] omega 1 の n 乗根．
     * @param[in] phi 1 の n 乗根の逆元．
     * @param[in] n 次数．
     * @param[in] n_inv 次数の逆元．
     */
    NttNative(ll mod, ll omega, ll phi, ll n, ll n_inv);

    /**
     * 次数を返す．
     *
     * @return ll 次数
     */
    virtual ll N() const { return n_; }

    /**
     * モジュラスを返す．
     *
     * @return ll モジュラス
     */
    virtual ll Mod() const { return mod_; }

    /**
     * 数列の離散フーリエ変換を計算して返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Dft(ll *a) const;

    /**
     * 数列の逆離散フーリエ変換を計算して返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Idft(ll *a) const;

    /**
     * べき乗を計算して返す．
     *
     * @param[in] x 基数
     * @param[in] k 指数
     * @return ll x の k 乗
     */
    virtual ll Pow(ll x, ll k) const;

private:
    /** モジュラス */
    ll mod_;

    /** 1 の n 乗根 */
    ll omega_;

    /** 1 の n 乗根の逆元 */
    ll phi_;

    /** 次数 */
    ll n_;

    /** 次数の逆元 */
    ll n_inv_;
};

/**
 * Number theoretic transform のための一般的な実装のためのクラス．
 */
class NttBase : public Ntt {

public:
    /**
     * コンストラクタ．
     *
     * @param[in] mod モジュラス．
     * @param[in] omega 1 の n 乗根．
     * @param[in] phi 1 の n 乗根の逆元．
     * @param[in] n 次数．
     * @param[in] n_inv 次数の逆元．
     * @param[in] log_n 次数が 2 の何乗か
     */
    NttBase(ll mod, ll omega, ll phi, ll n, ll n_inv, ll log_n);

    /**
     * 次数を返す．
     *
     * @return ll 次数
     */
    virtual ll N() const { return n_; }

    /**
     * モジュラスを返す．
     *
     * @return ll モジュラス
     */
    virtual ll Mod() const { return mod_; }

    /**
     * 数列の離散フーリエ変換を計算して返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Dft(ll *a) const;

    /**
     * 数列の逆離散フーリエ変換を計算して返す．
     *
     * @param[in, out] a 数列．変換後の数列を上書きして返す．
     */
    virtual void Idft(ll *a) const;

    /**
     * バタフライ演算を実行して結果を返す．
     *
     * @param[in, out] a 要素
     * @param[in, out] b 要素
     * @param[in] k 指数
     */
    virtual void Butterfly(ll& a, ll& b, ll k) const;

    /**
     * 逆離散フーリエ変換でのバタフライ演算を実行して結果を返す．
     *
     * @param[in, out] a 要素
     * @param[in, out] b 要素
     * @param[in] k 指数
     */
    virtual void ButterflyInv(ll& a, ll& b, ll k) const;

    /**
     * べき乗を計算して返す．
     *
     * @param[in] x 基数
     * @param[in] k 指数
     * @return ll x の k 乗
     */
    virtual ll Pow(ll x, ll k) const;

    /**
     * 1 の n 乗根のべき乗を計算して返す．
     *
     * @param[in] k 指数
     * @return ll 1 の n 乗根の k 乗
     */
    virtual ll PowOmega(ll k) const = 0;

    /**
     * 1 の n 乗根の逆元のべき乗を計算して返す．
     *
     * @param[in] k 指数
     * @return ll 1 の n 乗根の逆数の k 乗
     */
    virtual ll PowPhi(ll k) const = 0;

protected:
    /** モジュラス */
    ll mod_;

    /** 1 の n 乗根 */
    ll omega_;

    /** 1 の n 乗根の逆元 */
    ll phi_;

    /** 次数 */
    ll n_;

    /** 次数の逆元 */
    ll n_inv_;

    /** 次数が 2 の何乗か */
    ll log_n_;
};

/**
 * モジュラス 337, 次数 8 のNumber theoretic transform のためのクラス．
 */
class NttMod337Deg8 : public NttBase {

public:
    /** コンストラクタ */
    NttMod337Deg8();

    /**
     * 1 の n 乗根のべき乗を計算して返す．
     *
     * @param[in] k 指数
     * @return ll 1 の n 乗根の k 乗
     */
    virtual ll PowOmega(ll k) const;

    /**
     * 1 の n 乗根の逆元のべき乗を計算して返す．
     *
     * @param[in] k 指数
     * @return ll 1 の n 乗根の逆数の k 乗
     */
    virtual ll PowPhi(ll k) const;

private:
    /** モジュラス */
    static constexpr ll kMod = 337;

    /** 1 の n 乗根 */
    static constexpr ll kOmega = 85;

    /** 1 の n 乗根の逆元 */
    static constexpr ll kPhi = 226;

    /** 次数 */
    static constexpr ll kN = 8;

    /** 次数の逆元 */
    static constexpr ll kNInv = 295;

    /** 次数が 2 の何乗か */
    static constexpr ll kLogN = 3;

    /** 1 の n 乗根のべき乗リスト */
    ll omega_pows_[kN];

    /** 1 の n 乗根の逆数のべき乗リスト */
    ll phi_pows_[kN];
};

/**
 * モジュラス 337, 次数 8 のNumber theoretic transform の素朴な実装のためのクラス．
 */
class NttNativeMod337Deg8 : public NttNative {

public:
    /** コンストラクタ */
    NttNativeMod337Deg8();

private:
    /** モジュラス */
    static constexpr ll kMod = 337;

    /** 1 の n 乗根 */
    static constexpr ll kOmega = 85;

    /** 1 の n 乗根の逆元 */
    static constexpr ll kPhi = 226;

    /** 次数 */
    static constexpr ll kN = 8;

    /** 次数の逆元 */
    static constexpr ll kNInv = 295;
};

/**
 * モジュラス 19529729, 次数 131072 のNumber theoretic transform のためのクラス．
 */
class NttMod19529729Deg131072 : public NttBase {

public:
    /* コンストラクタ */
    NttMod19529729Deg131072();

    /**
     * 1 の n 乗根のべき乗を計算して返す．
     *
     * @param[in] k 指数
     * @return ll 1 の n 乗根の k 乗
     */
    virtual ll PowOmega(ll k) const;

    /**
     * 1 の n 乗根の逆元のべき乗を計算して返す．
     *
     * @param[in] k 指数
     * @return ll 1 の n 乗根の逆数の k 乗
     */
    virtual ll PowPhi(ll k) const;

private:
    /** モジュラス */
    static constexpr ll kMod = 19529729;

    /** 1 の n 乗根 */
    static constexpr ll kOmega = 770;

    /** 1 の n 乗根の逆元 */
    static constexpr ll kPhi = 16765131;

    /** 次数 */
    static constexpr ll kN = 131072;

    /** 次数の逆元 */
    static constexpr ll kNInv = 19529580;

    /** 次数が 2 の何乗か */
    static constexpr ll kLogN = 17;
};

/**
 * モジュラス 19529729, 次数 131072 のNumber theoretic transform の素朴な実装のためのクラス．
 */
class NttNativeMod19529729Deg131072 : public NttNative {

public:
    /* コンストラクタ */
    NttNativeMod19529729Deg131072();

private:
    /** モジュラス */
    static constexpr ll kMod = 19529729;

    /** 1 の n 乗根 */
    static constexpr ll kOmega = 770;

    /** 1 の n 乗根の逆元 */
    static constexpr ll kPhi = 16765131;

    /** 次数 */
    static constexpr ll kN = 131072;

    /** 次数の逆元 */
    static constexpr ll kNInv = 19529580;

    /** 次数が 2 の何乗か */
    static constexpr ll kLogN = 17;
};

} // namespace ntt

#endif // #ifndef FFT_NTT_HPP_
