#ifndef model_GAUSSLEGENDRE_2_16_H_
#define model_GAUSSLEGENDRE_2_16_H_

namespace model{

template <>
struct GaussLegendre<2,16>{
    // Dunavant p=8, ng=16 (weights sum to 0.5 on reference triangle)
    // Source: Dunavant 1985 tables / Burkardt triangle_dunavant_rule. :contentReference[oaicite:3]{index=3}
    static Eigen::Matrix<double,3,16> abcsissasAndWeights(){
        Eigen::Matrix<double,16,3> U;

        const double w0 = 0.5*0.144315607677787;

        const double w1 = 0.5*0.095091634267285;
        const double a1 = 0.081414823414554;
        const double b1 = 0.459292588292723;

        const double w2 = 0.5*0.103217370534718;
        const double a2 = 0.658861384496480;
        const double b2 = 0.170569307751760;

        const double w3 = 0.5*0.032458497623198;
        const double a3 = 0.898905543365938;
        const double b3 = 0.050547228317031;

        const double w4 = 0.5*0.027230314174435;
        const double a4 = 0.008394777409958;
        const double b4 = 0.263112829634638;
        const double c4 = 0.728492392955404;

        U <<
        // centroid
        0.333333333333333, 0.333333333333333, w0,

        // (a1,b1,b1)
        b1, b1, w1,
        a1, b1, w1,
        b1, a1, w1,

        // (a2,b2,b2)
        b2, b2, w2,
        a2, b2, w2,
        b2, a2, w2,

        // (a3,b3,b3)
        b3, b3, w3,
        a3, b3, w3,
        b3, a3, w3,

        // (a4,b4,c4) six permutations
        b4, c4, w4,
        c4, b4, w4,
        a4, c4, w4,
        c4, a4, w4,
        a4, b4, w4,
        b4, a4, w4;

        return U.transpose();
    }
};

} // namespace model

#endif