#ifndef model_GAUSSLEGENDRE_2_13_H_
#define model_GAUSSLEGENDRE_2_13_H_

namespace model{

template <>
struct GaussLegendre<2,13>{
    // Dunavant p=7, ng=13 (weights sum to 0.5 on reference triangle)
    // Source: Dunavant 1985 tables / Burkardt triangle_dunavant_rule. :contentReference[oaicite:2]{index=2}
    static Eigen::Matrix<double,3,13> abcsissasAndWeights(){
        Eigen::Matrix<double,13,3> U;

        const double w0 = -0.5*0.149570044467682;

        const double w1 = 0.5*0.175615257433208;
        const double a1 = 0.479308067841920;
        const double b1 = 0.260345966079040;

        const double w2 = 0.5*0.053347235608838;
        const double a2 = 0.869739794195568;
        const double b2 = 0.065130102902216;

        const double w3 = 0.5*0.077113760890257;
        const double a3 = 0.048690315425316;
        const double b3 = 0.312865496004874;
        const double c3 = 0.638444188569810;

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

        // (a3,b3,c3) six permutations
        b3, c3, w3,
        c3, b3, w3,
        a3, c3, w3,
        c3, a3, w3,
        a3, b3, w3,
        b3, a3, w3;

        return U.transpose();
    }
};

} // namespace model

#endif