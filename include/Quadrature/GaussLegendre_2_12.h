#ifndef model_GAUSSLEGENDRE_2_12_H_
#define model_GAUSSLEGENDRE_2_12_H_

namespace model{

template <>
struct GaussLegendre<2,12>{
    // Dunavant p=6, ng=12 (weights sum to 0.5 on reference triangle)
    // Source: Dunavant 1985 tables / Burkardt triangle_dunavant_rule. :contentReference[oaicite:1]{index=1}
    static Eigen::Matrix<double,3,12> abcsissasAndWeights(){
        Eigen::Matrix<double,12,3> U;

        const double w1 = 0.5*0.116786275726379;
        const double a1 = 0.501426509658179;
        const double b1 = 0.249286745170910;

        const double w2 = 0.5*0.050844906370207;
        const double a2 = 0.873821971016996;
        const double b2 = 0.063089014491502;

        const double w3 = 0.5*0.082851075618374;
        const double a3 = 0.053145049844817;
        const double b3 = 0.310352451033784;
        const double c3 = 0.636502499121399;

        U <<
        // (a1,b1,b1)
        b1, b1, w1,
        a1, b1, w1,
        b1, a1, w1,

        // (a2,b2,b2)
        b2, b2, w2,
        a2, b2, w2,
        b2, a2, w2,

        // (a3,b3,c3) six permutations in (u,v)=(beta,gamma)
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