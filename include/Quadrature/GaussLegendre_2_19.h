#ifndef model_GAUSSLEGENDRE_2_19_H_
#define model_GAUSSLEGENDRE_2_19_H_

namespace model{

template <>
struct GaussLegendre<2,19>{
    // Dunavant p=9, ng=19 (weights sum to 0.5 on reference triangle)
    static Eigen::Matrix<double,3,19> abcsissasAndWeights(){
        Eigen::Matrix<double,19,3> U;

        const double w0 = 0.5*0.097135796282799;

        const double w1 = 0.5*0.031334700227139;
        const double a1 = 0.020634961602525;
        const double b1 = 0.489682519198738;

        const double w2 = 0.5*0.077827541004774;
        const double a2 = 0.125820817014127;
        const double b2 = 0.437089591492937;

        const double w3 = 0.5*0.079647738927210;
        const double a3 = 0.623592928761935;
        const double b3 = 0.188203535619033;

        const double w4 = 0.5*0.025577675658698;
        const double a4 = 0.910540973211095;
        const double b4 = 0.044729513394453;

        const double w5 = 0.5*0.043283539377289;
        const double a5 = 0.036838412054736;
        const double b5 = 0.221962989160766;
        const double c5 = 0.741198598784498;

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

        // (a4,b4,b4)
        b4, b4, w4,
        a4, b4, w4,
        b4, a4, w4,

        // (a5,b5,c5) six permutations
        b5, c5, w5,
        c5, b5, w5,
        a5, c5, w5,
        c5, a5, w5,
        a5, b5, w5,
        b5, a5, w5;

        return U.transpose();
    }
};

} // namespace model

#endif