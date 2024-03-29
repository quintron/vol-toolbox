#include "quad.h"

namespace voltlbx
{

    #define M_SQRT2 1.414213562373095048801689
    #define M_INV_SQRT_PI 0.5641895835477562869480795

    template<>
    const double GaussHermite<2>::points[] =
    {
            -0.7071067811865475 * M_SQRT2,
            0.7071067811865475 * M_SQRT2,
    };

    template<>
    const double GaussHermite<2>::weights[] =
    {
            0.8862269254527579 * M_INV_SQRT_PI,
            0.8862269254527579 * M_INV_SQRT_PI,

    };

    template<>
    const double GaussHermite<3>::points[] =
    {
            -1.224744871391589 * M_SQRT2,
            0.0 * M_SQRT2,
            1.224744871391589 * M_SQRT2,
    };

    template<>
    const double GaussHermite<3>::weights[] =
    {
            0.2954089751509194 * M_INV_SQRT_PI,
            1.1816359006036772 * M_INV_SQRT_PI,
            0.2954089751509194 * M_INV_SQRT_PI,
    };

    template<>
    const double GaussHermite<4>::points[] = {
            -1.6506801238857845565883 * M_SQRT2,
            -0.5246476232752903178841 * M_SQRT2,
            +0.5246476232752903178841 * M_SQRT2,
            +1.650680123885784555883 * M_SQRT2,

    };

    template<>
    const double GaussHermite<4>::weights[] = {
            0.08131283544724517714304 * M_INV_SQRT_PI,
            0.8049140900055128365061 * M_INV_SQRT_PI,
            0.8049140900055128365061 * M_INV_SQRT_PI,
            0.08131283544724517714304 * M_INV_SQRT_PI,

    };

    template<>
    const double GaussHermite<5>::points[] = {
            -2.0201828704560856 * M_SQRT2,
            -0.9585724646138185 * M_SQRT2,
            0.0 * M_SQRT2,
            0.9585724646138185 * M_SQRT2,
            2.0201828704560856 * M_SQRT2,

    };

    template<>
    const double GaussHermite<5>::weights[] =

    {
            0.019953242059045917 * M_INV_SQRT_PI,
            0.3936193231522411 * M_INV_SQRT_PI,
            0.9453087204829418 * M_INV_SQRT_PI,
            0.3936193231522411 * M_INV_SQRT_PI,
            0.019953242059045917 * M_INV_SQRT_PI,

    };

    template<>
    const double GaussHermite<8>::points[] = {
            -2.9306374202572440192235027 * M_SQRT2,
            -1.9816567566958429258546306 * M_SQRT2,
            -1.1571937124467801947207658 * M_SQRT2,
            -0.38118699020732211685471889 * M_SQRT2,
            +0.38118699020732211685471889 * M_SQRT2,
            +1.1571937124467801947207658 * M_SQRT2,
            +1.9816567566958429258546306 * M_SQRT2,
            +2.9306374202572440192235027 * M_SQRT2,

    };

    template<>
    const double GaussHermite<8>::weights[] = {
            1.9960407221136761920609045E-4 * M_INV_SQRT_PI,
            0.017077983007413475456203056 * M_INV_SQRT_PI,
            0.2078023258148918795432586 * M_INV_SQRT_PI,
            0.66114701255824129103041597 * M_INV_SQRT_PI,
            0.66114701255824129103041597 * M_INV_SQRT_PI,
            0.207802325814891879543256862 * M_INV_SQRT_PI,
            0.017077983007413475456203056 * M_INV_SQRT_PI,
            1.996040722113676192060905E-4 * M_INV_SQRT_PI,

    };

    template<>
    const double GaussHermite<10>::points[] = {
            -3.4361591188377374 * M_SQRT2,
            -2.5327316742327897 * M_SQRT2,
            -1.7566836492998816 * M_SQRT2,
            -1.0366108297895136 * M_SQRT2,
            -0.3429013272237046 * M_SQRT2,
            0.3429013272237046 * M_SQRT2,
            1.0366108297895136 * M_SQRT2,
            1.7566836492998816 * M_SQRT2,
            2.5327316742327897 * M_SQRT2,
            3.4361591188377374 * M_SQRT2,
    };

    template<>
    const double GaussHermite<10>::weights[] = {
            7.640432855232641e-06 * M_INV_SQRT_PI,
            0.0013436457467812324 * M_INV_SQRT_PI,
            0.033874394455481106 * M_INV_SQRT_PI,
            0.2401386110823147 * M_INV_SQRT_PI,
            0.6108626337353258 * M_INV_SQRT_PI,
            0.6108626337353258 * M_INV_SQRT_PI,
            0.2401386110823147 * M_INV_SQRT_PI,
            0.033874394455481106 * M_INV_SQRT_PI,
            0.0013436457467812324 * M_INV_SQRT_PI,
            7.640432855232641e-06 * M_INV_SQRT_PI,
    };

    template<>
    const double GaussHermite<12>::points[] = {
            -3.889724897869782 * M_SQRT2,
            -3.0206370251208896 * M_SQRT2,
            -2.279507080501098 * M_SQRT2,
            -1.5976826351526048 * M_SQRT2,
            -0.9477883912401638 * M_SQRT2,
            -0.31424037625435913 * M_SQRT2,
            0.31424037625435913 * M_SQRT2,
            0.9477883912401638 * M_SQRT2,
            1.5976826351526048 * M_SQRT2,
            2.2795070805010598 * M_SQRT2,
            3.0206370251208896 * M_SQRT2,
            3.889724897869782 * M_SQRT2,
    };

    template<>
    const double GaussHermite<12>::weights[] = {
            2.6585516843563044e-07 * M_INV_SQRT_PI,
            8.573687043587868e-05 * M_INV_SQRT_PI,
            .00390639058462906 * M_INV_SQRT_PI,
            .05160798561588398 * M_INV_SQRT_PI,
            .2604923102641611 * M_INV_SQRT_PI,
            .5701352362624795 * M_INV_SQRT_PI,
            .5701352362624795 * M_INV_SQRT_PI,
            .2604923102641611 * M_INV_SQRT_PI,
            .05160798561588398 * M_INV_SQRT_PI,
            .00390639058462906 * M_INV_SQRT_PI,
            8.573687043587868e-05 * M_INV_SQRT_PI,
            2.6585516843563044e-07 * M_INV_SQRT_PI,
    };

    template<>
    const double GaussHermite<16>::points[] =

    {
            -4.6887389393058183646884986 * M_SQRT2,
            -3.8694479048601226987194241 * M_SQRT2,
            -3.1769991619799560268139946 * M_SQRT2,
            -2.5462021578474813621593287 * M_SQRT2,
            -1.9517879909162539774346554 * M_SQRT2,
            -1.3802585391988807963720897 * M_SQRT2,
            -0.8229514491446558925824545 * M_SQRT2,
            -0.2734810461381524521582804 * M_SQRT2,
            +0.2734810461381524521582804 * M_SQRT2,
            +0.8229514491446558925824545 * M_SQRT2,
            +1.3802585391988807963720897 * M_SQRT2,
            +1.9517879909162539774346554 * M_SQRT2,
            +2.5462021578474813621593287 * M_SQRT2,
            +3.1769991619799560268139946 * M_SQRT2,
            +3.8694479048601226987194241 * M_SQRT2,
            +4.6887389393058183646884986 * M_SQRT2,
    };

    template<>
    const double GaussHermite<16>::weights[] =
    {
            2.654807474011182244709264E-10 * M_INV_SQRT_PI,
            3.320980844865210653387494E-7 * M_INV_SQRT_PI,
            2.711860092537881512018914E-5 * M_INV_SQRT_PI,
            9.322840086241805299142773E-4 * M_INV_SQRT_PI,
            0.0128803115355099736834643 * M_INV_SQRT_PI,
            0.08381004139898582941542073 * M_INV_SQRT_PI,
            0.28064745852853367536946334 * M_INV_SQRT_PI,
            0.50792947901661374191351734 * M_INV_SQRT_PI,
            0.50792947901661374191351734 * M_INV_SQRT_PI,
            0.28064745852853367536946334 * M_INV_SQRT_PI,
            0.08381004139898582941542073 * M_INV_SQRT_PI,
            0.0128803115355099736834643 * M_INV_SQRT_PI,
            9.322840086241805299142773E-4 * M_INV_SQRT_PI,
            2.711860092537881512018914E-5 * M_INV_SQRT_PI,
            2.320980844865210653387494E-7 * M_INV_SQRT_PI,
            2.654807474011182244709264E-10 * M_INV_SQRT_PI,
    };


    template<>
    const double GaussLegendre<2>::points[] =
    {
            -0.5773502691896257,
            0.5773502691896257
    };
    template<>
    const double GaussLegendre<2>::weights[] = {
            1.0,
            1.0
    };


    template<>
    const double GaussLegendre<3>::points[] = {
            -0.7745966692414834,
            0.0,
            0.7745966692414834
    };

    template<>
    const double GaussLegendre<3>::weights[] = {
            0.5555555555555557,
            0.8888888888888888,
            0.5555555555555557

    };


    template<>
    const double GaussLegendre<4>::points[] = {
            -0.8611363115940526,
            -0.3399810435848563,
            +0.3399810435848563,
            +0.8611363115940526
    };
    template<>
    const double GaussLegendre<4>::weights[] =
    {
            +0.3478548451374538,
            +0.6521451548625461 ,
            +0.6521451548625461 ,
            +0.3478548451374538
    };

    template class GaussLegendre <4>;

    template<>
    const double GaussLegendre<8>::points[] =

    {
            -0.9602898564975363,
            -0.7966664774136267,
            -0.5255324099163290,
            -0.1834346424956498,
            +0.1834346424956498,
            +0.5255324099163290,
            +0.7966664774136267,
            +0.9602898564975363

    };

    template<>
    const double GaussLegendre<8>::weights[] =

    {
            +0.1012285362903761 ,
            +0.2223810344533745,
            +0.3137066458778873,
            +0.3626837833783620,
            +0.3626837833783620,
            +0.3137066458778873,
            +0.2223810344533745,
            +0.1012285362903761

    };
    template class GaussLegendre <8>;

    template<>
    const double GaussLegendre<10>::points[] =
    {
            -0.9739065285171717,
            -0.8650633666889845,
            -0.6794095682990244 ,
            -0.4333953941292472,
            -0.14887433898163122,
            0.14887433898163122,
            0.4333953941292472,
            0.6794095682990244,
            0.8650633666889845,
            0.9739065285171717,
    };

    template<>
    const double GaussLegendre<10>::weights[] =
    {
            .06667134430868807,
            .14945134915058036,
            .219086362515982,
            .2692667193099965,
            .295524224714753,
            .295524224714753,
            .2692667193099965,
            .219086362515982,
            .14945134915058036,
            .06667134430868807,
    };

    template class GaussLegendre <10>;

    template<>
    const double GaussLegendre<16>::points[] = {
            -0.9894009349916499,
            -0.9445750230732326,
            -0.8656312023878318,
            -0.7554044083550030,
            -0.6178762444026438,
            -0.4580167776572274,
            -0.2816035507792589,
            -0.9501250983763744E-01,
            +0.9501250983763744E-01,
            +0.28160355077925689,
            +0.4580167776572274,
            +0.6178762444026438,
            +0.7554044083550030,
            +0.8656312023878318,
            +0.9445750230732326,
            +0.9894009349916499
    };

    template<>
    const double GaussLegendre<16>::weights[] = {
            +0.2715245941175407E-01,
            +0.6225352393864787E-01,
            +0.9515851168249277E-01,
            +0.1246289712555339,
            +0.1495959888165767 ,
            +0.1691565193950026,
            +0.1826034150449236,
            +0.1894506104550686 ,
            +0.1894506104550686 ,
            +0.1826034150449236,
            +0.1691565193950026,
            +0.1495959888165767,
            +0.1246289712555339,
            +0.9515851168249277E-01,
            +0.6225352393864787E-01,
            +0.2715245941175407E-01
    };

    template class GaussLegendre<16>;

}
