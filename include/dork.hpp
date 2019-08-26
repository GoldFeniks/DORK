#pragma once

#include "solver.hpp"
#include "coefficients.hpp"

namespace dork {

    template<typename T>
    static runge_kutta_method<T, 2, 2> rkf1 = runge_kutta_method(coefficients::rkf1<T>);

    template<typename T>
    static runge_kutta_method<T, 2, 2> rkf1b = runge_kutta_method(coefficients::rkf1b<T>);

    template<typename T>
    static runge_kutta_method<T, 2> rk2 = runge_kutta_method(coefficients::rk2<T>);

    template<typename T>
    static runge_kutta_method<T, 4> rk4 = runge_kutta_method(coefficients::rk4<T>);

    template<typename T>
    static runge_kutta_method<T, 4, 0, 3> rk4_do3 = runge_kutta_method(coefficients::rk4_do3<T>);

    template<typename T>
    static runge_kutta_method<T, 7, 4> dopri5 = runge_kutta_method(coefficients::dopri5<T>);

    template<typename T>
    static runge_kutta_method<T, 7, 4, 5> dopri5_do5 = runge_kutta_method(coefficients::dopri5_do5<T>);

    template<typename T>
    static runge_kutta_method<T, 6, 4> fehlberg4 = runge_kutta_method(coefficients::fehlberg4<T>);

    template<typename T>
    static runge_kutta_method<T, 13, 7> fehlberg7 = runge_kutta_method(coefficients::fehlberg7<T>);

    template<typename T>
    static runge_kutta_method<T, 12, 6> dop853 = runge_kutta_method(coefficients::dop853<T>);

    template<typename T>
    static runge_kutta_method<T, 12, 6, 7, 4> dop853_do7 = runge_kutta_method(coefficients::dop853_do7<T>);

    template<typename T>
    static runge_kutta_method<T, 8, 5> dverk6 = runge_kutta_method(coefficients::dverk6<T>);

    template<typename T>
    static runge_kutta_method<T, 4, 2> ceschino2 = runge_kutta_method(coefficients::ceschino2<T>);

    template<typename T>
    static runge_kutta_method<T, 5, 3> merson4 = runge_kutta_method(coefficients::merson4<T>);

    template<typename T>
    static runge_kutta_method<T, 5, 3> zonneveld4 = runge_kutta_method(coefficients::zonneveld4<T>);

    template<typename T>
    static runge_kutta_method<T, 3, 2> rkf2 = runge_kutta_method(coefficients::rkf2<T>);

    template<typename T>
    static runge_kutta_method<T, 4, 2> rkf2b = runge_kutta_method(coefficients::rkf2b<T>);

}// namespace dork
