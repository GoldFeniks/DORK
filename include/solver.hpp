#pragma once

#include <tuple>
#include <utility>
#include <type_traits>
#include "utils/utils.hpp"
#include "coefficients.hpp"

namespace dork {

    namespace implementation {

        template<size_t N, size_t M>
        struct operations_implementation {

            operations_implementation() = delete;

            template<typename A>
            static inline void set_zero(A& a) {
                std::get<N>(a) = std::tuple_element_t<N, A>(0);
                operations_implementation<N + 1, M>::set_zero(a);
            }

            template<typename A, typename B, typename C>
            static inline void add_multiplied(const A& a, const B& b, C& c) {
                std::get<N>(c) += b * std::get<N>(a);
                operations_implementation<N + 1, M>::add_multiplied(a, b, c);
            }

            template<typename A, typename B, typename C>
            static inline void set_multiplied(const A& a, const B& b, C& c) {
                std::get<N>(c) = std::get<N>(a) + b * std::get<N>(c);
                operations_implementation<N + 1, M>::set_multiplied(a, b, c);
            }

            template<typename F, typename A, typename R>
            static inline void save_call(const F& functions, const A& arguments, R& result) {
                std::get<N>(result) = std::apply(std::get<N>(functions), arguments);
                operations_implementation<N + 1, M>::save_call(functions, arguments, result);
            }

            template<typename A, typename B>
            static inline void subtract(A& a, const B& b) {
                std::get<N>(a) -= std::get<N>(b);
                operations_implementation<N + 1, M>::subtract(a, b);
            }

            template<typename A, typename B, typename C>
            static inline void add_multiplied_map(const A& a, const B& b, C& c) {
                operations_implementation<0, std::tuple_size_v<C>>::add_multiplied(std::get<N>(a), std::get<N>(b), c);
                operations_implementation<N + 1, M>::add_multiplied_map(a, b, c);
            }

            template<typename A, typename B>
            static inline void set_powers(const A& a, B& b) {
                std::get<N>(b) = a;
                operations_implementation<N + 1, M>::set_powers(a * a, b);
            }

            template<typename A, typename B, typename C>
            static inline void dot(const A& a, const B& b, C& c) {
                c += std::get<N>(a) * std::get<N>(b);
                operations_implementation<N + 1, M>::dot(a, b, c);
            }

            template<typename A, typename B, typename C>
            static inline void err(const A& a0, const A& a1, const A& a2, const B& b0, const B& b1, C& c) {
                c += std::pow((std::get<N>(a1) - std::get<N>(a2)) / (b0 + b1 * std::max(std::abs(std::get<N>(a0)), std::abs(std::get<N>(a1)))), 2);
                operations_implementation<N + 1, M>::err(a0, a1, a2, b0, b1, c);
            }

        };

        template<size_t N>
        struct operations_implementation<N, N> {

            operations_implementation() = delete;

            template<typename A>
            static inline void set_zero(A&) {}

            template<typename A, typename B, typename C>
            static inline void add_multiplied(const A&, const B&, C&) {}

            template<typename A, typename B, typename C>
            static inline void set_multiplied(const A&, const B&, C&) {}

            template<typename F, typename A, typename R>
            static inline void save_call(const F&, const A&, R&) {}

            template<typename A, typename B>
            static inline void subtract(A&, const B&) {}

            template<typename A, typename B, typename C>
            static inline void add_multiplied_map(const A&, const B&, C&) {}

            template<typename A, typename B>
            static inline void set_powers(const A&, B&) {}

            template<typename A, typename B, typename C>
            static inline void dot(const A&, const B&, C&) {}

            template<typename A, typename B, typename C>
            static inline void err(const A&, const A&, const A&, const B&, const B&, C&) {}

        };

        template<size_t N>
        using operations = operations_implementation<0, N>;

        template<typename A, typename V, typename P>
        struct arguments_wrapper {

            arguments_wrapper() : arguments(utils::make_const_reference_tuple(x, values, parameters)) {}

            A x;
            V values;
            P parameters;

            utils::const_reference_tuple_t<A, V, P> arguments;

        };

        template<size_t N, size_t M>
        struct stepper {

            template<typename T, typename V, typename F, typename P>
            static inline void step(types::array1d_t<V, M>& k, const F& functions, const V& current, const T& x, const T& h,
                                    const types::array1d_t<T, M>& c, const types::array2d_t<T, M>& coefficients,
                                    arguments_wrapper<T, V, P>& arguments) {
                arguments.x = x + h * std::get<N>(c);
                step(k, std::get<N>(k), functions, current, h, std::get<N>(coefficients), arguments);
                stepper<N + 1, M>::step(k, functions, current, x, h, c, coefficients, arguments);
            }

            template<typename T, typename V, typename F, typename P>
            static inline void step(types::array1d_t<V, M>& k, V& nk, const F& functions, const V& current, const T& h,
                                    const types::array1d_t<T, M>& coefficients, arguments_wrapper<T, V, P>& arguments) {
                operations<std::tuple_size_v<V>>::set_zero(arguments.values);
                operations<N>::add_multiplied_map(k, coefficients, arguments.values);
                operations<std::tuple_size_v<V>>::set_multiplied(current, h, arguments.values);
                operations<std::tuple_size_v<V>>::save_call(functions, arguments.arguments, nk);
            }

            template<typename T, typename V, size_t K>
            static inline void step(const types::array1d_t<V, M>& k, const types::array1d_t<T, K>& powers,
                    const types::array_t<T, M, K>& coefficients, V& values) {
                T buffer = T(0);
                operations<K>::dot(powers, std::get<N>(coefficients), buffer);
                operations<std::tuple_size_v<V>>::add_multiplied(std::get<N>(k), buffer, values);
                stepper<N + 1, M>::step(k, powers, coefficients, values);
            }

        };

        template<size_t N>
        struct stepper<N, N> {

            template<typename T, typename V, typename F, typename P>
            static inline void step(types::array1d_t<V, N>&, const F&, const V&, const T&, const T&,
                                    const types::array1d_t<T, N>&, const types::array2d_t<T, N>&,
                                    arguments_wrapper<T, V, P>&) {}

            template<typename T, typename V, typename F, typename P>
            static inline void step(types::array1d_t<V, N>&, V&, const F&, const V&, const T&,
                                    const types::array1d_t<T, N>&, arguments_wrapper<T, V, P>&) {}

            template<typename T, typename V, size_t K>
            static inline void step(const types::array1d_t<V, N>&, const types::array1d_t<T, K>&,
                                    const types::array_t<T, N, K>&, V&) {}

        };

        template<typename T>
        class solver_iterator : std::input_iterator_tag {

        public:

            bool operator==(const solver_iterator& other) const {
                return _done == other._done && _owner == other._owner;
            }

            bool operator!=(const solver_iterator& other) const {
                return !(*this == other);
            }

            const typename T::return_value_t& operator*() const {
                return _owner->current();
            }

            const typename T::return_value_t* operator->() const {
                return &_owner->current();
            }

            solver_iterator& operator++() {
                if (!(_done = _owner->is_done()))
                    _owner->next();
                return *this;
            }

            solver_iterator operator++(int) {
                const auto res = *this;
                ++(*this);
                return res;
            }

        private:

            bool _done = false;
            T* _owner = nullptr;

            friend T;

            solver_iterator(T* owner, const bool done) : _owner(owner), _done(done) {}
            explicit solver_iterator(T* owner) : solver_iterator(owner, owner->is_done()) {}

        };

        template<typename B, typename F, typename V, typename P, typename = std::enable_if_t<std::tuple_size_v<F> == std::tuple_size_v<V>>>
        class solver {

        private:

            using coefficients_t = typename B::coefficients_t;
            static constexpr size_t n = coefficients_t::n;
            static constexpr size_t m = coefficients_t::m;
            static constexpr size_t k = coefficients_t::k;
            static constexpr bool embedded = B::embedded;

        public:

            using arg_t = typename B::arg_t;
            using iterator = solver_iterator<solver>;
            using return_value_t = utils::const_reference_tuple_t<arg_t, types::vector1d_t<V>>;

            static constexpr size_t steps = n + k;
            static constexpr size_t count = std::tuple_size_v<F>;

            solver(B base, F functions, types::vector1d_t<V> initial_values, types::vector1d_t<P> parameters) :
                    _base(std::move(base)), _functions(std::move(functions)), _initial_values(initial_values),
                    _buff_values(initial_values.size()), _last_values(initial_values.size()), _current_values(std::move(initial_values)),
                    _k(_current_values.size()), _parameters(std::move(parameters)),
                    _last_value(utils::make_const_reference_tuple(_x, _last_values)),
                    _current_value(utils::make_const_reference_tuple(_base.x(), _current_values)) {}

            const return_value_t& next() {
                _x = std::get<0>(_current_value);
                std::swap(_last_values, _current_values);
                start:
                _h = _base.h();
                for (size_t i = 0; i < _last_values.size(); ++i) {
                    _arguments.parameters = _parameters[i];
                    stepper<0, steps>::step(_k[i], _functions, _last_values[i], _x, _h,
                                            _base.coefficients().c_table, _base.coefficients().a_table, _arguments);
                    operations<count>::set_zero(_current_values[i]);
                    operations<n>::add_multiplied_map(_k[i], _base.coefficients().b_table[0], _current_values[i]);
                    operations<count>::set_multiplied(_last_values[i], _h, _current_values[i]);
                    if constexpr (embedded) {
                        operations<count>::set_zero(_buff_values[i]);
                        operations<n>::add_multiplied_map(_k[i], _base.coefficients().b_table[1], _buff_values[i]);
                        operations<count>::set_multiplied(_last_values[i], _h, _buff_values[i]);
                    }
                }
                if constexpr (embedded)
                    if (!_base.accept_step(_initial_values, _current_values, _buff_values))
                        goto start;
                _base.next();
                return _current_value;
            }

            utils::collect_types_t<arg_t, types::vector1d_t<V>> value(const arg_t& x) {
                if constexpr (m > 0) {
                    types::array1d_t<arg_t, m> powers;
                    operations<m>::set_powers(x, powers);
                    types::vector1d_t<V> res(_k.size());
                    for (size_t i = 0; i < res.size(); ++i) {
                        operations<count>::set_zero(res[i]);
                        stepper<0, steps>::step(_k[i], powers, _base.coefficients().d_table, res[i]);
                        operations<count>::set_multiplied(_last_values[i], _h, res[i]);
                    }
                    return std::make_tuple(_x + _h * x, res);
                } else {
                    auto res = _current_values;
                    for (size_t i = 0; i < res.size(); ++i) {
                        operations<count>::subtract(res[i], _last_values[i]);
                        operations<count>::set_multiplied(_last_values[i], x, res[i]);
                    }
                    return std::make_tuple(_x + _h * x, res);
                }
            }

            utils::collect_types_t<arg_t, types::vector1d_t<V>> operator()(const arg_t& x) {
                return value(x);
            }

            const return_value_t& last() const {
                return _last_value;
            }

            const return_value_t& current() const {
                return _current_value;
            }

            [[nodiscard]] bool is_done() const {
                return _base.is_done();
            }

            explicit operator bool() const {
                return !is_done();
            }

            iterator begin() {
                return iterator(this);
            }

            iterator end() {
                return iterator(this, true);
            }

        private:

            B _base;
            arg_t _x, _h;
            const F _functions;
            types::vector1d_t<P> _parameters;
            arguments_wrapper<arg_t, V, P> _arguments;
            const types::vector1d_t<V> _initial_values;
            std::conditional_t<embedded, types::vector1d_t<V>, size_t> _buff_values;
            types::vector1d_t<V> _last_values, _current_values;
            types::vector1d_t<types::array1d_t<V, steps>> _k;
            const return_value_t _last_value, _current_value;

        };

        template<typename B, typename F, typename V>
        class solver<B, F, V, std::tuple<>> {

        private:

            using coefficients_t = typename B::coefficients_t;
            static constexpr size_t n = coefficients_t::n;
            static constexpr size_t m = coefficients_t::m;
            static constexpr size_t k = coefficients_t::k;
            static constexpr bool embedded = B::embedded;

        public:

            using arg_t = typename B::arg_t;
            using iterator = solver_iterator<solver>;
            using return_value_t = utils::const_reference_tuple_t<arg_t, V>;

            static constexpr size_t steps = n + k;
            static constexpr size_t count = std::tuple_size_v<F>;

            solver(B base, F functions, V initial_values) :
                    _base(std::move(base)), _functions(std::move(functions)), _initial_values(initial_values),
                    _current_values(std::move(initial_values)),
                    _last_value(utils::make_const_reference_tuple(_x, _last_values)),
                    _current_value(utils::make_const_reference_tuple(_base.x(), _current_values)) {}

            const return_value_t& next() {
                _x = std::get<0>(_current_value);
                _last_values = std::move(_current_values);
                start:
                _h = _base.h();
                stepper<0, steps>::step(_k, _functions, _last_values, _x, _h,
                                        _base.coefficients().c_table, _base.coefficients().a_table, _arguments);
                operations<count>::set_zero(_current_values);
                operations<n>::add_multiplied_map(_k, _base.coefficients().b_table[0], _current_values);
                operations<count>::set_multiplied(_last_values, _h, _current_values);
                if constexpr (embedded) {
                    operations<count>::set_zero(_buff_values);
                    operations<n>::add_multiplied_map(_k, _base.coefficients().b_table[1], _buff_values);
                    operations<count>::set_multiplied(_last_values, _h, _buff_values);
                    if (!_base.accept_step(_initial_values, _current_values, _buff_values))
                        goto start;
                }
                _base.next();
                return _current_value;
            }

            utils::collect_types_t<arg_t, V> value(const arg_t& x) {
                if constexpr (m > 0) {
                    types::array1d_t<arg_t, m> powers;
                    operations<m>::set_powers(x, powers);
                    V res;
                    operations<count>::set_zero(res);
                    stepper<0, steps>::step(_k, powers, _base.coefficients().d_table, res);
                    operations<count>::set_multiplied(_last_values, _h, res);
                    return std::tuple_cat(std::make_tuple(_x + _h * x), res);
                } else {
                    auto res = _current_values;
                    operations<count>::subtract(res, _last_values);
                    operations<count>::set_multiplied(_last_values, x, res);
                    return std::tuple_cat(std::make_tuple(_x + _h * x), res);
                }
            }

            utils::collect_types_t<arg_t, V> operator()(const arg_t& x) {
                return value(x);
            }

            const return_value_t& last() const {
                return _last_value;
            }

            const return_value_t& current() const {
                return _current_value;
            }

            [[nodiscard]] bool is_done() const {
                return _base.is_done();
            }

            explicit operator bool() const {
                return !is_done();
            }

            iterator begin() {
                return iterator(this);
            }

            iterator end() {
                return iterator(this, true);
            }

        private:

            B _base;
            arg_t _x, _h;
            const F _functions;
            types::array1d_t<V, steps> _k;
            const V _initial_values;
            std::conditional_t<embedded, V, types::array1d_t<arg_t, 0>> _buff_values;
            V _last_values, _current_values;
            const return_value_t _last_value, _current_value;
            arguments_wrapper<arg_t, V, std::tuple<>> _arguments;

        };

        template<typename B, typename F>
        class functions {

        public:

            functions(B base, F functions) :
                    _base(std::move(base)), _functions(std::move(functions)) {}

            template<typename... V, typename... P>
            solver<B, F, std::tuple<V...>, std::tuple<P...>> initial_values(
                    const types::vector1d_t<std::tuple<V...>>& values, const types::vector1d_t<std::tuple<P...>>& parameters) const {
                return solver<B, F, std::tuple<V...>, std::tuple<P...>>(_base, _functions, values, parameters);
            }

            template<typename... V>
            solver<B, F, std::tuple<V...>, std::tuple<>> initial_values(const V&... values) const {
                return solver<B, F, std::tuple<V...>, std::tuple<>>(_base, _functions, std::make_tuple(values...));
            }

            template<typename... V, typename... P>
            solver<B, F, std::tuple<V...>, std::tuple<P...>> operator()(
                    const types::vector1d_t<std::tuple<V...>>& values, const types::vector1d_t<std::tuple<P...>>& parameters) const {
                return initial_values(values, parameters);
            }

            template<typename... V>
            solver<B, F, std::tuple<V...>, std::tuple<>> operator()(const V&... values) const {
                return initial_values(values...);
            }

        private:

            B _base;
            const F _functions;

        };

        template<typename D, typename T, size_t N, size_t Q = 0, size_t M = 0, size_t K = 0>
        class solver_base {

        public:

            using arg_t = T;
            using coefficients_t = dork::coefficients::coefficients<T, N, Q, M, K>;

            static constexpr bool embedded = false;

            explicit solver_base(const coefficients_t& coefficients) : _coefficients(coefficients) {}
            solver_base(const T& x, const coefficients_t& coefficients) : _x(x), _coefficients(coefficients) {}
            solver_base(const T& x, const T& h, const coefficients_t& coefficients) : _x(x), _h(h), _coefficients(coefficients) {}

            D& as_true_type() {
                return *static_cast<D*>(this);
            };

            const D& as_true_type() const {
                return *static_cast<const D*>(this);
            };

            template<typename... F>
            implementation::functions<D, std::tuple<const F...>> functions(const F&... functions) const {
                return implementation::functions<D, std::tuple<const F...>>(as_true_type(), std::make_tuple(functions...));
            }

            template<typename... F>
            implementation::functions<D, std::tuple<const F...>> operator()(const F&... functions) const {
                return this->functions<F...>(functions...);
            }

            const T& x() const {
                return _x;
            }

            const T& h() const {
                return _h;
            }

            const coefficients_t& coefficients() const {
                return _coefficients;
            }

        protected:

            T _x, _h;
            const coefficients_t& _coefficients;

        };

        template<typename T, size_t N, size_t Q = 0, size_t M = 0, size_t K = 0>
        class uniform_solver_base : public solver_base<uniform_solver_base<T, N, Q, M, K>, T, N, Q, M, K> {

        private:

            using base_t = solver_base<uniform_solver_base<T, N, Q, M, K>, T, N, Q, M, K>;

        public:

            using typename base_t::coefficients_t;

            uniform_solver_base(const T& x0, const T& x1, size_t n, const coefficients_t& coefficients)
                    : base_t(x0, (x1 - x0) / (n - 1), coefficients), _n(n) {}

            void next() {
                _x += _h;
                _index += 1;
            }

            [[nodiscard]] bool is_done() const {
                return _index >= _n;
            }

        private:

            using base_t::_x;
            using base_t::_h;
            const size_t _n;
            size_t _index = 1;

        };

        template<typename T, size_t N, size_t Q = 0, size_t M = 0, size_t K = 0>
        class vector_solver_base : public solver_base<vector_solver_base<T, N, Q, M, K>, T, N, Q, M, K>{

        private:

            using base_t = solver_base<vector_solver_base<T, N, Q, M, K>, T, N, Q, M, K>;

        public:

            using typename base_t::coefficients_t;

            vector_solver_base(types::vector1d_t<T> xs, const coefficients_t& coefficients)
                    : base_t(xs[0], xs[1] - xs[0], coefficients), _xs(std::move(xs)) {}

            void next() {
                _index += 1;
                _x = _xs[_index - 1];
                if (!is_done())
                    _h = _xs[_index] - _x;
            }

            [[nodiscard]] bool is_done() const {
                return _index >= _xs.size();
            }

        private:

            using base_t::_x;
            using base_t::_h;
            size_t _index = 1;
            const types::vector1d_t<T> _xs;

        };

        template<typename It, typename T, size_t N, size_t Q = 0, size_t M = 0, size_t K = 0>
        class iterator_solver_base : public solver_base<iterator_solver_base<It, T, N, Q, M, K>, T, N, Q, M, K> {

        private:

            using base_t = solver_base<iterator_solver_base<It, T, N, Q, M, K>, T, N, Q, M, K>;

        public:

            using typename base_t::coefficients_t;

            iterator_solver_base(It begin, It end, const coefficients_t& coefficients)
                    : base_t(*begin, coefficients), _current(std::move(begin)), _x1(*++_current), _last(end) {
                _h = _x1 - _x;
            }

            void next() {
                _x = _x1;
                ++_current;
                if (!is_done()) {
                    _x1 = *_current;
                    _h = _x1 - _x;
                }
            }

            [[nodiscard]] bool is_done() const {
                return _current == _last;
            }

        private:

            using base_t::_x;
            using base_t::_h;
            It _current;
            T _x1;
            const It _last;

        };

        template<typename T, size_t N, size_t Q, size_t M = 0, size_t K = 0>
        class step_control_solver_base : public solver_base<step_control_solver_base<T, N, Q, M, K>, T, N, Q, M, K> {

        private:

            using base_t = solver_base<step_control_solver_base<T, N, Q, M, K>, T, N, Q, M, K>;
            static constexpr T eps = T(1e-10);

        public:

            using typename base_t::coefficients_t;
            static constexpr bool embedded = true;

            step_control_solver_base(const T& x0, const T& x1, const T& h, const T& atol, const T& rtol,
                    const T& fac, const T& fac_min, const T& fac_max, const coefficients_t& coefficients)
                : base_t(x0, h, coefficients), _x1(x1), _fac(fac), _fac_min(fac_min), _fac_max(fac_max), _atol(atol), _rtol(rtol) {
                    static_assert(Q > 0, "Method does not support step control");
                }

            template<typename V>
            bool accept_step(const types::vector1d_t<V>& y0, const types::vector1d_t<V>& y1, const types::vector1d_t<V>& y2) {
                _lh = _h;
                T error;
                std::tie(error, _h) = _calc_step_error(y0.front(), y1.front(), y2.front());
                for (size_t i = 1; i < y0.size(); ++i) {
                    const auto& [err, h] = _calc_step_error(y0[i], y1[i], y2[i]);
                    error = std::max(error, err);
                    _h = std::min(_h, h);
                }
                return _check_step(error);
            }

            template<typename V>
            bool accept_step(const V& y0, const V& y1, const V& y2) {
                _lh = _h;
                T error;
                std::tie(error, _h) = _calc_step_error(y0, y1, y2);
                return _check_step(error);
            }

            void next() {
                _x += _lh;
            }

            [[nodiscard]] bool is_done() const {
                return _x1 - _x < eps;
            }

        private:

            using base_t::_x;
            using base_t::_h;
            T _lh;
            const T _x1, _fac, _fac_min, _fac_max, _atol, _rtol;

            template<typename V>
            std::tuple<T, T> _calc_step_error(const V& y0, const V& y1, const V& y2) {
                T err(0);
                operations<std::tuple_size_v<V>>::err(y0, y1, y2, _atol, _rtol, err);
                err = std::sqrt(err / std::tuple_size_v<V>);
                return { err, _h * std::min(_fac_max, std::max(_fac_min, _fac * std::pow(1 / err, T(1) / T(Q + 1)))) };
            }

            bool _check_step(const T& error) {
                if (T(1) - error > -eps) {
                    _h = std::min(std::abs(_x1 - _x - _lh), _h);
                    return true;
                }
                _h = std::min(std::abs(_x1 - _x), _h);
                return false;
            }

        };

        template<typename T, size_t N, size_t Q = 0, size_t M = 0, size_t K = 0>
        class endless_solver_base : public solver_base<endless_solver_base<T, N, Q, M, K>, T, N, Q, M, K> {

        private:

            using base_t = solver_base<endless_solver_base<T, N, Q, M, K>, T, N, Q, M, K>;

        public:

            using typename base_t::coefficients_t;

            endless_solver_base(const T& x0, const T& h, const coefficients_t& coefficients)
                    : base_t(x0, h, coefficients) {}

            void next() {
                _x += _h;
            }

            [[nodiscard]] bool is_done() const {
                return false;
            }

        private:

            using base_t::_x;
            using base_t::_h;

        };

    }// namespace implementation

    template<typename T, size_t N, size_t Q = 0, size_t M = 0, size_t K = 0>
    class runge_kutta_method {

        public:

            using coefficients_t = coefficients::coefficients<T, N, Q, M, K>;

            explicit runge_kutta_method(const coefficients_t& coefficients) : _coefficients(coefficients) {}

            implementation::uniform_solver_base<T, N, Q, M, K> solver(const T& x0, const T& x1, const size_t& n) {
                return implementation::uniform_solver_base<T, N, Q, M, K>(x0, x1, n, _coefficients);
            }

            implementation::uniform_solver_base<T, N, Q, M, K> operator()(const T& x0, const T& x1, const size_t& n) {
                return solver(x0, x1, n);
            }

            implementation::vector_solver_base<T, N, Q, M, K> solver(const types::vector1d_t<T>& xs) {
                return implementation::vector_solver_base<T, N, Q, M, K>(xs, _coefficients);
            }

            implementation::vector_solver_base<T, N, Q, M, K> operator()(const types::vector1d_t<T>& xs) {
                return solver(xs);
            }

            template<typename It>
            implementation::iterator_solver_base<It, T, N, Q, M, K> solver(It begin, It end) {
                return implementation::iterator_solver_base<It, T, N, Q, M, K>(begin, end, _coefficients);
            }

            template<typename It>
            implementation::iterator_solver_base<It, T, N, Q, M, K> operator()(It begin, It end) {
                return solver(begin, end);
            }

            implementation::step_control_solver_base<T, N, Q, M, K> solver(
                    const T& x0, const T& x1, const T& h, const T& atol, const T& rtol,
                    const T& fac = std::pow(T(1) / T(4), T(1) / T(Q + 1)),
                    const T& fac_min = T(1) / T(2),
                    const T& fac_max = T(3) / T(2)
            ) {
                return implementation::step_control_solver_base<T, N, Q, M, K>(x0, x1, h, atol, rtol, fac, fac_min, fac_max, _coefficients);
            }

            implementation::step_control_solver_base<T, N, Q, M, K> operator()(
                    const T& x0, const T& x1, const T& h, const T& atol, const T& rtol,
                    const T& fac = std::pow(T(1) / T(4), T(1) / T(Q + 1)),
                    const T& fac_min = T(1) / T(2),
                    const T& fac_max = T(3) / T(2)
            ) {
                return solver(x0, x1, h, atol, rtol, fac, fac_min, fac_max);
            }

            implementation::endless_solver_base<T, N, Q, M, K> solver(const T& x0, const T& h) {
                return implementation::endless_solver_base<T, N, Q, M, K>(x0, h, _coefficients);
            }

            implementation::endless_solver_base<T, N, Q, M, K> operator()(const T& x0, const T& h) {
                return solver(x0, h);
            }

        private:

            const coefficients_t& _coefficients;

    };

}// namespace dork
