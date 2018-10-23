#pragma once

#include <cstddef>
#include <functional>
#include <vector>
#include <type_traits>
#include <tuple>

namespace ssh {

    namespace utils {

        /**
         * Initializes vector values at indexes in range [start, last) by calling func(index)
         */
        template<typename Vector, typename Function, typename Index = size_t>
        void init_vector(Vector& vector, Index first, Index last, const Function& function) {
            for (Index i = first; i < last; ++i)
                vector[i] = function(i);
        }

        /**
         * Initializes vector values at indexes in range [0, last) by calling func(index)
         */
        template<typename Vector, typename Function, typename Index = size_t>
        void init_vector(Vector& vector, Index last, const Function& function) {
            init_vector<Vector, Function, Index>(vector, Index(0), last, function);
        }

        /**
         * Initializes vector values at indexes in range [0, vector.size()) by calling func(index)
         */
        template<typename Vector, typename Function>
        void init_vector(Vector& vector, const Function& function) {
            init_vector(vector, vector.size(), function);
        }

        /**
         * Creates Vector { a, a + d, a + 2 * d, ... , a + (n - 1) * d }
         */
        template<typename Vector, typename Value = typename Vector::value_type, typename Index = size_t>
        Vector create_mesh_step(const Value& a, const Value& d, const Index n) {
            Vector result(n);
            init_vector(result, [a,d](Index i) { return a + i * d; });
            return result;
        }

        /**
         * Creates Vector { a, a + d, a + 2 * d, ... , a + (n - 1) * d }, where d = (b - a) / (n - 1)
         */
        template<typename Vector, typename Value = typename Vector::value_type, typename Index = size_t>
        Vector create_mesh_bounds(const Value& a, const Value& b, const Index n) {
            return create_mesh_step<Vector>(a, (b - a) / (n - 1), n);
        }

        template<typename T, typename I = size_t>
        T factorial(const I n) {
            return n <= I(1) ? T(1) : factorial<T>(n - I(1)) * T(n);
        }

        template<typename... Args>
        struct arguments {

            using args_t = std::tuple<Args...>;

            const args_t args;
            
            arguments(const Args&... args) : args(std::tuple(args...)) {}

        };

        struct no_arguments {

            using args_t = std::tuple<>;
            
            const args_t args;

        };

        template<typename T, typename Tuple, typename = typename Tuple::args_t>
        struct is_constructible_from_tuple_types;

        template<typename T, typename Tuple, typename... Args>
        struct is_constructible_from_tuple_types<T, Tuple, std::tuple<Args...>> {

            static constexpr bool value = std::is_constructible_v<T, Args...>;

        };

        template<typename T, typename Tuple>
        inline constexpr bool is_constructible_from_tuple_types_v = is_constructible_from_tuple_types<T, Tuple>::value;

        template<typename T>
        class constructor {

            public:

                template<typename... Args>
                static auto construct_from_tuple(const std::tuple<Args...>& args) {
                    return choose_constructor<0, Args...>::type::construct(args);
                }

                template<typename... Args>
                static auto construct(const Args... args) {
                    return construct_from_tuple<Args...>(std::tuple<Args...>(args...));
                }

            private:

                template<size_t N>
                struct _constructor {

                    template<typename... Args>
                    static auto construct(const std::tuple<Args...>& args) {
                        return std::apply([](auto&&... args) { 
                            return T(std::forward<decltype(args)>(args)...); 
                        }, std::get<N>(args).args);
                    }

                };

                template<size_t N, typename TArgs, typename... Tail>
                struct choose_constructor {

                    using type = std::conditional_t<is_constructible_from_tuple_types_v<T, TArgs>,
                                                    _constructor<N>, typename choose_constructor<N + 1, Tail...>::type>;

                };

                template<size_t N, typename TArgs>
                struct choose_constructor<N, TArgs> {
                    
                    using type = std::enable_if_t<is_constructible_from_tuple_types_v<T, TArgs>, _constructor<N>>;

                };

        };

    }// namespace utils

}// namespace ssh
