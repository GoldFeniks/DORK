#pragma once

#include <cstddef>
#include <functional>
#include <vector>
#include <type_traits>
#include <tuple>
#include <experimental/type_traits>

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

        template<typename T1, typename T2, template<typename, typename> typename Assigner, template<typename, typename> typename... Assginers>
        struct choose_assigner {

            using type = std::conditional_t<Assigner<T1, T2>::can_assign, Assigner<T1, T2>, typename choose_assigner<T1, T2, Assginers...>::type>;

        };

        template<typename T1, typename T2, template<typename, typename> typename Assigner>
        struct choose_assigner<T1, T2, Assigner> {

            using type = std::enable_if_t<Assigner<T1, T2>::can_assign, Assigner<T1, T2>>;

        };

        template<typename T1, typename T2, template<typename, typename> typename... Assigners>
        void assign(T1& to, T2&& from) {
            choose_assigner<T1, T2, Assigners...>::type::assign(to, std::move(from));
        }

        template<typename T1, typename T2>
        struct move_assigner {

            static constexpr bool can_assign = std::is_assignable_v<T1&, T2&&>;

            static void assign(T1& to, T2&& from) {
                to = std::move(from);
            }

        };

        template<typename T1, typename T2>
        struct copy_assigner {

            static constexpr bool can_assign = std::is_assignable_v<T1&, const T2&>;

            static void assign(T1& to, const T2& from) {
                to = from;
            }

        };

        template<typename T>
        using index_operator_t = decltype(std::declval<T&>()[1]);

        template<typename T>
        using size_method_t = decltype(std::declval<T&>().size());

        template<typename T, template<typename> typename Method>
        struct has_method {

            static constexpr bool value = std::experimental::is_detected<Method, T>::value;

        };

        template<typename T, template<typename> typename Method>
        inline constexpr bool has_method_v = has_method<T, Method>::value;

        template<typename T>
        using is_indexable = has_method<T, index_operator_t>;

        template<typename T>
        inline constexpr bool is_indexable_v = is_indexable<T>::value;

        template<typename T1, typename T2>
        struct loop_assigner {

            static constexpr bool can_assign = is_indexable_v<T1> && is_indexable_v<T2> && has_method_v<T2, size_method_t>;

            static void assign(T1& to, T2&& from) {
                to = constructor<T1>::construct(
                    arguments(from.size()),
                    no_arguments()
                );
                for (size_t i = 0; i < from.size(); ++i) {
                    assign(to[i], std::move(from[i]));
                }
            }

        };


    }// namespace utils

}// namespace ssh
