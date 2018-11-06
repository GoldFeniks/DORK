#pragma once

#include <tuple>

namespace ssh {

    namespace utils {

        template<typename T, typename Arguments, typename = typename Arguments::arguments_tuple_t>
        struct is_constructible_from_arguments;

        template<typename T, typename Arguments, typename... Types>
        struct is_constructible_from_arguments<T, Arguments, std::tuple<Types...>> {

            static constexpr bool value = std::is_constructible_v<T, Types...>;

        };

        template<typename T, typename Arguments>
        inline constexpr bool is_constructible_from_arguments_v = 
            is_constructible_from_arguments<T, Arguments>::value;
        
        template<typename... Types>
        struct arguments {

            using arguments_tuple_t = std::tuple<Types...>;

            const arguments_tuple_t arguments_tuple;

            arguments(const Types&... values) : arguments_tuple(std::tuple<Types...>(values...)) {}
            arguments(const arguments&) = default;
            arguments(arguments&&) = default;
            ~arguments() = default;

        };

        using no_arguments = arguments<>;

        namespace {

            template<typename T, size_t N>
            struct constructor_implementation {

                template<typename... Arguments>
                static auto construct(const std::tuple<Arguments...>& arguments) {
                    return std::apply([](auto&&... arguments) { 
                        return T(std::forward<decltype(arguments)>(arguments)...); 
                    }, std::get<N>(arguments).arguments_tuple);
                }

            };
            
        }// namespace

        template<typename T, size_t N, typename Arguments, typename... Tail>
        struct choose_constructor {

            using type = std::conditional_t<is_constructible_from_arguments_v<T, Arguments>,
                                            constructor_implementation<T, N>, 
                                            typename choose_constructor<T, N + 1, Tail...>::type>;

        };

        template<typename T, size_t N, typename Arguments>
        struct choose_constructor<T, N, Arguments> {
            
            using type = std::conditional_t<is_constructible_from_arguments_v<T, Arguments>, 
                                          constructor_implementation<T, N>, void>;

        };

        template<typename T, typename... Arguments>
        using choose_constructor_t = typename choose_constructor<T, 0, Arguments...>::type;

        template<typename T>
        struct constructor {

            template<typename... Arguments>
            static auto construct_from_tuple(const std::tuple<Arguments...>& arguments) {
                return choose_constructor_t<T, Arguments...>::construct(arguments);
            }

            template<typename... Arguments>
            static auto construct(const Arguments... arguments) {
                return construct_from_tuple<Arguments...>(std::tuple(arguments...));
            }

        };

    }// namespace utils

}// namespace ssh