#pragma once

#include <type_traits>
#include <experimental/type_traits>
#include "constructor.hpp"

namespace dork {

    namespace utils {

        template<typename T, template<typename> typename Method>
        struct has_method {

            static constexpr bool value = std::experimental::is_detected<Method, T>::value;

        };

        template<typename T, template<typename> typename Method>
        inline constexpr bool has_method_v = has_method<T, Method>::value;

        template<typename T>
        using index_operator_t = decltype(std::declval<T&>()[1]);

        template<typename T>
        using size_method_t = decltype(std::declval<T&>().size());

        template<typename T>
        using is_indexable = has_method<T, index_operator_t>;

        template<typename T>
        inline constexpr bool is_indexable_v = is_indexable<T>::value;

        template<
            typename T1, typename T2, 
            template<typename, typename> typename Assigner, 
            template<typename, typename> typename... Assigners>
        struct choose_assigner {

            using type = std::conditional_t<Assigner<T1, T2>::can_assign, 
                                            Assigner<T1, T2>, 
                                            typename choose_assigner<T1, T2, Assigners...>::type>;

        };

        template<typename T1, typename T2, template<typename, typename> typename Assigner>
        struct choose_assigner<T1, T2, Assigner> {

            using type = std::conditional_t<Assigner<T1, T2>::can_assign, Assigner<T1, T2>, void>;

        };

        template<typename T1, typename T2, template<typename, typename> typename... Assigners>
        using choose_assigner_t = typename choose_assigner<T1, T2, Assigners...>::type;

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

        template<typename T1, typename T2>
        struct loop_assigner {

            static constexpr bool can_assign = is_indexable_v<T1> && 
                                               is_indexable_v<T2> && 
                                               has_method_v<T2, size_method_t>;

            static void assign(T1& to, T2&& from) {
                to = constructor<T1>::construct(
                    arguments(from.size()),
                    no_arguments()
                );
                for (size_t i = 0; i < from.size(); ++i) {
                    to[i] = std::move(from[i]);
                }
            }

        };

        template<typename T1, typename T2, template<typename, typename> typename... Assigners>
        inline void assign(T1& to, T2&& from) {
            choose_assigner_t<T1, T2, Assigners...>::assign(to, std::move(from));
        }

        template<typename T1, typename T2>
        void default_assign(T1& to, T2&& from) {
            utils::assign<T1, T2, utils::move_assigner, utils::copy_assigner, utils::loop_assigner>(to, std::move(from));
        }

    }// namespace utils
    
}// namespace dork
