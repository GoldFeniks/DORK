#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace ssh {

    namespace utils {

        /**
         * Initializes vector values at indexes in range [start, last) by calling func(index)
         */
        template<typename Vector, typename Function, typename Index = size_t>
        void init_vector(Vector& vector, Index first, Index last, Function function) {
            for (Index i = first; i < last; ++i)
                vector[i] = function(i);
        }

        /**
         * Initializes vector values at indexes in range [0, last) by calling func(index)
         */
        template<typename Vector, typename Function, typename Index = size_t>
        void init_vector(Vector& vector, Index last, Function function) {
            init_vector<Vector, Function, Index>(vector, Index(0), last, function);
        }

        /**
         * Initializes vector values at indexes in range [0, vector.size()) by calling func(index)
         */
        template<typename Vector, typename Function>
        void init_vector(Vector& vector, Function function) {
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
            return create_mesh_step(a, (b - a) / (n - 1), n);
        }

    }// namespace utils

}// namespace ssh