#ifndef ACTIVATION_SIGMOID_H_
#define ACTIVATION_SIGMOID_H_

#include <Eigen/Core>
#include "../Config.h"

namespace MiniDNN
{


///
/// \ingroup Activations
///
/// The sigmoid activation function
///
class Sigmoid
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    public:
        // a = activation(z) = 1 / (1 + exp(-z))
        // Z = [z1, ..., zn], A = [a1, ..., an], n observations
        static inline void activate(const Matrix& Z, Matrix& A)
        {
            A.array() = Scalar(1) / (Scalar(1) + (-Z.array()).exp());
        }

        // Apply the Jacobian matrix J to a vector f
        // J = d_a / d_z = diag(a .* (1 - a))
        // g = J * f = a .* (1 - a) .* f
        // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
        // Note: When entering this function, Z and G may point to the same matrix
        static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                          const Matrix& F, Matrix& G)
        {
            G.array() = A.array() * (Scalar(1) - A.array()) * F.array();
        }

        static std::string return_type()
        {
            return "Sigmoid";
        }
};


} // namespace MiniDNN


#endif /* ACTIVATION_SIGMOID_H_ */
