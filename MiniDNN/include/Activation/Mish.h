#ifndef ACTIVATION_MISH_H_
#define ACTIVATION_MISH_H_

#include <Eigen/Core>
#include "../Config.h"

namespace MiniDNN
{


///
/// \ingroup Activations
///
/// The Mish activation function
///
/// From: https://arxiv.org/abs/1908.08681
///
class Mish
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    public:
        // Mish(x) = x * tanh(softplus(x))
        // softplus(x) = log(1 + exp(x))
        // a = activation(z) = Mish(z)
        // Z = [z1, ..., zn], A = [a1, ..., an], n observations
        static inline void activate(const Matrix& Z, Matrix& A)
        {
            // h(x) = tanh(softplus(x)) = (1 + exp(x))^2 - 1
            //                            ------------------
            //                            (1 + exp(x))^2 + 1
            // Let s = exp(-abs(x)), t = 1 + s
            // If x >= 0, then h(x) = (t^2 - s^2) / (t^2 + s^2)
            // If x <= 0, then h(x) = (t^2 - 1) / (t^2 + 1)
            Matrix S = (-Z.array().abs()).exp();
            A.array() = (S.array() + Scalar(1)).square();  // t^2
            S.noalias() = (Z.array() >= Scalar(0)).select(S.cwiseAbs2(), Scalar(1));  // s^2 or 1
            A.array() = (A.array() - S.array()) / (A.array() + S.array());
            A.array() *= Z.array();
        }

        // Apply the Jacobian matrix J to a vector f
        // J = d_a / d_z = diag(Mish'(z))
        // g = J * f = Mish'(z) .* f
        // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
        // Note: When entering this function, Z and G may point to the same matrix
        static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                          const Matrix& F, Matrix& G)
        {
            // Let h(x) = tanh(softplus(x))
            // Mish'(x) = h(x) + x * h'(x)
            // h'(x) = tanh'(softplus(x)) * softplus'(x)
            //       = [1 - h(x)^2] * exp(x) / (1 + exp(x))
            //       = [1 - h(x)^2] / (1 + exp(-x))
            // Mish'(x) = h(x) + [x - Mish(x) * h(x)] / (1 + exp(-x))
            // A = Mish(Z) = Z .* h(Z) => h(Z) = A ./ Z, h(0) = 0.6
            G.noalias() = (Z.array() == Scalar(0)).select(Scalar(0.6), A.cwiseQuotient(Z));
            G.array() += (Z.array() - A.array() * G.array()) / (Scalar(1) + (-Z).array().exp());
            G.array() *= F.array();
        }

        static std::string return_type()
        {
            return "Mish";
        }
};


} // namespace MiniDNN


#endif /* ACTIVATION_MISH_H_ */
