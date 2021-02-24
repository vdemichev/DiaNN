#ifndef ACTIVATION_TANH_H_
#define ACTIVATION_TANH_H_

#include <Eigen/Core>
#include "../Config.h"

namespace MiniDNN
{


///
/// \ingroup Activations
///
/// The tanh activation function
///
class Tanh
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    public:
        // a = activation(z) = tanh(z)
        // Z = [z1, ..., zn], A = [a1, ..., an], n observations
        static inline void activate(const Matrix& Z, Matrix& A)
        {
            A.array() = Z.array().tanh();
        }

        // Apply the Jacobian matrix J to a vector f
        // tanh'(x) = 1 - tanh(x)^2
        // J = d_a / d_z = diag(1 - a^2)
        // g = J * f = (1 - a^2) .* f
        // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
        // Note: When entering this function, Z and G may point to the same matrix
        static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                          const Matrix& F, Matrix& G)
        {
            G.array() = (Scalar(1) - A.array().square()) * F.array();
        }

        static std::string return_type()
        {
            return "Tanh";
        }
};


} // namespace MiniDNN


#endif /* ACTIVATION_SIGMOID_H_ */
