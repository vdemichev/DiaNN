#ifndef ACTIVATION_IDENTITY_H_
#define ACTIVATION_IDENTITY_H_

#include <Eigen/Core>
#include "../Config.h"

namespace MiniDNN
{


///
/// \defgroup Activations Activation Functions
///

///
/// \ingroup Activations
///
/// The identity activation function
///
class Identity
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    public:
        // a = activation(z) = z
        // Z = [z1, ..., zn], A = [a1, ..., an], n observations
        static inline void activate(const Matrix& Z, Matrix& A)
        {
            A.noalias() = Z;
        }

        // Apply the Jacobian matrix J to a vector f
        // J = d_a / d_z = I
        // g = J * f = f
        // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
        // Note: When entering this function, Z and G may point to the same matrix
        static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                          const Matrix& F, Matrix& G)
        {
            G.noalias() = F;
        }

        static std::string return_type()
        {
            return "Identity";
        }
};


} // namespace MiniDNN


#endif /* ACTIVATION_IDENTITY_H_ */
