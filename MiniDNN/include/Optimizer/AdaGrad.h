#ifndef OPTIMIZER_ADAGRAD_H_
#define OPTIMIZER_ADAGRAD_H_

#include <Eigen/Core>
#include <map>
#include "../Config.h"
#include "../Optimizer.h"

namespace MiniDNN
{


///
/// \ingroup Optimizers
///
/// The AdaGrad algorithm
///
class AdaGrad: public Optimizer
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;
        typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
        typedef Vector::AlignedMapType AlignedMapVec;

        std::map<const Scalar*, Array> m_history;

    public:
        Scalar m_lrate;
        Scalar m_eps;

        AdaGrad(const Scalar& lrate = Scalar(0.001), const Scalar& eps = Scalar(1e-6)) :
            m_lrate(lrate), m_eps(eps)
        {}

        void reset()
        {
            m_history.clear();
        }

        void update(ConstAlignedMapVec& dvec, AlignedMapVec& vec)
        {
            // Get the accumulated squared gradient associated with this gradient
            Array& grad_square = m_history[dvec.data()];

            // If length is zero, initialize it
            if (grad_square.size() == 0)
            {
                grad_square.resize(dvec.size());
                grad_square.setZero();
            }

            // Update accumulated squared gradient
            grad_square += dvec.array().square();
            // Update parameters
            vec.array() -= m_lrate * dvec.array() / (grad_square.sqrt() + m_eps);
        }
};


} // namespace MiniDNN


#endif /* OPTIMIZER_ADAGRAD_H_ */
