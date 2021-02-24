#ifndef CALLBACK_VERBOSECALLBACK_H_
#define CALLBACK_VERBOSECALLBACK_H_

#include <Eigen/Core>
#include <iostream>
#include "../Config.h"
#include "../Callback.h"
#include "../Network.h"

namespace MiniDNN
{


///
/// \ingroup Callbacks
///
/// Callback function that prints the loss function value in each mini-batch training
///
class VerboseCallback: public Callback
{
    public:
        void post_training_batch(const Network* net, const Matrix& x, const Matrix& y)
        {
            const Scalar loss = net->get_output()->loss();
            std::cout << "[Epoch " << m_epoch_id << ", batch " << m_batch_id << "] Loss = "
                      << loss << std::endl;
        }

        void post_training_batch(const Network* net, const Matrix& x,
                                 const IntegerVector& y)
        {
            Scalar loss = net->get_output()->loss();
            std::cout << "[Epoch " << m_epoch_id << ", batch " << m_batch_id << "] Loss = "
                      << loss << std::endl;
        }
};


} // namespace MiniDNN


#endif /* CALLBACK_VERBOSECALLBACK_H_ */
