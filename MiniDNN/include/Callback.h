#ifndef CALLBACK_H_
#define CALLBACK_H_

#include <Eigen/Core>
#include "Config.h"

namespace MiniDNN
{


class Network;

///
/// \defgroup Callbacks Callback Functions
///

///
/// \ingroup Callbacks
///
/// The interface and default implementation of the callback function during
/// model fitting. The purpose of this class is to allow users printing some
/// messages in each epoch or mini-batch training, for example the time spent,
/// the loss function values, etc.
///
/// This default implementation is a silent version of the callback function
/// that basically does nothing. See the VerboseCallback class for a verbose
/// version that prints the loss function value in each mini-batch.
///
class Callback
{
    protected:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::RowVectorXi IntegerVector;

    public:
        // Public members that will be set by the network during the training process
        int m_nbatch;   // Number of total batches
        int m_batch_id; // The index for the current mini-batch (0, 1, ..., m_nbatch-1)
        int m_nepoch;   // Total number of epochs (one run on the whole data set) in the training process
        int m_epoch_id; // The index for the current epoch (0, 1, ..., m_nepoch-1)

        Callback() :
            m_nbatch(0), m_batch_id(0), m_nepoch(0), m_epoch_id(0)
        {}

        virtual ~Callback() {}

        // Before training a mini-batch
        virtual void pre_training_batch(const Network* net, const Matrix& x,
                                        const Matrix& y) {}
        //virtual void pre_training_batch(const Network* net, const Matrix& x, const IntegerVector& y) {}

        // After a mini-batch is trained
        virtual void post_training_batch(const Network* net, const Matrix& x,
                                         const Matrix& y) {}
        //virtual void post_training_batch(const Network* net, const Matrix& x, const IntegerVector& y) {}
};


} // namespace MiniDNN


#endif /* CALLBACK_H_ */
