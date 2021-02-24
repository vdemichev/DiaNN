#ifndef UTILS_FACTORY_H_
#define UTILS_FACTORY_H_

#include <string>
#include <map>
#include <stdexcept>
#include "../Config.h"
#include "IO.h"
#include "Enum.h"

#include "../Layer.h"
#include "../Layer/FullyConnected.h"
#include "../Layer/Convolutional.h"
#include "../Layer/MaxPooling.h"

#include "../Activation/Identity.h"
#include "../Activation/ReLU.h"
#include "../Activation/Sigmoid.h"
#include "../Activation/Softmax.h"
#include "../Activation/Tanh.h"
#include "../Activation/Mish.h"

#include "../Output.h"
#include "../Output/RegressionMSE.h"
#include "../Output/BinaryClassEntropy.h"
#include "../Output/MultiClassEntropy.h"

namespace MiniDNN
{

namespace internal
{


// Create a layer from the network meta information and the index of the layer
inline Layer* create_layer(const std::map<std::string, int>& map, int index)
{
    std::string ind = internal::to_string(index);
    const int lay_id = map.find("Layer" + ind)->second;
    const int act_id = map.find("Activation" + ind)->second;
    Layer* layer;

    if (lay_id == FULLY_CONNECTED)
    {
        const int in_size = map.find("in_size" + ind)->second;
        const int out_size = map.find("out_size" + ind)->second;

        switch (act_id)
        {
        case IDENTITY:
            layer = new FullyConnected<Identity>(in_size, out_size);
            break;
        case RELU:
            layer = new FullyConnected<ReLU>(in_size, out_size);
            break;
        case SIGMOID:
            layer = new FullyConnected<Sigmoid>(in_size, out_size);
            break;
        case SOFTMAX:
            layer = new FullyConnected<Softmax>(in_size, out_size);
            break;
        case TANH:
            layer = new FullyConnected<Tanh>(in_size, out_size);
            break;
        case MISH:
            layer = new FullyConnected<Mish>(in_size, out_size);
            break;
        default:
            throw std::invalid_argument("[function create_layer]: Activation is not of a known type");
        }

    } else if (lay_id == CONVOLUTIONAL) {
        const int in_width = map.find("in_width" + ind)->second;
        const int in_height = map.find("in_height" + ind)->second;
        const int in_channels = map.find("in_channels" + ind)->second;
        const int out_channels = map.find("out_channels" + ind)->second;
        const int window_width = map.find("window_width" + ind)->second;
        const int window_height = map.find("window_height" + ind)->second;

        switch(act_id)
        {
        case IDENTITY:
            layer = new Convolutional<Identity>(in_width, in_height, in_channels,
                                                out_channels, window_width, window_height);
            break;
        case RELU:
            layer = new Convolutional<ReLU>(in_width, in_height, in_channels,
                                            out_channels, window_width, window_height);
            break;
        case SIGMOID:
            layer = new Convolutional<Sigmoid>(in_width, in_height, in_channels,
                                               out_channels, window_width, window_height);
            break;
        case SOFTMAX:
            layer = new Convolutional<Softmax>(in_width, in_height, in_channels,
                                               out_channels, window_width, window_height);
            break;
        case TANH:
            layer = new Convolutional<Tanh>(in_width, in_height, in_channels,
                                            out_channels, window_width, window_height);
            break;
        case MISH:
            layer = new Convolutional<Mish>(in_width, in_height, in_channels,
                                            out_channels, window_width, window_height);
            break;
        default:
            throw std::invalid_argument("[function create_layer]: Activation is not of a known type");
        }

    } else if (lay_id == MAX_POOLING) {
        const int in_width = map.find("in_width" + ind)->second;
        const int in_height = map.find("in_height" + ind)->second;
        const int in_channels = map.find("in_channels" + ind)->second;
        const int pooling_width = map.find("pooling_width" + ind)->second;
        const int pooling_height = map.find("pooling_height" + ind)->second;

        switch (act_id)
        {
        case IDENTITY:
            layer = new MaxPooling<Identity>(in_width, in_height, in_channels,
                                             pooling_width, pooling_height);
            break;
        case RELU:
            layer = new MaxPooling<ReLU>(in_width, in_height, in_channels,
                                         pooling_width, pooling_height);
            break;
        case SIGMOID:
            layer = new MaxPooling<Sigmoid>(in_width, in_height, in_channels,
                                            pooling_width, pooling_height);
            break;
        case SOFTMAX:
            layer = new MaxPooling<Softmax>(in_width, in_height, in_channels,
                                            pooling_width, pooling_height);
            break;
        case TANH:
            layer = new MaxPooling<Tanh>(in_width, in_height, in_channels,
                                         pooling_width, pooling_height);
            break;
        case MISH:
            layer = new MaxPooling<Mish>(in_width, in_height, in_channels,
                                         pooling_width, pooling_height);
            break;
        default:
            throw std::invalid_argument("[function create_layer]: Activation is not of a known type");
        }

    } else {

        throw std::invalid_argument("[function create_layer]: Layer is not of a known type");
    }

    layer->init();
    return layer;
}

// Create an output layer from the network meta information
inline Output* create_output(const std::map<std::string, int>& map)
{
    Output* output;
    int out_id = map.find("OutputLayer")->second;

    switch (out_id)
    {
    case REGRESSION_MSE:
        return new RegressionMSE();
    case BINARY_CLASS_ENTROPY:
        return new BinaryClassEntropy();
    case MULTI_CLASS_ENTROPY:
        return new MultiClassEntropy();
    default:
        throw std::invalid_argument("[function create_output]: Output is not of a known type");
    }

    return output;
}


} // namespace internal

} // namespace MiniDNN


#endif /* UTILS_FACTORY_H_ */
