#include "std_includes.h"
#include "matrix.h"
#include "function.h"
#include "layer.h"
#include "network.h"

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

// types of loss functions available to optimize under
typedef enum LOSS_FUNCTION_ {
    CROSS_ENTROPY_LOSS, 
    MEAN_SQUARED_ERROR
} LOSS_FUNCTION;

// convenience struct for easier parameter filling
class NN {
public:
    Network* network = NULL;
    DataSet* data = NULL;
    DataSet* classes = NULL;
    LOSS_FUNCTION lossFunction = CROSS_ENTROPY_LOSS;
    size_t batchSize = 1;
    float learningRate = 0.003;
    float searchTime = 0;
    float regularizationStrength = 0.000000001;
    float B1 = 0.9, B2 = 0.999;
    int maxIters = 1;
    int shuffle = 1;
    int verbose = 0;

	int id; // for bagging
	std::minstd_rand random;

	NN() {}

	void seed(int _id) {
		id = _id;
		random.seed(id + 1);
	}

	void optimise() {
		assert(network->layers[0]->size == data->cols);
		assert(data->rows == classes->rows);
		assert(network->layers[network->numLayers - 1]->size == classes->cols);
		assert(batchSize <= data->rows);
		assert(maxIters >= 1);

		int i, j, k;

		// these will be reused per training instance
		Matrix** errori = (Matrix**)alloca(network->numLayers * sizeof(Matrix*));
		Matrix** dWi = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** dbi = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** regi = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix* beforeOutputT = createMatrixZeroes(network->layers[network->numLayers - 2]->size, 1);
		for (i = 0; i < network->numConnections; i++) {
			errori[i] = createMatrixZeroes(1, network->layers[i]->size);
			dWi[i] = createMatrixZeroes(network->connections[i]->weights->rows, network->connections[i]->weights->cols);
			dbi[i] = createMatrixZeroes(1, network->connections[i]->bias->cols);
			regi[i] = createMatrixZeroes(network->connections[i]->weights->rows, network->connections[i]->weights->cols);
		}
		errori[i] = createMatrixZeroes(1, network->layers[i]->size);

		// these will be reused per training instance if network has hidden layers
		int numHidden = network->numLayers - 2;
		Matrix** WTi = NULL, ** errorLastTi = NULL, ** fprimei = NULL, ** inputTi = NULL;
		if (numHidden > 0) {
			WTi = (Matrix**)malloc(sizeof(Matrix*) * numHidden);
			errorLastTi = (Matrix**)malloc(sizeof(Matrix*) * numHidden);
			fprimei = (Matrix**)malloc(sizeof(Matrix*) * numHidden);
			inputTi = (Matrix**)malloc(sizeof(Matrix*) * numHidden);
			for (k = 0; k < numHidden; k++) {
				WTi[k] = createMatrixZeroes(network->connections[k + 1]->weights->cols, network->connections[k + 1]->weights->rows);
				errorLastTi[k] = createMatrixZeroes(1, WTi[k]->cols);
				fprimei[k] = createMatrixZeroes(1, network->connections[k]->to->size);
				inputTi[k] = createMatrixZeroes(network->connections[k]->from->size, 1);
			}
		}

		// these will be reused per epoch
		Matrix** dWi_avg = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** dbi_avg = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** Mw = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** Mb = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** Vw = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		Matrix** Vb = (Matrix**)alloca(network->numConnections * sizeof(Matrix*));
		for (i = 0; i < network->numConnections; i++) {
			dWi_avg[i] = createMatrixZeroes(network->connections[i]->weights->rows, network->connections[i]->weights->cols);
			dbi_avg[i] = createMatrixZeroes(1, network->connections[i]->bias->cols);
			Mw[i] = createMatrixZeroes(network->connections[i]->weights->rows, network->connections[i]->weights->cols);
			Mb[i] = createMatrixZeroes(1, network->connections[i]->bias->cols);
			Vw[i] = createMatrixZeroes(network->connections[i]->weights->rows, network->connections[i]->weights->cols);
			Vb[i] = createMatrixZeroes(1, network->connections[i]->bias->cols);
		}

		int numBatches = (data->rows / batchSize) + (data->rows % batchSize != 0 ? 1 : 0);
		int training, batch, epoch, layer;
		DataSet** dataBatches, ** classBatches;
		epoch = 1;
		float B1t = B1, B2t = B2;
		while (epoch <= maxIters) {
			// shuffle all data and classes but maintain training/class alignment
			if (shuffle != 0) {
				shuffleTogether(data, classes, random);
			}

			// split into overall batches
			dataBatches = createBatches(data, numBatches);
			classBatches = createBatches(classes, numBatches);
			for (batch = 0; batch < numBatches && epoch <= maxIters; batch++, epoch++) {
				// find current batch
				int curBatchSize = dataBatches[batch]->rows;
				DataSet* batchTraining = dataBatches[batch];
				DataSet* batchClasses = classBatches[batch];
				Matrix** splitTraining = splitRows(batchTraining);
				Matrix** splitClasses = splitRows(batchClasses);
				for (training = 0; training < curBatchSize; training++) {
					// current data point to train on
					Matrix* example = splitTraining[training];
					Matrix* target = splitClasses[training];

					// pass error forward
					forwardPass(network, example);

					// calculate each iteration of backpropagation
					for (layer = network->numLayers - 1; layer > 0; layer--) {
						Layer* to = network->layers[layer];
						Connection* con = network->connections[layer - 1];
						if (layer == network->numLayers - 1) {
							// calculate output layer's error
							copyValuesInto(to->input, errori[layer]);
							if (lossFunction == CROSS_ENTROPY_LOSS) {
								for (j = 0; j < errori[layer]->cols; j++) {
									errori[layer]->data[j] -= target->data[j];
								}
							}
							else {
								for (j = 0; j < errori[layer]->cols; j++) {
									errori[layer]->data[j] -= target->data[j];
								}
							}

							// calculate dWi and dbi
							transposeInto(con->from->input, beforeOutputT);
							multiplyInto(beforeOutputT, errori[layer], dWi[layer - 1]);
							copyValuesInto(errori[layer], dbi[layer - 1]);
						}
						else {
							// calculate error term for hidden layer
							int hiddenLayer = layer - 1;
							transposeInto(network->connections[layer]->weights, WTi[hiddenLayer]);
							multiplyInto(errori[layer + 1], WTi[hiddenLayer], errorLastTi[hiddenLayer]);
							copyValuesInto(con->to->input, fprimei[hiddenLayer]);
							float(*derivative)(float) = activationDerivative(con->to->activation);
							for (j = 0; j < fprimei[hiddenLayer]->cols; j++) {
								fprimei[hiddenLayer]->data[j] = derivative(fprimei[hiddenLayer]->data[j]);
							}
							hadamardInto(errorLastTi[hiddenLayer], fprimei[hiddenLayer], errori[layer]);

							// calculate dWi and dbi
							transposeInto(con->from->input, inputTi[hiddenLayer]);
							multiplyInto(inputTi[hiddenLayer], errori[layer], dWi[layer - 1]);
							copyValuesInto(errori[layer], dbi[layer - 1]);
						}
					}

					// add one example's contribution to total gradient
					for (i = 0; i < network->numConnections; i++) {
						addTo(dWi[i], dWi_avg[i]);
						addTo(dbi[i], dbi_avg[i]);
					}

					// zero out reusable matrices
					zeroMatrix(beforeOutputT);
					for (i = 0; i < network->numConnections; i++) {
						zeroMatrix(errori[i]);
						zeroMatrix(dWi[i]);
						zeroMatrix(dbi[i]);
					}
					zeroMatrix(errori[i]);
					if (numHidden > 0) {
						for (i = 0; i < numHidden; i++) {
							zeroMatrix(WTi[i]);
							zeroMatrix(errorLastTi[i]);
							zeroMatrix(fprimei[i]);
							zeroMatrix(inputTi[i]);
						}
					}
				}

				// add regularization
				for (i = 0; i < network->numConnections; i++) {
					copyValuesInto(network->connections[i]->weights, regi[i]);
					scalarMultiply(regi[i], regularizationStrength);
					addTo(regi[i], dWi_avg[i]);
				}

				// ADAM
				for (i = 0; i < network->numConnections; i++) {
					matrix_mul_add(Mw[i], B1, dWi_avg[i], 1.0 - B1);
					matrix_mul_add(Mb[i], B1, dbi_avg[i], 1.0 - B1);
					matrix_mul_add(Vw[i], B2, dWi_avg[i], 1.0 - B2, true);
					matrix_mul_add(Vb[i], B2, dbi_avg[i], 1.0 - B2, true);
					float rate = learningRate * (sqrt(1.0 - B2t) / (1.0 - B1t));
					adam_update(dWi_avg[i], Mw[i], Vw[i], rate);
					adam_update(dbi_avg[i], Mb[i], Vb[i], rate);
					B1t *= B1, B2t *= B2;
				}

				// adjust weights and bias
				for (i = 0; i < network->numConnections; i++) {
					scalarMultiply(dWi_avg[i], -1);
					scalarMultiply(dbi_avg[i], -1);
					addTo(dWi_avg[i], network->connections[i]->weights);
					addTo(dbi_avg[i], network->connections[i]->bias);
				}

				// zero out reusable average matrices and regularization matrices
				for (i = 0; i < network->numConnections; i++) {
					zeroMatrix(dWi_avg[i]);
					zeroMatrix(dbi_avg[i]);
					zeroMatrix(regi[i]);
				}

				// free list of training examples
				for (i = 0; i < curBatchSize; i++) {
					free(splitTraining[i]);
					free(splitClasses[i]);
				}
				free(splitTraining);
				free(splitClasses);

				// if verbose is set, print loss every 100 epochs
				if (verbose != 0) {
					if (epoch % ((data->rows * 2) / batchSize) == 0 || epoch == 1) {
						forwardPassDataSet(network, data);
						if (lossFunction == CROSS_ENTROPY_LOSS) {
							printf("Batch %d: loss is %f\n", epoch, crossEntropyLoss(network, getOuput(network), classes, regularizationStrength));
						}
						else {
							printf("Batch %d: loss is %f\n", epoch, meanSquaredError(network, getOuput(network), classes, regularizationStrength));
						}
					}
				}
			}

			// free batch data for next batch creation
			for (i = 0; i < numBatches; i++) {
				free(dataBatches[i]);
				free(classBatches[i]);
			}
			free(dataBatches);
			free(classBatches);
		}

		// free all reusable matrices
		destroyMatrix(beforeOutputT);
		for (i = 0; i < network->numConnections; i++) {
			destroyMatrix(errori[i]);
			destroyMatrix(dWi[i]);
			destroyMatrix(dbi[i]);
		}
		destroyMatrix(errori[i]);

		for (i = 0; i < network->numConnections; i++) {
			destroyMatrix(dWi_avg[i]);
			destroyMatrix(dbi_avg[i]);
			destroyMatrix(Mw[i]);
			destroyMatrix(Mb[i]);
			destroyMatrix(Vw[i]);
			destroyMatrix(Vb[i]);
			destroyMatrix(regi[i]);
		}

		if (numHidden > 0) {
			for (i = 0; i < numHidden; i++) {
				destroyMatrix(WTi[i]);
				destroyMatrix(errorLastTi[i]);
				destroyMatrix(fprimei[i]);
				destroyMatrix(inputTi[i]);
			}
			free(WTi);
			free(errorLastTi);
			free(fprimei);
			free(inputTi);
		}
	}
};

#endif