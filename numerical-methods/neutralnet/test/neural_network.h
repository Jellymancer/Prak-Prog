#ifndef NEURAL_NETWORK_H_
#define NEURAL_NETWORK_H_ 

typedef struct {
		int size; 
		double (*function)(double); 
		gsl_vector* data;
		} neural_network;

neural_network* neural_network_alloc(int number_of_hidden_neurons, double(*activation_function)(double));
void neural_network_free(neural_network* NN);
double neural_network_feed_forward(neural_network* network, double x);
void neural_network_trainer(neural_network* NN, gsl_vector* xlist, gsl_vector* ylist);

#endif
