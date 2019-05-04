#include <Arduino.h>
#include <ArduinoSTL.h>

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <sstream>

using namespace std;

void printDouble( double val, unsigned int precision){
// prints val with number of decimal places determine by precision
// NOTE: precision is 1 followed by the number of zeros for the desired number of decimial places
// example: printDouble( 3.1415, 100); // prints 3.14 (two decimal places)

   	Serial.print (int(val));  //prints the int part
   	Serial.print("."); // print the decimal point
   	unsigned int frac;
   	if(val >= 0)
		frac = (val - int(val)) * precision;
   	else
		frac = (int(val)- val ) * precision;
   	int frac1 = frac;
   	while( frac1 /= 10 )
		precision /= 10;
	precision /= 10;
   	while(  precision /= 10)
		Serial.print("0");

   	Serial.println(frac,DEC) ;
}

class TrainingData
{
public:
	TrainingData() = default;
	bool isEof(void)
	{
		return output_counter == output_len;
	}
	void getTopology(vector<unsigned> &topology);

	// Returns the number of input values read from the file:
	unsigned getNextInputs(vector<double> &inputVals);
	unsigned getTargetOutputs(vector<double> &targetOutputVals);
	const unsigned input_sample_size = 2;
	const unsigned output_sample_size = 1;

private:
	const unsigned topology_arr[3] = { 2,4,1 };
	const double input_arr[100] = { 1.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,1.0 };
	const double output_arr[49] = { 1.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0 };
	const unsigned input_len = sizeof(input_arr)/sizeof(double);
	const unsigned output_len = sizeof(output_arr)/sizeof(double);
	const unsigned topology_len = sizeof(topology_arr)/sizeof(unsigned);
	unsigned input_counter = 0;
	unsigned output_counter = 0;
};

void TrainingData::getTopology(vector<unsigned> &topology)
{
	// string line;
	// string label;

	// getline(m_trainingDataFile, line);
	// stringstream ss(line);
	// ss >> label;
	// if(this->isEof() || label.compare("topology:") != 0)
	// {
	// 	abort();
	// }

	// while(!ss.eof())
	// {
	// 	unsigned n;
	// 	ss >> n;
	// 	topology.push_back(n);
	// }
	// return;
	for (unsigned n: TrainingData::topology_arr) {
		Serial.print("topology: ");
		Serial.print(n);
		Serial.println();
		delay(200);
		topology.push_back(n);
	}
	return;
}

// TrainingData::TrainingData(const char* data)
// {
// 	m_trainingDataFile(data);
// }


unsigned TrainingData::getNextInputs(vector<double> &inputVals)
{
    inputVals.clear();

    // string line;
    // getline(m_trainingDataFile, line);
    // stringstream ss(line);

    // string label;
    // ss >> label;
    // if (label.compare("in:") == 0) {
    //     double oneValue;
    //     while (ss >> oneValue) {
    //         inputVals.push_back(oneValue);
    //     }
    // }

	if ((input_counter + input_sample_size) <= input_len) {
		for (unsigned i=0; i != input_sample_size; ++i) {
			Serial.print("input data: ");
			Serial.print(input_arr[input_counter]);
			Serial.println();
			delay(200);
			inputVals.push_back(input_arr[input_counter]);
			++input_counter;
		}
	}

    return inputVals.size();
}

unsigned TrainingData::getTargetOutputs(vector<double> &targetOutputVals)
{
    targetOutputVals.clear();

    // string line;
    // getline(m_trainingDataFile, line);
    // stringstream ss(line);

    // string label;
    // ss>> label;
    // if (label.compare("out:") == 0) {
    //     double oneValue;
    //     while (ss >> oneValue) {
    //         targetOutputVals.push_back(oneValue);
    //     }
    // }

	if ((output_counter + output_sample_size) <= output_len) {
		for (unsigned i=0; i != output_sample_size; ++i) {
			Serial.print("output data: ");
			Serial.print(output_arr[output_counter]);
			Serial.println();
			delay(200);
			targetOutputVals.push_back(output_arr[output_counter]);
			++output_counter;
		}
	}

    return targetOutputVals.size();
}

struct Connection
{
	double weight;
	double deltaWeight;
};

class Neuron;

typedef vector<Neuron> Layer;

// ****************** class Neuron ******************

class Neuron
{
public:
	Neuron(unsigned numOutputs, unsigned myIndex);
	void setOutputVal(double val) { m_outputVal = val; }
	double getOutputVal(void) const { return m_outputVal; }
	void feedForward(const Layer &prevLayer);
	void calcOutputGradients(double targetVals);
	void calcHiddenGradients(const Layer &nextLayer);
	void updateInputWeights(Layer &prevLayer);
private:
	static double eta; // [0.0...1.0] overall net training rate
	static double alpha; // [0.0...n] multiplier of last weight change [momentum]
	static double transferFunction(double x);
	static double transferFunctionDerivative(double x);
	// randomWeight: 0 - 1
	static double randomWeight(void) { return rand() / double(RAND_MAX); }
	double sumDOW(const Layer &nextLayer) const;
	double m_outputVal;
	vector<Connection> m_outputWeights;
	unsigned m_myIndex;
	double m_gradient;
};

double Neuron::eta = 0.15; // overall net learning rate
double Neuron::alpha = 0.5; // momentum, multiplier of last deltaWeight, [0.0..n]


void Neuron::updateInputWeights(Layer &prevLayer)
{
	// The weights to be updated are in the Connection container
	// in the nuerons in the preceding layer

	for(unsigned n = 0; n < prevLayer.size(); ++n)
	{
		Neuron &neuron = prevLayer[n];
		double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;

		double newDeltaWeight = 
				// Individual input, magnified by the gradient and train rate:
				eta
				* neuron.getOutputVal()
				* m_gradient
				// Also add momentum = a fraction of the previous delta weight
				+ alpha
				* oldDeltaWeight;
		neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
		neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
	}
}
double Neuron::sumDOW(const Layer &nextLayer) const
{
	double sum = 0.0;

	// Sum our contributions of the errors at the nodes we feed

	for (unsigned n = 0; n < nextLayer.size() - 1; ++n)
	{
		sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
	}

	return sum;
}

void Neuron::calcHiddenGradients(const Layer &nextLayer)
{
	double dow = sumDOW(nextLayer);
	m_gradient = dow * Neuron::transferFunctionDerivative(m_outputVal);
}
void Neuron::calcOutputGradients(double targetVals)
{
	double delta = targetVals - m_outputVal;
	m_gradient = delta * Neuron::transferFunctionDerivative(m_outputVal);
}

double Neuron::transferFunction(double x)
{
	// tanh - output range [-1.0..1.0]
	return tanh(x);
}

double Neuron::transferFunctionDerivative(double x)
{
	// tanh derivative
	return 1.0 - x * x;
}

void Neuron::feedForward(const Layer &prevLayer)
{
	double sum = 0.0;

	// Sum the previous layer's outputs (which are our inputs)
    // Include the bias node from the previous layer.

	for(unsigned n = 0 ; n < prevLayer.size(); ++n)
	{
		sum += prevLayer[n].getOutputVal() * 
				 prevLayer[n].m_outputWeights[m_myIndex].weight;
	}

	m_outputVal = Neuron::transferFunction(sum);
}

Neuron::Neuron(unsigned numOutputs, unsigned myIndex)
{
	Serial.print("numOutputs:");
	Serial.print(numOutputs);
	Serial.println();
	delay(200);
	for(unsigned c = 0; c < numOutputs; ++c){
		Serial.print("push_back connection");
		Serial.println();
		delay(200);
		m_outputWeights.push_back(Connection());
		Serial.print("save random weight");
		Serial.println();
		delay(200);

		double randnum = randomWeight();
		Serial.print("randnum: ");
		printDouble(randnum, 3);
		Serial.println();
		delay(200);
		m_outputWeights.back().weight = randnum;
		Serial.print("finish loop ");
		Serial.print(c);
		Serial.println();
		delay(200);
	}

	m_myIndex = myIndex;
}
// ****************** class Net ******************
class Net
{
public:
	Net(const vector<unsigned> &topology);
	void feedForward(const vector<double> &inputVals);
	void backProp(const vector<double> &targetVals);
	void getResults(vector<double> &resultVals) const;
	double getRecentAverageError(void) const { return m_recentAverageError; }

private:
	vector<Layer> m_layers; //m_layers[layerNum][neuronNum]
	double m_error;
	double m_recentAverageError;
	static double m_recentAverageSmoothingFactor;
};

double Net::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over

void Net::getResults(vector<double> &resultVals) const
{
	resultVals.clear();

	for(unsigned n = 0; n < m_layers.back().size() - 1; ++n)
	{
		resultVals.push_back(m_layers.back()[n].getOutputVal());
	}
}

void Net::backProp(const std::vector<double> &targetVals)
{
	// Calculate overal net error (RMS of output neuron errors)

	Layer &outputLayer = m_layers.back();
	m_error = 0.0;

	for(unsigned n = 0; n < outputLayer.size() - 1; ++n)
	{
		double delta = targetVals[n] - outputLayer[n].getOutputVal();
		m_error += delta *delta;
	}
	m_error /= outputLayer.size() - 1; // get average error squared
	m_error = sqrt(m_error); // RMS

	// Implement a recent average measurement:

	m_recentAverageError = 
			(m_recentAverageError * m_recentAverageSmoothingFactor + m_error)
			/ (m_recentAverageSmoothingFactor + 1.0);
	// Calculate output layer gradients

	for(unsigned n = 0; n < outputLayer.size() - 1; ++n)
	{
		outputLayer[n].calcOutputGradients(targetVals[n]);
	}
	// Calculate gradients on hidden layers

	for(unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum)
	{
		Layer &hiddenLayer = m_layers[layerNum];
		Layer &nextLayer = m_layers[layerNum + 1];

		for(unsigned n = 0; n < hiddenLayer.size(); ++n)
		{
			hiddenLayer[n].calcHiddenGradients(nextLayer);
		}
	}

	// For all layers from outputs to first hidden layer,
	// update connection weights

	for(unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum)
	{
		Layer &layer = m_layers[layerNum];
		Layer &prevLayer = m_layers[layerNum - 1];

		for(unsigned n = 0; n < layer.size() - 1; ++n)
		{
			layer[n].updateInputWeights(prevLayer);
		}
	}
}

void Net::feedForward(const vector<double> &inputVals)
{
	// Check the num of inputVals euqal to neuronnum expect bias
	assert(inputVals.size() == m_layers[0].size() - 1);

	// Assign {latch} the input values into the input neurons
	for(unsigned i = 0; i < inputVals.size(); ++i){
		m_layers[0][i].setOutputVal(inputVals[i]); 
	}

	// Forward propagate
	for(unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum){
		Layer &prevLayer = m_layers[layerNum - 1];
		for(unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n){
			m_layers[layerNum][n].feedForward(prevLayer);
		}
	}
}
Net::Net(const vector<unsigned> &topology)
{
	unsigned numLayers = topology.size();
	delay(200);
	for(unsigned layerNum = 0; layerNum < numLayers; ++layerNum){
		m_layers.push_back(Layer());
		// numOutputs of layer[i] is the numInputs of layer[i+1]
		// numOutputs of last layer is 0
		Serial.print("numLayers: ");
		Serial.print(topology.size());
		Serial.println();
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 :topology[layerNum + 1];
		Serial.print("numOutputs: ");
		Serial.print(numOutputs);
		Serial.println();
		delay(200);

		// We have made a new Layer, now fill it ith neurons, and
		// add a bias neuron to the layer:
		Serial.print("layerNum: ");
		Serial.print(layerNum);
		Serial.println();
		delay(200);
		Serial.print("topology[layerNum]: ");
		Serial.print(topology[layerNum]);
		Serial.println();
		delay(200);
		for(unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum){
			Serial.print("Start pushback");
			Serial.println();
			delay(200);
			m_layers.back().push_back(Neuron(numOutputs, neuronNum));
			// cout << "Mad a Neuron!" << endl;
			Serial.print("Made a Neuron!");
			Serial.println();
			delay(200);
		}

		// Force the bias node's output value to 1.0. It's the last neuron created above
		Serial.print("Force the bias node's output value to 1.0.");
		Serial.println();
		delay(200);
		m_layers.back().back().setOutputVal(1.0);
	}
}

void showVectorVals(string label, vector<double> &v, unsigned size)
{
	Serial.print("label");
	Serial.println();
	delay(200);
	// cout << label << " ";
	// for(unsigned i = 0; i < v.size(); ++i)
	// {
	// 	// cout << v[i] << " ";
	// 	Serial.print(v[i]);
	// 	Serial.println();
	// 	delay(200);
	// }
	// cout << endl;
	printDouble(v[0], 3);
	for(unsigned i = 0; i < size; ++i)
	{
		// cout << v[i] << " ";
		printDouble(v[i], 3);
		Serial.println();
		delay(200);
	}
	Serial.println();
	delay(200);
}

void setup() {
	Serial.begin(9600);
}

void loop()
{
	Serial.print("Start");
	Serial.println();
	delay(200);
	TrainingData trainData = TrainingData();
	//e.g., {3, 2, 1 }
	vector<unsigned> topology;
	//topology.push_back(3);
	//topology.push_back(2);
	//topology.push_back(1);

	trainData.getTopology(topology);
	Net myNet(topology);

	vector<double> inputVals, targetVals, resultVals;
	int trainingPass = 0;
	Serial.print("Start While");
	Serial.println();
	delay(200);
	while(!trainData.isEof())
	{
		++trainingPass;
		Serial.print("Pass");
		Serial.print(trainingPass);
		Serial.println();
		delay(200);
		// cout << endl << "Pass" << trainingPass;

		// Get new input data and feed it forward:
		if(trainData.getNextInputs(inputVals) != topology[0])
			break;
		showVectorVals(": Inputs :", inputVals, trainData.input_sample_size);
		myNet.feedForward(inputVals);

		// Collect the net's actual results:
		myNet.getResults(resultVals);
		showVectorVals("Outputs:", resultVals, trainData.output_sample_size);

		// Train the net what the outputs should have been:
		trainData.getTargetOutputs(targetVals);
		showVectorVals("Targets:", targetVals, trainData.output_sample_size);
		assert(targetVals.size() == topology.back());

		myNet.backProp(targetVals);

		// Report how well the training is working, average over recnet
		Serial.print("Net recent average error: ");
		printDouble(myNet.getRecentAverageError(), 3);
		Serial.println();
		delay(200);
		// cout << "Net recent average error: "
		//      << myNet.getRecentAverageError() << endl;
	}

	// cout << endl << "Done" << endl;
	Serial.print("DONE");
	Serial.println();
	delay(200);

}
