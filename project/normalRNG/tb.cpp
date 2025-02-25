#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "dut.hpp"

// Kernel function declaration
extern "C" void dut(DT temp_inv[NCHAINS],
                            DT sigma[NCHAINS],
                            DT sample_output[NSAMPLES_MAX][3],
                            unsigned int nSamples,
							DT y[NCHAINS][NSAMPLES_MAX],
							DT data[Observation][2]);

int main() {
    std::cout << "*************" << std::endl;
    std::cout << "MCMC Demo v1.0" << std::endl;
    std::cout << "*************" << std::endl;
    std::cout << std::endl;


    xf::fintech::MT19937 uniformRNG(42);
    int a = 3;
    int b = 20;
    int c = 5;
    DT data[Observation][2];
    for (int i = 0; i < Observation; i++) {
        DT x = xf::fintech::inverseCumulativeNormalAcklam<DT>(uniformRNG.next()) * 30;
        DT e = xf::fintech::inverseCumulativeNormalAcklam<DT>(uniformRNG.next());
        DT y = a * x + b + e;
        data[i][0] = x;
        data[i][1] = y;
    }

    /*for (unsigned int i = 0; i < 2; ++i) {
            	std::cout << "Chain " << " samples:" << std::endl;
            	for (unsigned int j = 0; j < Observation; ++j) {
            		std::cout << data[j][i] << std::endl;
            	}
            	std::cout << std::endl; // Separate chains
            }*/


    unsigned int num_samples = 20000;  // Example number of samples
    //unsigned int num_burn = 2000;      // Example number of burn-in samples

    std::vector<DT> temp_inv(NCHAINS);
    std::vector<DT> temp(NCHAINS);
    std::vector<DT> sigma(NCHAINS);
    //std::vector<DT> sample(num_samples);
    DT sample_output[NSAMPLES_MAX][3];

    DT y[NCHAINS][NSAMPLES_MAX];



    for (unsigned int n = 0; n < NCHAINS; n++) {
        sigma[n] = 0.4;
        temp_inv[0] = 1;
        for (int i = 1; i < NCHAINS; i++) {
        	temp_inv[i] = temp_inv[i - 1] - 1.0 / i;
        }
    }


    // Call the kernel function
    //dut(temp_inv.data(), sigma.data(), sample.data(), num_samples);
    //dut(temp_inv.data(), sigma.data(), sample_output, num_samples, y);
    dut(temp_inv.data(), sigma.data(), sample_output, num_samples, y, data);

    /*for (unsigned int i = 0; i < num_samples; ++i) {
        std::cout << sample[i] << std::endl;
    }*/


    // Output results
    /*for (unsigned int i = 0; i < 3; ++i) {
    	std::cout << "Chain " << i << " samples:" << std::endl;
    	for (unsigned int j = num_burn; j < num_samples; ++j) {
    		std::cout << sample_output[j][i] << std::endl;
    	}
    	std::cout << std::endl; // Separate chains
    }
    */



    /*for (unsigned int chain = 0; chain < NCHAINS; ++chain) {
        	std::cout << "Chain " << chain << " samples:" << std::endl;
        	for (unsigned int i = 0; i < num_samples; ++i) {
        		std::cout << y[chain][i] << std::endl;
        	}
        	std::cout << std::endl; // Separate chains
        }
        */

    std::cout << "######## theta1 #######" << std::endl;
    for (int i = 0; i < num_samples; i++) {
        std::cout << sample_output[i][0] << std::endl;
    }

    std::cout << "######## theta2 #######" << std::endl;
    for (int i = 0; i < num_samples; i++) {
        std::cout << sample_output[i][1] << std::endl;
    }

    std::cout << "######## theta3 #######" << std::endl;
    for (int i = 0; i < num_samples; i++) {
        std::cout << sample_output[i][2] << std::endl;
    }



    return 0;
}

