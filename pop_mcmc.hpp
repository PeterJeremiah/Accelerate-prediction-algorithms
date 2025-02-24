/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**********
 * Copyright (c) 2019, Xilinx, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * **********/
/**
 *  @file pop_mcmc.hpp
 *  @brief  Implementation of Population Markov Chain Monte Carlo (MCMC)
 *
 *  $DateTime: 2019/07/24 12:00:00 $
 */

#ifndef _MCMC_CORE_
#define _MCMC_CORE_

#include "xf_fintech/rng.hpp"
#include "hls_math.h"
#include <hls_stream.h>
#include <random>

namespace xf {
namespace fintech {
namespace internal {

//proposal提议
template <typename DT>
void Proposal (DT prec_theta[3], DT search_width, xf::fintech::MT19937& uniformRNG, DT new_theta[3]) {
    for (int i = 0; i < 3; i++) {
        new_theta[i] = prec_theta[i] + search_width * xf::fintech::inverseCumulativeNormalAcklam<DT>(uniformRNG.next());
    }    
}


//似然函数

template <typename DT, unsigned int Observation>
DT Likelihood (DT data[Observation][2], DT theta[3]) {
    DT Likelihood_out = 0;
    for (int i = 0; i < Observation; i++) {
        DT xs = data[i][0];
        DT ys = data[i][1];

        Likelihood_out += -0.5 * std::log(2.0 * M_PI * theta[2] * theta[2]) - (ys - (theta[0] * xs + theta[1])) * (ys - (theta[0] * xs + theta[1])) / (2 * theta[2] * theta[2]);
    }
    return Likelihood_out;
}


//先验概率
template <typename DT>
DT Prior (DT theta[3]) {
    DT Prior_out = -0.5 * std::log((2 * M_PI) * (2 * M_PI) * 100.0 * 100.0 * 100.0) - 0.005 * (theta[0] * theta[0] + theta[1] * theta[1] + theta[2] * theta[2]);
    return Prior_out;

}


template <typename DT>
DT TargetDist(DT x[3], DT temp_inv) {
#pragma HLS inline
    DT result;
    DT exponent = -0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    DT norm_factor = std::sqrt(2 * M_PI);
    result = std::log(std::exp(exponent) / norm_factor) * temp_inv;
    return result;
}

/**
* @brief Calculates final transformation of Gaussian Sample.
*
*@tparam DT         Data type used in whole function (double by default)
*@param[in] in    - Sample from Uniform Distribution \n
*@param[in] mu    - Expected value for Normal Distribution  \n
*@param[in] sigma - Sigma for Proposal generation \n
*@return          - Generated Sample
*/
/*template <typename DT>
DT GaussTransform(DT in, DT mu, DT sigma) {
    DT result = in * sigma + mu;
    return result;
}
*/
/**
* @brief Probability evaluation function. \n
* It Generates samples for all chains. Metropolis sampler is used in this function. \n
* Fully pipelined for chains. \n
* During Probability evaluation gauss sample for next sample is generated in parallel, \n
* this allows to save half of the time for probability evaluation. \n
* Part of the dataflow streaming region.
*
*@tparam DT data type used in whole function (double by default)
*@tparam NCHAINS Number of chains
*@param[in] chain_in    - Previous samples for each chains \n
*@param[in] gauss       - Gaussian sample proposal on [0:1] for current sample (1/Temp) \n
*@param[out] gauss_next - Gaussian sample proposal on [0:1] for next sample \n
*@param[out] chain_out  - Samples streaming output  \n
*@param[in] uniformRNG  - Pointer to Uniform RNG for Accept/Reject \n
*@param[in] temp_inv    - Array of Inverted temperatures of the chain that density is generate for (1/Temp) \n
*@param[in] sigma       - Array of sigmas for Proposal generation for each chain  \n
*/

template <typename DT, unsigned int NCHAINS, unsigned int Observation>
void ProbEval(DT chain_in[NCHAINS][3],
              DT chain_out[NCHAINS][3],
              DT gauss[NCHAINS],
              DT gauss_next[NCHAINS],
              xf::fintech::MT19937& uniformRNG,
              DT temp_inv[NCHAINS],
              DT sigma[NCHAINS],
              DT data[Observation][2],
              DT width) {



#pragma HLS inline off
    DT xStar;
    DT alpha;
    DT u;
    DT sample_buff;
    DT theta[3];
    DT new_theta[3];

PROB_EVALUATION_LOOP:
    for (int n = 0; n < NCHAINS; n++) {
#pragma HLS pipeline

        for (int i = 0; i < 3; i++) {
            theta[i] = chain_in[n][i];
        }
        internal::Proposal<DT>(theta, width, uniformRNG, new_theta);
        DT likelihood_new = internal::Likelihood<DT, Observation>(data, new_theta);
        DT likelihood_old = internal::Likelihood<DT, Observation>(data, theta);
        DT prior_new = internal::Prior<DT>(new_theta);
        DT prior_old = internal::Prior<DT>(theta);
        DT likelihood_prior_proposal_ratio = likelihood_new - likelihood_old + prior_new - prior_old;
        // CALCULATE THE ACCEPTANCE PROBABILITY
        //xStar = GaussTransform<DT>(gauss[n], chain_in[n], sigma[n]);
        
        //alpha = TargetDist<DT>(xStar, temp_inv[n]) - TargetDist<DT>(chain_in[n], temp_inv[n]);
        //DT in = uniformRNG_2.next();
        //gauss_next[n] = xf::fintech::inverseCumulativeNormalAcklam<DT>(in);
        // ACCEPT OR REJECT?
        u = uniformRNG.next();
        if (std::log(u) < likelihood_prior_proposal_ratio) {
            for (int i = 0; i < 3; i++) {
                chain_out[n][i] = new_theta[i];
            }

        } else {
            for (int i = 0; i < 3; i++) {
                chain_out[n][i] = chain_in[n][i];
            }
        }
    } // end of PROB_EVALUATION_LOOP
}


/**
* @brief Chain Exchange function. \n
* Calculates exchange ratio and exchanges chains if needed. \n
* Fully pipelined for chains. \n
* Part of the dataflow streaming region.
*
*@tparam DT data type used in whole function (double by default)
*@tparam NCHAINS Number of chains
*@param[in]  chain_in    - Current sample streaming input interface \n
*@param[in]  chain_out   - Array of generated samples for each chain samples.  \n
*@param[in]  temp_inv    - Array of Inverted temperatures of the chain that density is generate for (1/Temp) \n
*/

template <typename DT, unsigned int NCHAINS>
void ChainExchange(DT chain_in[NCHAINS][3], DT chain_out[NCHAINS][3], DT temp_inv[NCHAINS]) {
    DT chain_buff[NCHAINS][3];

    static bool even;

    bool last_read = 0;
    unsigned int sequence = 0;
    DT u;
    DT alpha_ex;
    static xf::fintech::MT19937 uniformRNG_ex(71);
    // Do first read if even pairs are exchanging because first chain is not reached then.
    if (even) {
        //chain_buff[0] = chain_in.read();
        for (int i = 0; i < 3; i++) {
            chain_buff[0][i] = chain_in[sequence][i];
        }
        sequence++;
    }

EXCHANGE_LOOP:
    for (int n = 1 + even; n < NCHAINS; n = n + 2) { // exchange loop
#pragma HLS LOOP_TRIPCOUNT min = 4 max = 5
#pragma HLS pipeline II = 2
        for (int i = 0; i < 3; i++) {
            chain_buff[n - 1][i] = chain_in[sequence][i];
            chain_buff[n][i] = chain_in[sequence + 1][i];
        }
        sequence += 2;
        // alpha_ex =
        // TargetDist<DT>(chain_buff[n],temp_inv[n-1])*TargetDist<DT>(chain_buff[n-1],temp_inv[n])/(TargetDist<DT>(chain_buff[n],temp_inv[n])*TargetDist<DT>(chain_buff[n-1],temp_inv[n-1]));
        alpha_ex = TargetDist<DT>(chain_buff[n], temp_inv[n - 1]) + TargetDist<DT>(chain_buff[n - 1], temp_inv[n]) -
                   (TargetDist<DT>(chain_buff[n], temp_inv[n]) + TargetDist<DT>(chain_buff[n - 1], temp_inv[n - 1]));

        u = uniformRNG_ex.next();

        /*if (std::log(u) < alpha_ex) {
            chain_out[n - 1] = chain_buff[n];
            chain_out[n] = chain_buff[n - 1];
        } else {
            chain_out[n - 1] = chain_buff[n - 1];
            chain_out[n] = chain_buff[n];
        }*/
        if (std::log(u) < alpha_ex) {
            for (int i = 0; i < 3; i++) {
                chain_out[n - 1][i] = chain_buff[n][i];
                chain_out[n][i] = chain_buff[n - 1][i];
            }
        } else {
            for (int i = 0; i < 3; i++) {
                chain_out[n - 1][i] = chain_buff[n - 1][i];
                chain_out[n][i] = chain_buff[n][i];
            }
        }

        if (NCHAINS - n == 2) {
            last_read = 1;
        }
    } // end of EXCHANGE_LOOP
    // Do last read from stream if last iteration hasn't reached last chain.
    /*if (last_read) {
        chain_out[NCHAINS - 1] = chain_in.read();
    }*/
    if (last_read) {
        for (int i = 0; i < 3; i++) {
            chain_out[NCHAINS - 1][i] = chain_in[sequence][i];
        }
    }
    // Echagning odd or even pairs of chains
    even = !even;
}

/**
* @brief Wraping function for dataflow region. /n
*
*@tparam DT data type used in whole function (double by default)
*@tparam NCHAINS Number of chains
*@param[in] chain       - Previous samples for each chains \n
*@param[in] gauss       - Gaussian sample proposal on [0:1] for current sample (1/Temp) \n
*@param[out] gauss_next - Gaussian sample proposal on [0:1] for next sample \n
*@param[out] chain_out  - Array of generated samples  \n
*@param[in] uniformRNG  - Pointer to Uniform RNG for Accept/Reject \n
*@param[in] temp_inv    - Array of Inverted temperatures of the chain that density is generate for (1/Temp) \n
*@param[in] sigma       - Array of sigmas for Proposal generation for each chain  \n
*/

} // internal


/**
* @brief Top level Kernel function. Consists of INIT_LOOP and main sample loop: SAMPLES_LOOP \n
* \n
* Generates sample from target distribution function.\n
* Uses multiple Markov Chains to allow drawing samples from multi mode target distribution functions. \n
* Proposal is generated ussing Normal Distribution  \n
*@tparam DT             - Data type used in whole function (double by default)
*@tparam NCHAINS        - Number of chains
*@tparam NSAMPLES_MAX   - Maximum Number of chains for synthesis purpose
*@param[in] temp_inv    - Array of Inverted temperatures of the chain that density is generate for (1/Temp) \n
*@param[in] sigma       - Array of sigmas for Proposal generation for each chain  \n
*@param[in] nSamples    - Number of samples to generate  \n
*@param[out] x          - Sample output  \n
*/


template <typename DT, unsigned int NCHAINS, unsigned int NSAMPLES_MAX, unsigned int Observation>
void McmcCore(DT temp_inv[NCHAINS], DT sigma[NCHAINS], DT output[NSAMPLES_MAX][3], unsigned int nSamples, DT y[NCHAINS][NSAMPLES_MAX], DT data[Observation][2]) {
    
    DT width = 0.2;
    //在0-1之间生成一个均匀分布的随机数赋值给theta的初始值
    DT chain_in[NCHAINS][3];
    DT chain_out[NCHAINS][3];
    xf::fintech::MT19937 uniformRNG(42);
    DT theta[3];
    DT u;
    for (int i = 0; i < 3; i++) {
        theta[i] = uniformRNG.next();
        //std::cout << "theta[" << i << "] = " << theta[i] << std::endl;
    }
    for (int i = 0; i < NCHAINS; i++) {
        for (int j = 0; j < 3; j++) {
            chain_in[i][j] = uniformRNG.next();
        }
    }
    unsigned int accepted = 0;
    unsigned int rejected = 0;
    for (int i = 0; i < nSamples; i ++){
        internal::ProbEval<DT, NCHAINS, Observation>(chain_in, chain_out, theta, theta, uniformRNG, temp_inv, sigma, data, width);
        for (int j = 0; j < NCHAINS; j++) {
            for (int k = 0; k < 3; k++) {
                chain_in[j][k] = chain_out[j][k];
            }
        }
        internal::ChainExchange<DT, NCHAINS>(chain_in, chain_out, temp_inv);
        for (int j = 0; j < NCHAINS; j++) {
            if (j == 0) {
                for (int k = 0; k < 3; k++) {
                //output[accepted][k] = chain_out[j][k];
                //std::cout << "output[" << accepted << "][" << k << "] = " << chain_out[j][k] << std::endl;
                if (k == 0) {
                    std:cout << chain_out[j][k] << std::endl;
                }
                }            
                accepted++;
            }
           
            /*if (j == 0) {
                for (int k = 0; k < 3; k++) {
                    std::cout << "output[" << accepted << "][" << k << "] = " << chain_out[accepted][k] << std::endl;
                }
            }*/

            for (int k = 0; k < 3; k++) {
                chain_in[j][k] = chain_out[j][k];
            }
        }
        
    }

}

} // namespace solver
} // namespace xf

#endif
