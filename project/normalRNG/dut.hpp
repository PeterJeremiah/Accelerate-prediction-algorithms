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
/**
 *  @file mcmc_kernel.hpp
 *  @brief  Header file for kernel wrapper
 */

#ifndef _DUT_HPP_
#define _DUT_HPP_

#include "xf_fintech/rng.hpp"
#include <ap_fixed.h>
/// @brief Specific implementation of this kernel
#define NCHAINS 10
#define NSAMPLES_MAX 20000
#define Observation 200
//#define FIXED_DOUBLE ap_fixed<64, 11, AP_TRN, AP_SAT>
#define DT double

/**
* @brief Top level Kernel function.  \n
*
*@tparam DT                 - Data type used in whole kernel (double by default)
*@tparam NCHAINS            - Number of chains
*@tparam NSAMPLES_MAX       - Maximum Number of chains for synthesis purpose
*@param[in] temp_inv        - Array of Inverted temperatures of the chain that density is generate for (1/Temp)
*@param[in] sigma           - Array of sigmas for Proposal generation for each chain
*@param[in] nSamples        - Number of samples to generate
*@param[out] sample_output  - Sample output
*/
extern "C" void dut(DT temp_inv[NCHAINS],
                            DT sigma[NCHAINS],
                            DT sample_output[NSAMPLES_MAX][3],
                            unsigned int nSamples,
							DT y[NCHAINS][NSAMPLES_MAX],
							DT data[Observation][2]);

#endif
