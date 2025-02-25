################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
D:/zu7/normalRNG/dut.cpp 

OBJS += \
./source/dut.o 

CPP_DEPS += \
./source/dut.d 


# Each subdirectory must supply rules for building sources it contributes
source/dut.o: D:/zu7/normalRNG/dut.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ID:/tmp/Vitis_Libraries/quantitative_finance//L1/include -ID:/Xilinx/Vitis_HLS/2023.2/include/ap_sysc -ID:/Xilinx/Vitis_HLS/2023.2/win64/tools/auto_cc/include -ID:/Xilinx/Vitis_HLS/2023.2/win64/tools/systemc/include -ID:/Xilinx/Vitis_HLS/2023.2/include -ID:/Xilinx/Vitis_HLS/2023.2/include/etc -ID:/zu7/normalRNG -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"source/dut.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


