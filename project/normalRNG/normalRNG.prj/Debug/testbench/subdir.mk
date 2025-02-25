################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
D:/zu7/normalRNG/tb.cpp 

OBJS += \
./testbench/tb.o 

CPP_DEPS += \
./testbench/tb.d 


# Each subdirectory must supply rules for building sources it contributes
testbench/tb.o: D:/zu7/normalRNG/tb.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ID:/tmp/Vitis_Libraries/quantitative_finance//L1/include -ID:/Xilinx/Vitis_HLS/2023.2/include/ap_sysc -ID:/Xilinx/Vitis_HLS/2023.2/win64/tools/auto_cc/include -ID:/Xilinx/Vitis_HLS/2023.2/win64/tools/systemc/include -ID:/Xilinx/Vitis_HLS/2023.2/include -ID:/Xilinx/Vitis_HLS/2023.2/include/etc -ID:/zu7/normalRNG -O0 -g3 -Wall -Wno-unknown-pragmas -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"testbench/tb.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


