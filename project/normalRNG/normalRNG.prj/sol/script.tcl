############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
## Copyright 2022-2023 Advanced Micro Devices, Inc. All Rights Reserved.
############################################################
open_project normalRNG.prj
set_top dut
add_files dut.cpp -cflags "-ID:/tmp/Vitis_Libraries/quantitative_finance//L1/include"
add_files -tb tb.cpp -cflags "-ID:/tmp/Vitis_Libraries/quantitative_finance//L1/include -Wno-unknown-pragmas"
open_solution "sol" -flow_target vitis
set_part {xczu7eg-ffvc1156-2-i}
create_clock -period 300MHz -name default
config_export -format ip_catalog -rtl verilog
config_interface -m_axi_alignment_byte_size 64 -m_axi_latency 64 -m_axi_max_widen_bitwidth 512
config_rtl -register_reset_num 3
source "./normalRNG.prj/sol/directives.tcl"
csim_design
csynth_design
cosim_design
export_design -rtl verilog -format ip_catalog
