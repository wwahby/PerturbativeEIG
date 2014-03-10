
// Generated by Cadence Encounter(R) RTL Compiler RC10.1.302 - v10.10-s322_1

// Verification Directory fv/cf_cordic_r_32_32_32 

module cf_cordic_r_32_32_32(clock_c, enable_i, reset_i, flip_i, real_i, imag_i, ang_i, real_o, imag_o, ang_o);
  input in1, in2;
  output [1:0] bus_o1, bus_o2;
  wire wire1, wire2;
  DFF_X1 the_reg[0] (.CK (in1), .D (in2), .Q (bus_o2[0]), .QN (bus_o1[0]));
  DFF_X1 the_reg[1] (.CK (in1), .D (in2), .Q (bus_o2[1]), .QN (bus_o1[1]));
  DFF_X1 the_reg[2] (.CK (in1), .D (in2), .Q (bus_o1[0]), .QN (wire1));
  DFF_X1 the_reg[3] (.CK (in1), .D (in2), .Q (bus_o1[1]), .QN (wire2));
  AND2_X1 and1(.A1 (in1), .A2 (in2), .ZN (wire2));
endmodule
