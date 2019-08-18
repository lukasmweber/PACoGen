`timescale 1ns / 1ps
//(* use_dsp = "no" *)
module posit_mult(clk, in1, in2, start, out, done);

function [31:0] log2;
input reg [31:0] value;
	begin
	value = value-1;
	for (log2=0; value>0; log2=log2+1)
        	value = value>>1;
      	end
endfunction

parameter N = 36;
parameter Bs = log2(N); 
parameter es = 5;

input clk;
input [N-1:0] in1, in2;
input start; 
output [N-1:0] out;
output done;

wire start0= start;
wire s1 = in1[N-1];
wire s2 = in2[N-1];
wire zero_tmp1 = |in1[N-2:0];
wire zero_tmp2 = |in2[N-2:0];
wire inf1 = in1[N-1] & (~zero_tmp1),
	inf2 = in2[N-1] & (~zero_tmp2);
wire zero1 = ~(in1[N-1] | zero_tmp1),
	zero2 = ~(in2[N-1] | zero_tmp2);
wire inf = inf1 | inf2,
	zero = zero1 & zero2;

//Data Extraction
wire rc1, rc2;
wire [Bs-1:0] regime1, regime2;
wire [es-1:0] e1, e2;
wire [N-es-1:0] mant1, mant2;
wire [N-1:0] xin1 = s1 ? -in1 : in1;
wire [N-1:0] xin2 = s2 ? -in2 : in2;
data_extract_v1 #(.N(N),.es(es)) uut_de1(.in(xin1), .rc(rc1), .regime(regime1), .exp(e1), .mant(mant1));
data_extract_v1 #(.N(N),.es(es)) uut_de2(.in(xin2), .rc(rc2), .regime(regime2), .exp(e2), .mant(mant2));

wire [N-es:0] m1 = {zero_tmp1,mant1}, 
	m2 = {zero_tmp2,mant2};
wire mult_s = s1 ^ s2;


wire rc1_6, rc2_6, start_6, inf_6, zero_6, mult_s_6;
wire [Bs-1:0] regime1_6, regime2_6;
wire [es-1:0] e1_6, e2_6;

SHIFT_REGISTER #(.SIZE(1), .LEN(6)) sr_rc1 (.clk(clk), .in(rc1), .out(rc1_6));
SHIFT_REGISTER #(.SIZE(1), .LEN(6)) sr_rc2 (.clk(clk), .in(rc2), .out(rc2_6));
SHIFT_REGISTER #(.SIZE(1), .LEN(6)) sr_inf(.clk(clk), .in(inf), .out(inf_6));
SHIFT_REGISTER #(.SIZE(1), .LEN(6)) sr_zero (.clk(clk), .in(zero), .out(zero_6));
SHIFT_REGISTER #(.SIZE(1), .LEN(6)) sr_mult_s (.clk(clk), .in(mult_s), .out(mult_s_6));
SHIFT_REGISTER #(.SIZE(1), .LEN(6)) sr_start (.clk(clk), .in(start0), .out(start_6));

SHIFT_REGISTER #(.SIZE(Bs), .LEN(6)) sr_regime1 (.clk(clk), .in(regime1), .out(regime1_6));
SHIFT_REGISTER #(.SIZE(Bs), .LEN(6)) sr_regime2(.clk(clk), .in(regime2), .out(regime2_6));

SHIFT_REGISTER #(.SIZE(es), .LEN(6)) sr_e1 (.clk(clk), .in(e1), .out(e1_6));
SHIFT_REGISTER #(.SIZE(es), .LEN(6)) sr_e2 (.clk(clk), .in(e2), .out(e2_6));


reg [N-es:0] m1_1, m2_1;
always @(posedge clk) begin
    m1_1 <= m1;
    m2_1 <= m2;
end


//Sign, Exponent and Mantissa Computation
wire [2*(N-es)+1:0] mult_m;
wire rst = 0;
DSPMultv2 mult (.clock(clk), .reset(rst), .io_a(m1_1), .io_b(m2_1), .io_r(mult_m)); 

wire mult_m_ovf = mult_m[2*(N-es)+1];
wire [2*(N-es)+1:0] mult_mN = ~mult_m_ovf ? mult_m << 1'b1 : mult_m;

wire [Bs+1:0] r1 = rc1_6 ? {2'b0,regime1_6} : -regime1_6;
wire [Bs+1:0] r2 = rc2_6 ? {2'b0,regime2_6} : -regime2_6;
wire [Bs+es+1:0] mult_e;
add_N_Cin #(.N(Bs+es+1)) uut_add_exp ({r1,e1_6}, {r2,e2_6}, mult_m_ovf, mult_e);

//Exponent and Regime Computation
wire [es-1:0] e_o;
wire [Bs:0] r_o;
reg_exp_op #(.es(es), .Bs(Bs)) uut_reg_ro (mult_e[es+Bs+1:0], e_o, r_o);

//Exponent, Mantissa and GRS Packing
wire [2*N-1+3:0] tmp_o = {{N{~mult_e[es+Bs+1]}},mult_e[es+Bs+1],e_o,mult_mN[2*(N-es):2*(N-es)-(N-es-1)+1], mult_mN[2*(N-es)-(N-es-1):2*(N-es)-(N-es-1)-1], |mult_mN[2*(N-es)-(N-es-1)-2:0] }; 

reg start_7, inf_7, zero_7;
reg [es-1:0] e_o_7;
reg [Bs:0] r_o_7;
reg mult_s_7;
reg [2*N-1+3:0] tmp_o_7;
reg [2*(N-es)+1:0] mult_mN_7;

always @(posedge clk) begin
    start_7 <= start_6;
    inf_7 <= inf_6;
    zero_7 <= zero_6;
    mult_s_7 <= mult_s_6;
    tmp_o_7 <= tmp_o;
    e_o_7 <= e_o;
    r_o_7 <= r_o;
    mult_mN_7 <= mult_mN;
end



//Including Regime bits in Exponent-Mantissa Packing
wire [3*N-1+3:0] tmp1_o;
DSR_right_N_S #(.N(3*N+3), .S(Bs+1)) dsr2 (.a({tmp_o_7,{N{1'b0}}}), .b(r_o_7[Bs] ? {Bs{1'b1}} : r_o_7), .c(tmp1_o));

reg [3*N-1+3:0] tmp1_o_8;
reg start_8, inf_8, zero_8;
reg [Bs:0] r_o_8;
reg mult_s_8;
reg [2*(N-es)+1:0] mult_mN_8;

always @(posedge clk) begin
    start_8 <= start_7;
    inf_8 <= inf_7;
    zero_8 <= zero_7;
    mult_s_8 <= mult_s_7;
    r_o_8 <= r_o_7;
    mult_mN_8 <= mult_mN_7;
    tmp1_o_8 <= tmp1_o;
end



//Rounding RNE : ulp_add = G.(R + S) + L.G.(~(R+S))
wire L = tmp1_o_8[N+4], G = tmp1_o_8[N+3], R = tmp1_o_8[N+2], St = |tmp1_o_8[N+1:0],
     ulp = ((G & (R | St)) | (L & G & ~(R | St)));
wire [N-1:0] rnd_ulp = {{N-1{1'b0}},ulp};

wire [N:0] tmp1_o_rnd_ulp;
add_N #(.N(N)) uut_add_ulp (tmp1_o_8[2*N-1+3:N+3], rnd_ulp, tmp1_o_rnd_ulp);
wire [N-1:0] tmp1_o_rnd = (r_o_8 < N-es-2) ? tmp1_o_rnd_ulp[N-1:0] : tmp1_o_8[2*N-1+3:N+3];


//Final Output
wire [N-1:0] tmp1_oN = mult_s_8 ? -tmp1_o_rnd : tmp1_o_rnd;
wire done_w = start_8;
wire [N-1:0] output_w = inf_8|zero_8|(~mult_mN_8[2*(N-es)+1]) ? {inf_8,{N-1{1'b0}}} : {mult_s_8, tmp1_oN[N-1:1]};
	
reg [N-1:0] out_9;
reg done_9;

always @(posedge clk) begin
    out_9 <= output_w;
    done_9 <= done_w;
end

assign out = out_9;
assign done = done_9;

endmodule

/////////////////////////
module data_extract_v1(in, rc, regime, exp, mant);

function [31:0] log2;
input reg [31:0] value;
	begin
	value = value-1;
	for (log2=0; value>0; log2=log2+1)
        	value = value>>1;
      	end
endfunction

parameter N=16;
parameter Bs=log2(N);
parameter es = 2;
input [N-1:0] in;
output rc;
output [Bs-1:0] regime;
output [es-1:0] exp;
output [N-es-1:0] mant;

wire [N-1:0] xin = in;
assign rc = xin[N-2];

wire [N-1:0] xin_r = rc ? ~xin : xin;

wire [Bs-1:0] k;
LOD_N #(.N(N)) xinst_k(.in({xin_r[N-2:0],rc^1'b0}), .out(k));

assign regime = rc ? k-1 : k;

wire [N-1:0] xin_tmp;
DSR_left_N_S #(.N(N), .S(Bs)) ls (.a({xin[N-3:0],2'b0}),.b(k),.c(xin_tmp));

assign exp= xin_tmp[N-1:N-es];
assign mant= xin_tmp[N-es-1:0];

endmodule


/////////////////
module sub_N (a,b,c);
parameter N=10;
input [N-1:0] a,b;
output [N:0] c;
assign c = {1'b0,a} - {1'b0,b};
endmodule

/////////////////////////
module add_N (a,b,c);
parameter N=10;
input [N-1:0] a,b;
output [N:0] c;
assign c = {1'b0,a} + {1'b0,b};
endmodule

/////////////////////////
module add_N_Cin (a,b,cin,c);
parameter N=10;
input [N:0] a,b;
input cin;
output [N:0] c;
assign c = a + b + cin;
endmodule


/////////////////////////
module add_1 (a,mant_ovf,c);
parameter N=10;
input [N:0] a;
input mant_ovf;
output [N:0] c;
assign c = a + mant_ovf;
endmodule

/////////////////////////
module conv_2c (a,c);
parameter N=10;
input [N:0] a;
output [N:0] c;
assign c = a + 1'b1;
endmodule

/////////////////////////
module reg_exp_op (exp_o, e_o, r_o);
parameter es=3;
parameter Bs=5;
input [es+Bs+1:0] exp_o;
output [es-1:0] e_o;
output [Bs:0] r_o;

assign e_o = exp_o[es-1:0];

wire [es+Bs:0] exp_oN_tmp;
conv_2c #(.N(es+Bs)) uut_conv_2c1 (~exp_o[es+Bs:0],exp_oN_tmp);
wire [es+Bs:0] exp_oN = exp_o[es+Bs+1] ? exp_oN_tmp[es+Bs:0] : exp_o[es+Bs:0];

assign r_o = (~exp_o[es+Bs+1] || |(exp_oN[es-1:0])) ? exp_oN[es+Bs:es] + 1 : exp_oN[es+Bs:es];
endmodule

/////////////////////////
module DSR_left_N_S(a,b,c);
        parameter N=16;
        parameter S=4;
        input [N-1:0] a;
        input [S-1:0] b;
        output [N-1:0] c;

wire [N-1:0] tmp [S-1:0];
assign tmp[0]  = b[0] ? a << 7'd1  : a; 
genvar i;
generate
	for (i=1; i<S; i=i+1)begin:loop_blk
		assign tmp[i] = b[i] ? tmp[i-1] << 2**i : tmp[i-1];
	end
endgenerate
assign c = tmp[S-1];

endmodule


/////////////////////////
module DSR_right_N_S(a,b,c);
        parameter N=16;
        parameter S=4;
        input [N-1:0] a;
        input [S-1:0] b;
        output [N-1:0] c;

wire [N-1:0] tmp [S-1:0];
assign tmp[0]  = b[0] ? a >> 7'd1  : a; 
genvar i;
generate
	for (i=1; i<S; i=i+1)begin:loop_blk
		assign tmp[i] = b[i] ? tmp[i-1] >> 2**i : tmp[i-1];
	end
endgenerate
assign c = tmp[S-1];

endmodule

/////////////////////////

module SHIFT_REGISTER (clk, in, out);

parameter SIZE = 32;
parameter LEN = 4;

input clk;
input [SIZE-1:0] in;
output [SIZE-1:0] out;

reg [SIZE-1:0] shiftr [LEN-1:0];

integer index;
always @(posedge clk) begin
    shiftr[0] <= in;

    for(index= 0; index < LEN - 1; index = index + 1) begin
        shiftr[index + 1] <= shiftr[index];
    end
end

assign out = shiftr[LEN-1];

endmodule


module LOD_N (in, out);

  function [31:0] log2;
    input reg [31:0] value;
    begin
      value = value-1;
      for (log2=0; value>0; log2=log2+1)
	value = value>>1;
    end
  endfunction

parameter N = 64;
parameter S = log2(N); 
input [N-1:0] in;
output [S-1:0] out;

wire vld;
LOD #(.N(N)) l1 (in, out, vld);
endmodule


module LOD (in, out, vld);

  function [31:0] log2;
    input reg [31:0] value;
    begin
      value = value-1;
      for (log2=0; value>0; log2=log2+1)
	value = value>>1;
    end
  endfunction


parameter N = 64;
parameter S = log2(N);

   input [N-1:0] in;
   output [S-1:0] out;
   output vld;

  generate
    if (N == 2)
      begin
	assign vld = |in;
	assign out = ~in[1] & in[0];
      end
    else if (N & (N-1))
      //LOD #(1<<S) LOD ({1<<S {1'b0}} | in,out,vld);
      LOD #(1<<S) LOD ({in,{((1<<S) - N) {1'b0}}},out,vld);
    else
      begin
	wire [S-2:0] out_l, out_h;
	wire out_vl, out_vh;
	LOD #(N>>1) l(in[(N>>1)-1:0],out_l,out_vl);
	LOD #(N>>1) h(in[N-1:N>>1],out_h,out_vh);
	assign vld = out_vl | out_vh;
	assign out = out_vh ? {1'b0,out_h} : {out_vl,out_l};
      end
  endgenerate
endmodule

module DSPMult24x17( // @[:@3.2]
  input         clock, // @[:@4.4]
  input  [23:0] io_a, // @[:@6.4]
  input  [16:0] io_b, // @[:@6.4]
  output [40:0] io_r // @[:@6.4]
);
  reg [40:0] rOut; // @[DSPMultv2.scala 79:21:@9.4]
  reg [63:0] _RAND_0;
  wire [23:0] _GEN_0; // @[DSPMultv2.scala 79:27:@8.4]
  wire [40:0] _T_11; // @[DSPMultv2.scala 79:27:@8.4]
  assign _GEN_0 = {{7'd0}, io_b}; // @[DSPMultv2.scala 79:27:@8.4]
  assign _T_11 = io_a * _GEN_0; // @[DSPMultv2.scala 79:27:@8.4]
  assign io_r = rOut;
`ifdef RANDOMIZE_GARBAGE_ASSIGN
`define RANDOMIZE
`endif
`ifdef RANDOMIZE_INVALID_ASSIGN
`define RANDOMIZE
`endif
`ifdef RANDOMIZE_REG_INIT
`define RANDOMIZE
`endif
`ifdef RANDOMIZE_MEM_INIT
`define RANDOMIZE
`endif
`ifdef RANDOMIZE
  integer initvar;
  initial begin
    `ifndef verilator
      #0.002 begin end
    `endif
  `ifdef RANDOMIZE_REG_INIT
  _RAND_0 = {2{$random}};
  rOut = _RAND_0[40:0];
  `endif // RANDOMIZE_REG_INIT
  end
`endif // RANDOMIZE
  always @(posedge clock) begin
    rOut <= _T_11;
  end
endmodule
module PipelinedAdderv2( // @[:@43.2]
  input         clock, // @[:@44.4]
  input  [63:0] io_a, // @[:@46.4]
  input  [63:0] io_b, // @[:@46.4]
  output [64:0] io_r // @[:@46.4]
);
  reg [31:0] aDelayed_1; // @[Reg.scala 11:16:@52.4]
  reg [31:0] _RAND_0;
  reg [31:0] bDelayed_1; // @[Reg.scala 11:16:@56.4]
  reg [31:0] _RAND_1;
  reg [32:0] _T_19; // @[PipelinedAdderv2.scala 82:25:@63.4]
  reg [63:0] _RAND_2;
  reg [32:0] _T_24; // @[PipelinedAdderv2.scala 82:25:@70.4]
  reg [63:0] _RAND_3;
  reg [31:0] flushRes_0; // @[Reg.scala 11:16:@74.4]
  reg [31:0] _RAND_4;
  wire [31:0] a_0; // @[PipelinedAdderv2.scala 39:53:@48.4]
  wire [31:0] a_1; // @[PipelinedAdderv2.scala 39:53:@49.4]
  wire [31:0] b_0; // @[PipelinedAdderv2.scala 40:53:@50.4]
  wire [31:0] b_1; // @[PipelinedAdderv2.scala 40:53:@51.4]
  wire [32:0] _T_15; // @[PipelinedAdderv2.scala 82:28:@60.4]
  wire [33:0] _T_16; // @[PipelinedAdderv2.scala 82:33:@61.4]
  wire [32:0] _T_17; // @[PipelinedAdderv2.scala 82:33:@62.4]
  wire [31:0] resultSignals_0; // @[PipelinedAdderv2.scala 83:12:@65.4]
  wire  carrySignals_1; // @[PipelinedAdderv2.scala 83:39:@66.4]
  wire [32:0] _T_20; // @[PipelinedAdderv2.scala 82:28:@67.4]
  wire [32:0] _GEN_3; // @[PipelinedAdderv2.scala 82:33:@68.4]
  wire [33:0] _T_21; // @[PipelinedAdderv2.scala 82:33:@68.4]
  wire [32:0] _T_22; // @[PipelinedAdderv2.scala 82:33:@69.4]
  wire [31:0] flushRes_1; // @[PipelinedAdderv2.scala 83:12:@72.4]
  wire  carrySignals_2; // @[PipelinedAdderv2.scala 83:39:@73.4]
  wire [63:0] _T_27; // @[Cat.scala 30:58:@78.4]
  wire [64:0] _T_28; // @[Cat.scala 30:58:@79.4]
  assign a_0 = io_a[31:0]; // @[PipelinedAdderv2.scala 39:53:@48.4]
  assign a_1 = io_a[63:32]; // @[PipelinedAdderv2.scala 39:53:@49.4]
  assign b_0 = io_b[31:0]; // @[PipelinedAdderv2.scala 40:53:@50.4]
  assign b_1 = io_b[63:32]; // @[PipelinedAdderv2.scala 40:53:@51.4]
  assign _T_15 = a_0 + b_0; // @[PipelinedAdderv2.scala 82:28:@60.4]
  assign _T_16 = _T_15 + 33'h0; // @[PipelinedAdderv2.scala 82:33:@61.4]
  assign _T_17 = _T_16[32:0]; // @[PipelinedAdderv2.scala 82:33:@62.4]
  assign resultSignals_0 = _T_19[31:0]; // @[PipelinedAdderv2.scala 83:12:@65.4]
  assign carrySignals_1 = _T_19[32]; // @[PipelinedAdderv2.scala 83:39:@66.4]
  assign _T_20 = aDelayed_1 + bDelayed_1; // @[PipelinedAdderv2.scala 82:28:@67.4]
  assign _GEN_3 = {{32'd0}, carrySignals_1}; // @[PipelinedAdderv2.scala 82:33:@68.4]
  assign _T_21 = _T_20 + _GEN_3; // @[PipelinedAdderv2.scala 82:33:@68.4]
  assign _T_22 = _T_21[32:0]; // @[PipelinedAdderv2.scala 82:33:@69.4]
  assign flushRes_1 = _T_24[31:0]; // @[PipelinedAdderv2.scala 83:12:@72.4]
  assign carrySignals_2 = _T_24[32]; // @[PipelinedAdderv2.scala 83:39:@73.4]
  assign _T_27 = {flushRes_1,flushRes_0}; // @[Cat.scala 30:58:@78.4]
  assign _T_28 = {carrySignals_2,_T_27}; // @[Cat.scala 30:58:@79.4]
  assign io_r = _T_28;
`ifdef RANDOMIZE_GARBAGE_ASSIGN
`define RANDOMIZE
`endif
`ifdef RANDOMIZE_INVALID_ASSIGN
`define RANDOMIZE
`endif
`ifdef RANDOMIZE_REG_INIT
`define RANDOMIZE
`endif
`ifdef RANDOMIZE_MEM_INIT
`define RANDOMIZE
`endif
`ifdef RANDOMIZE
  integer initvar;
  initial begin
    `ifndef verilator
      #0.002 begin end
    `endif
  `ifdef RANDOMIZE_REG_INIT
  _RAND_0 = {1{$random}};
  aDelayed_1 = _RAND_0[31:0];
  `endif // RANDOMIZE_REG_INIT
  `ifdef RANDOMIZE_REG_INIT
  _RAND_1 = {1{$random}};
  bDelayed_1 = _RAND_1[31:0];
  `endif // RANDOMIZE_REG_INIT
  `ifdef RANDOMIZE_REG_INIT
  _RAND_2 = {2{$random}};
  _T_19 = _RAND_2[32:0];
  `endif // RANDOMIZE_REG_INIT
  `ifdef RANDOMIZE_REG_INIT
  _RAND_3 = {2{$random}};
  _T_24 = _RAND_3[32:0];
  `endif // RANDOMIZE_REG_INIT
  `ifdef RANDOMIZE_REG_INIT
  _RAND_4 = {1{$random}};
  flushRes_0 = _RAND_4[31:0];
  `endif // RANDOMIZE_REG_INIT
  end
`endif // RANDOMIZE
  always @(posedge clock) begin
    aDelayed_1 <= a_1;
    bDelayed_1 <= b_1;
    _T_19 <= _T_17;
    _T_24 <= _T_22;
    flushRes_0 <= resultSignals_0;
  end
endmodule
module DSPMultv2( // @[:@160.2]
  input         clock, // @[:@161.4]
  input         reset, // @[:@162.4]
  input  [31:0] io_a, // @[:@163.4]
  input  [31:0] io_b, // @[:@163.4]
  output [63:0] io_r // @[:@163.4]
);
  wire  DSPMult24x17_clock; // @[DSPMultv2.scala 186:22:@167.4]
  wire [23:0] DSPMult24x17_io_a; // @[DSPMultv2.scala 186:22:@167.4]
  wire [16:0] DSPMult24x17_io_b; // @[DSPMultv2.scala 186:22:@167.4]
  wire [40:0] DSPMult24x17_io_r; // @[DSPMultv2.scala 186:22:@167.4]
  wire  DSPMult24x17_1_clock; // @[DSPMultv2.scala 186:22:@175.4]
  wire [23:0] DSPMult24x17_1_io_a; // @[DSPMultv2.scala 186:22:@175.4]
  wire [16:0] DSPMult24x17_1_io_b; // @[DSPMultv2.scala 186:22:@175.4]
  wire [40:0] DSPMult24x17_1_io_r; // @[DSPMultv2.scala 186:22:@175.4]
  wire  DSPMult24x17_2_clock; // @[DSPMultv2.scala 186:22:@183.4]
  wire [23:0] DSPMult24x17_2_io_a; // @[DSPMultv2.scala 186:22:@183.4]
  wire [16:0] DSPMult24x17_2_io_b; // @[DSPMultv2.scala 186:22:@183.4]
  wire [40:0] DSPMult24x17_2_io_r; // @[DSPMultv2.scala 186:22:@183.4]
  wire  DSPMult24x17_3_clock; // @[DSPMultv2.scala 186:22:@190.4]
  wire [23:0] DSPMult24x17_3_io_a; // @[DSPMultv2.scala 186:22:@190.4]
  wire [16:0] DSPMult24x17_3_io_b; // @[DSPMultv2.scala 186:22:@190.4]
  wire [40:0] DSPMult24x17_3_io_r; // @[DSPMultv2.scala 186:22:@190.4]
  wire  PipelinedAdderv2_clock; // @[PipelinedAdderv2.scala 98:21:@196.4]
  wire [63:0] PipelinedAdderv2_io_a; // @[PipelinedAdderv2.scala 98:21:@196.4]
  wire [63:0] PipelinedAdderv2_io_b; // @[PipelinedAdderv2.scala 98:21:@196.4]
  wire [64:0] PipelinedAdderv2_io_r; // @[PipelinedAdderv2.scala 98:21:@196.4]
  wire  PipelinedAdderv2_1_clock; // @[PipelinedAdderv2.scala 98:21:@201.4]
  wire [63:0] PipelinedAdderv2_1_io_a; // @[PipelinedAdderv2.scala 98:21:@201.4]
  wire [63:0] PipelinedAdderv2_1_io_b; // @[PipelinedAdderv2.scala 98:21:@201.4]
  wire [64:0] PipelinedAdderv2_1_io_r; // @[PipelinedAdderv2.scala 98:21:@201.4]
  wire  PipelinedAdderv2_2_clock; // @[PipelinedAdderv2.scala 98:21:@206.4]
  wire [63:0] PipelinedAdderv2_2_io_a; // @[PipelinedAdderv2.scala 98:21:@206.4]
  wire [63:0] PipelinedAdderv2_2_io_b; // @[PipelinedAdderv2.scala 98:21:@206.4]
  wire [64:0] PipelinedAdderv2_2_io_r; // @[PipelinedAdderv2.scala 98:21:@206.4]
  wire [16:0] _T_11; // @[DSPMultv2.scala 121:16:@165.4]
  wire [7:0] _T_12; // @[DSPMultv2.scala 122:16:@166.4]
  wire [64:0] _T_14; // @[Cat.scala 30:58:@172.4]
  wire [14:0] _T_15; // @[DSPMultv2.scala 121:16:@173.4]
  wire [14:0] _T_16; // @[DSPMultv2.scala 122:16:@174.4]
  wire [74:0] _T_18; // @[Cat.scala 30:58:@180.4]
  wire [23:0] _T_20; // @[DSPMultv2.scala 122:16:@182.4]
  wire [16:0] _T_22; // @[DSPMultv2.scala 122:16:@189.4]
  wire [57:0] _T_24; // @[Cat.scala 30:58:@195.4]
  DSPMult24x17 DSPMult24x17 ( // @[DSPMultv2.scala 186:22:@167.4]
    .clock(DSPMult24x17_clock),
    .io_a(DSPMult24x17_io_a),
    .io_b(DSPMult24x17_io_b),
    .io_r(DSPMult24x17_io_r)
  );
  DSPMult24x17 DSPMult24x17_1 ( // @[DSPMultv2.scala 186:22:@175.4]
    .clock(DSPMult24x17_1_clock),
    .io_a(DSPMult24x17_1_io_a),
    .io_b(DSPMult24x17_1_io_b),
    .io_r(DSPMult24x17_1_io_r)
  );
  DSPMult24x17 DSPMult24x17_2 ( // @[DSPMultv2.scala 186:22:@183.4]
    .clock(DSPMult24x17_2_clock),
    .io_a(DSPMult24x17_2_io_a),
    .io_b(DSPMult24x17_2_io_b),
    .io_r(DSPMult24x17_2_io_r)
  );
  DSPMult24x17 DSPMult24x17_3 ( // @[DSPMultv2.scala 186:22:@190.4]
    .clock(DSPMult24x17_3_clock),
    .io_a(DSPMult24x17_3_io_a),
    .io_b(DSPMult24x17_3_io_b),
    .io_r(DSPMult24x17_3_io_r)
  );
  PipelinedAdderv2 PipelinedAdderv2 ( // @[PipelinedAdderv2.scala 98:21:@196.4]
    .clock(PipelinedAdderv2_clock),
    .io_a(PipelinedAdderv2_io_a),
    .io_b(PipelinedAdderv2_io_b),
    .io_r(PipelinedAdderv2_io_r)
  );
  PipelinedAdderv2 PipelinedAdderv2_1 ( // @[PipelinedAdderv2.scala 98:21:@201.4]
    .clock(PipelinedAdderv2_1_clock),
    .io_a(PipelinedAdderv2_1_io_a),
    .io_b(PipelinedAdderv2_1_io_b),
    .io_r(PipelinedAdderv2_1_io_r)
  );
  PipelinedAdderv2 PipelinedAdderv2_2 ( // @[PipelinedAdderv2.scala 98:21:@206.4]
    .clock(PipelinedAdderv2_2_clock),
    .io_a(PipelinedAdderv2_2_io_a),
    .io_b(PipelinedAdderv2_2_io_b),
    .io_r(PipelinedAdderv2_2_io_r)
  );
  assign _T_11 = io_a[16:0]; // @[DSPMultv2.scala 121:16:@165.4]
  assign _T_12 = io_b[31:24]; // @[DSPMultv2.scala 122:16:@166.4]
  assign _T_14 = {DSPMult24x17_io_r,24'h0}; // @[Cat.scala 30:58:@172.4]
  assign _T_15 = io_a[31:17]; // @[DSPMultv2.scala 121:16:@173.4]
  assign _T_16 = io_b[31:17]; // @[DSPMultv2.scala 122:16:@174.4]
  assign _T_18 = {DSPMult24x17_1_io_r,34'h0}; // @[Cat.scala 30:58:@180.4]
  assign _T_20 = io_b[23:0]; // @[DSPMultv2.scala 122:16:@182.4]
  assign _T_22 = io_b[16:0]; // @[DSPMultv2.scala 122:16:@189.4]
  assign _T_24 = {DSPMult24x17_3_io_r,17'h0}; // @[Cat.scala 30:58:@195.4]
  assign io_r = PipelinedAdderv2_2_io_r[63:0];
  assign DSPMult24x17_clock = clock;
  assign DSPMult24x17_io_a = {{16'd0}, _T_12};
  assign DSPMult24x17_io_b = _T_11;
  assign DSPMult24x17_1_clock = clock;
  assign DSPMult24x17_1_io_a = {{9'd0}, _T_15};
  assign DSPMult24x17_1_io_b = {{2'd0}, _T_16};
  assign DSPMult24x17_2_clock = clock;
  assign DSPMult24x17_2_io_a = _T_20;
  assign DSPMult24x17_2_io_b = _T_11;
  assign DSPMult24x17_3_clock = clock;
  assign DSPMult24x17_3_io_a = {{9'd0}, _T_15};
  assign DSPMult24x17_3_io_b = _T_22;
  assign PipelinedAdderv2_clock = clock;
  assign PipelinedAdderv2_io_a = _T_14[63:0];
  assign PipelinedAdderv2_io_b = _T_18[63:0];
  assign PipelinedAdderv2_1_clock = clock;
  assign PipelinedAdderv2_1_io_a = {{23'd0}, DSPMult24x17_2_io_r};
  assign PipelinedAdderv2_1_io_b = {{6'd0}, _T_24};
  assign PipelinedAdderv2_2_clock = clock;
  assign PipelinedAdderv2_2_io_a = PipelinedAdderv2_io_r[63:0];
  assign PipelinedAdderv2_2_io_b = PipelinedAdderv2_1_io_r[63:0];
endmodule

