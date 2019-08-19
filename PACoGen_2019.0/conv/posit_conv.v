`timescale 1ns / 1ps
module posit_conv(clk, in1, start, out, done);

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

input [N-1:0] in1;
input start, clk;
output [63:0] out;
output done;

wire start0 = start;
wire s1 = in1[N-1];
wire rc1;
wire [Bs-1:0] regime1;
wire [es-1:0] e1;
wire [N-es-1:0] mant1;
wire zero_tmp1 = |in1[N-2:0];
wire [N-1:0] xin1 = s1 ? -in1 : in1;

data_extract_v1 #(.N(N),.es(es)) uut_de1(.in(xin1), .rc(rc1), .regime(regime1), .exp(e1), .mant(mant1));

wire [N-es:0] m1 = {zero_tmp1,mant1};
reg rc;
reg sign_1, start_1;
reg [Bs-1:0] regime_1;
reg [es-1:0] e1_1;
reg [N-es-1:0] mant1_1;

always @(posedge clk) begin
    rc <= rc1;
    regime_1 <= regime1;
    e1_1 <= e1;
    mant1_1 <= m1;
    sign_1 <= s1;
    start_1 <= start0;
end


wire [10:0] bias = 11'b01111111111;
wire [10:0] regi = {regime_1, {es{1'b0}}};
wire [10:0] exp = e1_1;
wire [10:0] fullexp = bias - regi + exp;

wire [51-(N-es):0] padding = 0;
wire [51:0] fullmant = {mant1_1, padding};

reg [10:0] fullexp_2;
reg [52:0] mant1_2;
reg sign_2, start_2;

always @(posedge clk) begin
    fullexp_2 <= fullexp;
    mant1_2 <= fullmant;
    sign_2 <= sign_1;
    start_2 <= start_1;
end

assign out[63] = sign_2;
assign out[62:52] = fullexp_2;
assign out[51:0] = mant1_2;
assign done = start_2;






endmodule



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