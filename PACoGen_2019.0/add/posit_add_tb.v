`timescale 1ns/1ps

module posit_add_tb;

parameter N=33;
parameter es=5;

reg [N-1:0] in1, in2;
reg start;

wire [N-1:0] out;
wire inf, zero, done;

posit_add #(.N(N), .es(es)) add (in1, in2, start, out, inf, zero, done);

initial begin
  in1 = 0;
  in2 = 0;
  start = 0;
  
  #10

  in1 = {10'b0001_0000_01, {(N-10){1'b0}}};
  in2 = {10'b0001_0000_01, {(N-10){1'b0}}};
  start = 1'b1;

  #100

  $finish;

end

endmodule