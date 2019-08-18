`timescale 1ns/1ps

module posit_add_random_tb;

function [31:0] log2;
input reg [31:0] value;
	begin
	value = value-1;
	for (log2=0; value>0; log2=log2+1)
        	value = value>>1;
      	end
endfunction

parameter N=36;
parameter es=5;
parameter samples=10000;
parameter sample_bits = log2(samples);

reg[N-1:0] in1, in2;
reg start;

wire[N-1:0] out;
wire inf, zero, done;

reg[sample_bits-1:0] in_counter, out_counter;
reg clk, cen;

reg [N-1:0] data1 [1:samples];
reg [N-1:0] data2 [1:samples];
reg [N-1:0] ref   [1:samples];

initial $readmemb("Pin1_36bit.txt", data1);
initial $readmemb("Pin2_36bit.txt", data2);
initial $readmemb("Pout_36bit_MUL.txt", ref);

posit_mult  #(.N(N), .es(es)) add (clk, in1, in2, start, out, done);

initial begin
  in1 = 0;
  in2 = 0;
  start = 0;
  in_counter = 1;
  out_counter = 1;
  clk = 0;
  cen = 0;
  
  #100;

  cen = 1;
  //#10000 $finish ;

end

always begin
  clk = ~clk;
  #5;
end

always @(posedge clk) begin
  if(cen) begin  
    if(in_counter <= samples) begin
      in1 = data1[in_counter];
      in2 = data2[in_counter];
      in_counter = in_counter + 1;
      start = 1;
    end
    else begin
      in1 = 0;
      in2 = 0;
      start = 0;
    end
  end
end

reg [N-1:0] error;
always @(negedge clk) begin
  if(out_counter <= samples) begin
    if(done) begin
      error = (ref[out_counter] > out) ? ref[out_counter] - out : out - ref[out_counter];
      if(error>1) begin
        $display("Computation error %h vs %h\n", out, ref[out_counter]);
      end
      out_counter = out_counter + 1;
    end
  end
  else begin
    $finish;
  end
end
endmodule
