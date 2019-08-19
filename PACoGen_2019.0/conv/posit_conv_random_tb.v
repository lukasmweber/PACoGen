`timescale 1ns/1ps

module posit_conv_random_tb;

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

reg[N-1:0] in1;
reg start;

wire[63:0] out;
wire done;

reg[sample_bits-1:0] in_counter, out_counter;
reg clk, cen;

reg [N-1:0] data1 [1:samples];
reg [63:0] ref   [1:samples];

initial $readmemb("Pin1_36bit.txt", data1);
initial $readmemb("Pout_36bit_CONV.txt", ref);

posit_conv  #(.N(N), .es(es)) add (clk, in1, start, out, done);

initial begin
  in1 = 0;
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
      in_counter = in_counter + 1;
      start = 1;
    end
    else begin
      in1 = 0;
      start = 0;
    end
  end
end

reg [N-1:0] error;
wire [24:0] maxErrBits = {24{1'b1}};
always @(negedge clk) begin
  if(out_counter <= samples) begin
    if(done) begin
      error = (ref[out_counter] > out) ? ref[out_counter] - out : out - ref[out_counter];
      if(error>maxErrBits) begin
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
