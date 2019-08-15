vlib work

vlog posit_add_tb.v
vlog posit_add.v

vsim -voptargs=+acc -t ps posit_add_tb

add wave *
run -all
