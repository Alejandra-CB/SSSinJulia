function mpc = bs_2area4gen
mpc.version = '2';
mpc.baseMVA = 100.0;

%% area data
%	area	refbus
mpc.areas = [
	1	 4;
];

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.03000	    20.2000	 20.00	 1	    1.10000	    0.90000;
	2	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.01000	    10.5000	 20.00	 1	    1.10000	    0.90000;
	3	 3	 0.0	 0.0	 0.0	 0.0	 1	    1.03000	    -6.8000	 20.00	 1	    1.10000	    0.90000;
	4	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.01000	    -17.000	 20.00	 1	    1.10000	    0.90000;
	5	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	6	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	7	 1	 967.0	 100.0	 0.0	 200.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	8	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	9	 1	 1767.0	 100.0	 0.0	 350.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	10	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	11	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
    1	 700.0	 185.0	 300.0	 -300.0	 1.03	 100.0	 1	 800.0	 100.0;
	2	 700.0	 235.0	 300.0	 -300.0	 1.01	 100.0	 1	 800.0	 100.0;
	3	 719.0	 176.0	 300.0	 -300.0	 1.03	 100.0	 1	 800.0	 100.0;
	4	 700.0	 202.0	 300.0	 -300.0	 1.01	 100.0	 1	 800.0	 100.0;
];

%% generator cost data -> Not used in the PF formulation (I am not working with OPF)
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 0;
	2	 0.0	 0.0	 0;
	2	 0.0	 0.0	 0;
	2	 0.0	 0.0	 0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 5	 0.00000	 0.0166	 0.00000	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	2	 6	 0.00000	 0.0166	 0.00000	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	3	 11	 0.00000	 0.0166	 0.00000	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	4	 10	 0.00000	 0.0166	 0.00000	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	5	 6	 0.00250	 0.0250	 0.04375	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	11	 10	 0.00250	 0.0250	 0.04375	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	6	 7	 0.00100	 0.0100	 0.01750	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	10	 9	 0.00100	 0.0100	 0.01750	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	8	 9	 0.01100	 0.1100	 0.19250	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	8	 9	 0.01100	 0.1100	 0.19250	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	7	 8	 0.01100	 0.1100	 0.19250	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
	7	 8	 0.01100	 0.1100	 0.19250	 1e6	 1e6	 1e6	 1.0	 0.0	 1	 -30.0	 30.0;
];
