% Sat 31 Aug 12:49:37 CEST 2024
% Karl Kastner, Berlin

	% addpath_recursive('lib/');
	mkdir('mat/')

	% model parameters
	param        = struct();
	% advection in x-direction
	param.pmu.vx = [0, 0, 0];
	% diffusion in x-direction
	param.pmu.ex = [0.1, 0.1, 100];
	% precipitation
	param.pmu.R  = 0.75;
	% mean infiltration coefficient
	param.pmu.a = 0.2;
	% spatial heterogeneity of parameter a
	param.pss.a  = 0;
	% boundary condition
	param.boundary_condition    = {'circular','ciruclar'};
	% reload values of intermediate time steps
	param.opt.loadfinal = false;
	% domain size
	param.L  = 1024*[1,1];
	% number of grid cells
	param.nx = param.L(1)*[1,1]/2;
	% final time
	param.T  = 365*100;
	% time step
	param.opt.dt = 1;
	% keep time step constant
	param.opt.adapt_time_step = 0;
	% time step for writing output files
	param.opt.dto = 365;
	% data type of output file
	param.opt.compute_class = @single;
	% output directory
	param.opt.path_str = 'mat/';
	% solve using a splitting scheme
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_implicit_q_fft';
	
	% random initial condition
	%param.initial_condition = 'obj.random_state()';
	param.initial_condition = 'obj.ic_single_patch()';
	rk = Rietkerk(param);

	[t,y]	  = rk.run();

	[b,w,h] = rk.extract2(y(end,:));

	imagesc(double(b))

