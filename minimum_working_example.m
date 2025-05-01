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
	param.L  = 256*[1,1];
	% number of grid cells
	param.nx = param.L(1)*[1,1]/2;
	% final time
	param.T  = 365*10;
	% time step
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
	param.opt.dt_min       = 1/400;
	param.opt.dt_max       = 0.5;
	param.opt.outer_abstol = 1e-4;
	param.opt.outer_reltol = 1e-2;
	param.opt.dt_max_scale_up   = sqrt(2);
	param.opt.dt_min_scale_down = 0;
	% keep time step constant
	%param.opt.adapt_time_step = 0;
	% time step for writing output files
	param.opt.dto = 365;
	% data type of output file
	param.opt.compute_class = @single;
	% output directory
	param.opt.path_str = 'mat/';
	% solve using a splitting scheme
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_spectral';
	
	% random initial condition
	%param.initial_condition = 'obj.random_state()';
	param.initial_condition = 'obj.ic_single_patch()';
	rad = Rietkerk(param);

	[t,y]	  = rad.run();

	[b,w,h] = rad.extract2(y(end,:));

	subplot(2,3,1)
	imagesc(double(b))

	% analysis
	sp = Spatial_Pattern();
	sp.opt.suppress_low_frequency_components = 0;
	sp.L = rad.L;
	sp.b = b;
	
	% sp.source = a;
	sp.analyze_grid();
	%sp.fit_parametric_densities();
	% sp(kdx).predict_pattern();
	
	subplot(2,3,2);
	sp.plot('S.radial.con');

	subplot(2,3,3);
	sp.plot('R.radial.con');
	
