% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% generated 2D patterns with the rietkerk model for varying degrees of
% spatial heterogeneity of the bare soil infiltration
%
%% parameters for rietkerk model run
function [param,vp,p_noise,sa_plot,nkc]  = rk_2d_heterogeneity_setup_complete(aniso)

	param   = struct();
	vp = struct();

	param.L       = 1024*[1, 1];
	%param.T       = 5e5;
	param.T       = 5e4; % 5e5
	param.opt.dto = inf;
	param.opt.rms_delta_zo_rel_max = 0.2;
	dx            = [1, 1];
	% the initial time step is limited by the reaction term, as the initial
	% contitions is not smooth, the time step is adapted to about 0.5
	% shortly after the initial state has diffused
	param.nx      = (param.L./dx);
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_spectral';

	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
	param.opt.dt_min       = 1/400;
	param.opt.dt_max       = 0.5;
	param.opt.outer_abstol = 1e-4;
	param.opt.outer_reltol = 1e-2;
	param.opt.dt_max_scale_up   = sqrt(2);
	param.opt.dt_min_scale_down = 0;
	% param.opt.rms_delta_zo_rel_max=0.2;
	% param.opt.dto = inf;

	param.pmu.a  = 0.2;
	% coefficient of variation of a at unit grid cell size
	%s_a = (0.00:0.01:0.50);

	% standard deviation of a at unit grid cell size
	%vp.sa_times_a = param.pmu.a*s_a;

	% seed of random number generator
	vp.seed = 0; % (0:2);

	if (aniso)
		% sa for which spectra is plotted
		sa_plot = [0,0.01,0.10,0.13];
	else
		sa_plot = [0,0.3,0.44];
	end

	param.opt.loadfinal = true;
        param.initial_condition = 'obj.random_state(4,[],[],8)';
	param.opt.output_class = @half;
	param.opt.compute_class = @single;
	param.boundary_condition = {'circular','circular'};
	if (aniso)
		param.pmu.vx = [0; 0; 10];
		param.pmu.vy = [0; 0;  0];
		param.pmu.ex = [0.1,0.1,10];
		vp.pmu.ey = {... [0.1,0.1,100];
                                      [0.1,0.1, 20];
			};
		vp.pmu.R  = 0.9;
		nkc = 2.5;
		%param.psl.a = 0.35;
		%param.psdist.a = 'geometric-pink';
		%param.psl.a = 2;
	else
		param.pmu.vx  = [0,0,0];
		param.pmu.vy  = [0,0,0];
		param.pmu.ex  = [0.1,0.1,100];
		vp.pmu.ey        = {[0.1,0.1,100]};
		%vp.pmu.R      = [0.7,0.85,1];
		% local 1.15
		hostname = char(java.net.InetAddress.getLocalHost.getHostName);
		switch (hostname)
		case {'balmung'}
			% even
			%s_a = (0.00:0.02:0.50);
			%vp.pmu.R      = [0.8,1.15];
			%vp.pmu.R = 1.0;
			%s_a = [ 0    0.01   0.02    0.03    0.08    0.09    0.42];
			vp.pmu.R = 0.8;
			s_a = [0.09 0.14 0.26 0.32 0.43 0.48];
		case {'joyeuse'}
			% odd
			%s_a = (0.01:0.02:0.50);
			%vp.pmu.R      = [1.00,1.15];
			%vp.pmu.R = 1.15;
			%s_a = [0.05 0.06 0.19 0.25 0.47];
			%vp.pmu.R = 1;
			%s_a = [0.01   0.08];
			vp.pmu.R = 1.15;
			s_a = [0.19];
		otherwise
			%vp.pmu.R      = 1.15;
		end
		nkc           = 3;
	end
	vp.sa_times_a = param.pmu.a*s_a;

	p_noise = 2;
	if (true)
		param.psdist.a = 'geometric-ornstein-uhlenbeck';
		vp.psl.a       = [min(256,param.L(1)/4)];
	else
		param.psdist.a = 'geometric-pink';
		param.psl.a = 2;
	end

end

