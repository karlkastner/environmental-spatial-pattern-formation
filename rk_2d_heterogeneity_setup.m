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
function [param,vp,p_noise,cva_plot,nkc]  = rk_2d_heterogeneity_setup(aniso)

	param   = struct();
	vp = struct();

	% spatial extent
	param.L       = 1024*[1, 1];
	% final time
	param.T       = 5e4; % 5e5
	% output after fixed time step
	param.opt.dto = inf;
	% output after relative change of state
	param.opt.rms_delta_zo_rel_max = 0.2;
	% patial resolution
	dx            = [1, 1];
	% number of grid points
	param.nx      = (param.L./dx);
	%
	% mean of infiltration coefficient
	param.pmu.a  = 0.2;
	% coefficient of variation of a at unit grid cell size
	vp.cva = (0.00:0.01:0.50);
	% process generating exogenous heterogeneity
	param.psdist.a = 'geometric-ornstein-uhlenbeck';
	% correlation length of exogenous heterogeneity
	vp.psl.a       = 512;

	% seed of random number generator
	vp.seed = 0; % (0:2);

	% solver settings
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_spectral';

	% time step settings
	% the initial time step is limited by the reaction term, as the initial
	% contitions is not smooth, the time step is adapted to about 0.5
	% shortly after the initial state has diffused
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
	param.opt.dt_min       = 1/400;
	param.opt.dt_max       = 0.5;
	param.opt.outer_abstol = 1e-4;
	param.opt.outer_reltol = 1e-2;
	param.opt.dt_max_scale_up   = sqrt(2);
	param.opt.dt_min_scale_down = 0;

	if (aniso)
		% sa for which spectra is plotted
		cva_plot = [0,0.01,0.10,0.13];
	else
		cva_plot = [0,0.3,0.44];
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
		vp.pmu.ey = {
                             [0.1,0.1, 20];
			};
		vp.pmu.R  = 0.9;
		nkc = 2.5;
	else
		param.pmu.vx  = [0,0,0];
		param.pmu.vy  = [0,0,0];
		param.pmu.ex  = [0.1,0.1,100];
		vp.pmu.ey     = {[0.1,0.1,100]};
		vp.pmu.R      = [0.8, 1, 1.15];
		if (0)
		% for parallel server runs
		hostname = char(java.net.InetAddress.getLocalHost.getHostName);
		switch (hostname)
		case {'balmung'}
			vp.cva = [0.00:0.02:0.50];
		case {'joyeuse'}
			vp.cva = [0.01:0.02:0.50];
		otherwise
		end
		end
		nkc           = 3;
	end
end % rk_2d_heterogeneity_setup()

