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
function [param,vp,p_noise,sap,nkc]  = setup1(aniso)

	param   = struct();
	vp = struct();

	param.L     = 1024*[1, 1];
	param.T     = 5e5;
	param.opt.dto = 1e4;
	dx    = [1, 1];
	nx = param.L/dx;
	param.nx = nx;
	param.opt.dt = 0.5; % was 1

	p_noise = 2; % 2
	%s_a = 0:0.01:0.5;
	param.pmu.a  = 0.2;
	%s_a = 0.2:0.1:0.5;
	% 0.0,0.1, 
	if (aniso)
		s_a = [0.00:0.01:0.50];
	else
		s_a = 0:0.01:0.5;
	end

	vp.sa_times_a = param.pmu.a*s_a;
	vp.seed = 0;

	%s_a = [0,0.3,0.44];
	if (aniso)
		% sa for which spectra is plotted
		%sap  = [0,0.25,0.5];
		sap = [0,0.01,0.10,0.13];
		%sap = 0.01;
		%sap = 0.09;
	else
		sap = [0,0.3,0.44];
		%s_a = sap;
	end
	%sap = s_a;

	param.opt.loadfinal = true;
        param.initial_condition = 'obj.random_state(4,[],[],8)';
	% TODO there seems to be a matlab but that spoils half output
	param.opt.output_class = 'half';
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
		param.psdist.a = 'geometric-pink';
		param.psl.a = 2;
	else
		param.pmu.vx  = [0,0,0];
		param.pmu.vy  = [0,0,0];
		param.pmu.ex  = [0.1,0.1,100];
		vp.pmu.ey        = {[0.1,0.1,100]};
		vp.pmu.R         = 0.7;
		nkc           = 3;
		if (p_noise == 1)
			param.psdist.a = 'geometric-ornstein-uhlenbeck';
			param.psl.a   = 0.35;
		else
			param.psdist.a = 'geometric-pink';
			param.psl.a = 2;
		end
	end
	param.opt.solver = 'solve_split';
	%param.opt.rng = 0;

end

