% Mon 31 May 20:20:46 CEST 2021
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
%% metadata for analysis and model runs of environmental spatial patterns
%
function meta = pattern_formation_metadata()
	
	meta.url          = 'https://github.com/karlkastner/';

	meta.filename.observed_patterns    = 'mat/observed-patterns-formation.mat';
	meta.filename.dependencies = 'dependencies.csv';
	meta.reload = true;	

	%	
	% parameters for example rietkerk model run
	% and corresponding density plots
	%

	% TODO most of these values have moved to the scripts
	% and are superfluous	

	% transect length (m)
	meta.example.L  = 8e3;
	% grid interval (m)
	meta.example.dx = 2;
	% central wrequency = 1/wavelength (1/m)
	meta.example.fc = 0.01;
	% regularity
	meta.example.Sc_fc = 0.67;
	% flow velocity (~ sqrt(hill slope))
	meta.example.vh = 10;
	% spatial heterogeneity
	meta.example.sd_a = 0.11;

	%
	% parameters for rietkerk series run
	%

	meta.mod2d.L  = 500;	
	meta.mod2d.dx = 2;	
	meta.mod2d.Ti = 1e3;
	meta.mod2d.To = 1e4; % 5e4
	meta.mod2d.dt = 1; % ?
	meta.mod2d.dto = 10;
	meta.mod2d.eh = 100;
	meta.mod2d.vh = 0;

	% transect length
	meta.mopt.L    = 1.6e4;
	% grid interval (m)
	meta.mopt.dx   = 1;
	% duration of initial run
	meta.mopt.Ti   = 1e4;
	% duration of second run
	meta.mopt.To   = 2e4;
	% output interval
	meta.mopt.dto  = 50;
	meta.mopt.dt   = 0.5;
	% spatial heterogeneity
	meta.mopt.sd_a = 0.00:0.01:0.20;
	% surface water diffusivity
	meta.mopt.eh   = 1;
	% surface water velocity
	meta.mopt.vh   = 10; %1:1:20;
	meta.mopt.solver    = @Rietkerk.solve_split;

	% recommended geomean or median as these yield consistent results for wavelength and wavenumber(frequency),
	% 	i.e. geomean(wavelength) = 1/geomean(f_c)
	% (arithmetic)-mean and harmonic-mean do not yield consistent results:
	%	i.e. mean(wavelength) != 1/mean(f_c)
	meta.mfun  = @median;

	% plot options
	meta.pflag = false;
	meta.dflag = false;
	meta.analyze = false;
	meta.visible = true;
	meta.colororder = [0,0,0; 0.9,0,0; 0,0.25,0.8];
	meta.areacol  = [0.55,0.75,1];
	meta.aspect   = 4/3;
	meta.plotscale = 4;
	meta.pattern.xlim = [0,3.5]; % 4.5
	meta.pattern.xlabel = '$x / \lambda_c$';
	meta.pattern.ylabel = '';
	meta.periodogram.xlim = [0,3.5]; % 4.5
	meta.periodogram.ylim = [0,2.8];
	meta.periodogram.ytick = (0:4);
	meta.acf.ytick = -0.5:0.5:1;
	meta.pattern.ytick = 0:4;
	meta.periodogram.xlabel = '$k / k_c$';
	meta.periodogram.ylabel = ''; %'$S \frac{k_\mu}{2 \pi}\;\;\;\;$',
	meta.acf.xlim = [0,1.5]; % 2.5
	meta.acf.ylim = [-0.5,1.05];
	meta.acf.xlabel = '$x / \lambda_c$';
	meta.acf.ylabel = '';

	% colormap for line plots
	meta.colormap = [0,0,0;
		       0.8,0,0;
		       0,0.2,0.8];
	%meta.colormap_b = flipud(colormap_vegetation(256));
	meta.colormap_b = flipud(gray);
	q = 2/3;
	meta.fcmap = @(n) flipud(q*colormap(gray(n)) + (1-q)*colormap_vegetation(n));
end

