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

	meta.pflag = false;

	% plot options for series runs
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
end % pattern_formation_metadata

