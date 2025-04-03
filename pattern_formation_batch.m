% 2021-11-26 19:40:22.060116517 +0100
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
%% batch scripts for generating figures in the associated manuscript
%%
%%
%% note : generation of the 2d-patterns can take up to two weeks,
%%       the run time can be cut by reducing the number of cases and resolution
%%	 in rk_2d_heterogeneity_experiments
%%
%%	dependencies_determine('dependencies.csv','mat/profile-3.mat',{'mimimum_working_example'})

	% set to true to save fitures to files 
	pflag = false;

	meta = pattern_formation_metadata();
	meta.pflag = pflag;

	mkdir('mat/');
	mkdir('img/');
	mkdir('lib/');
	mkdir('lib/auxiliar');

	% Toolbox check
	toolbox_C = {
		'image_toolbox',           'Image Processing Toolbox'
		'signal_toolbox',          'Signal Processing Toolbox'
		'statistics_toolbox',      'Statistics and Machine Learning Toolbox'
		... % 'symbolic_toolbox',        'Symbolic Math Toolbox'
	};

	for idx=1:size(toolbox_C,1)
		if (~license('test',toolbox_C{idx,1})
			printf('%s is missing, execution will likely fail at a later point.\',toolbox_C{idx,2});
		end
	end
	url  = 'https://raw.githubusercontent.com/karlkastner/auxiliar/master/dependencies_fetch.m';
	dest = './lib/auxiliar/dependencies_fetch.m';
	urlwrite(url,dest);

	% fetch library files
	% dependencies_determine(meta.filename.dependencies,meta.filename.profile,{'pattern_analysis_batch','pdfprint'});
	dependencies_fetch(meta.url,meta.filename.dependencies);

	addpath_recursive('lib/');

	minimum_working_example();

	% Figure 01 : natural regular patterns from aerial images
	close all;
	plot_patterns_observed_regular(meta);

	% Figure 01 : natural regular patterns from aerial images
	close all;
	plot_patterns_observed_irregular(meta);

	% Figure 02
	close all;
	plot_filter_schematic_2d(pflag);

	% Figure
	close all;
	plot_exogenous_heterogeneity();

	% Figure 03 (thresholding)
	experiment_rk_aridity_transition();
	
	% Figure 03
	grazing_model_experiment();

	% Figure 04, 07 : computer generated patterns, regularity, corr(a,b)
	close all;
	% TODO rename exp in series simulate 
	rk_2d_heterogeneity_experiment(meta);
	close all;
	rk_2d_series_postprocess_aniso();
	close all;
	rk_2d_heterogeneity_postprocess_iso();
	close all;
	rk_2d_hetero_series_plot_isotropic_spectral_coherence();

	% figure 05: bandpass-like frequency response of the isotropci RK-model
	close all;
	rk_1d_frequency_response(pflag);

	% figure 06 : bandpass bandpass generated patterns and density
	close all;
	plot_bandpass_2d(meta);

	% Figure 08 : phase-noise-integrating property of the anisotropic RK-model
	close all;
	rk_1d_phase_shift_experiment_single_bump(meta);

	% Figure 08
	close all
	rk_1d_experiment_phase_noise_integration(meta);

	% Figure 09
	close all
	%phase_noise_illustration();
	example_phase_noise();

	% Figure 10 : phase-noise-integration patterns and density
	close all;
	plot_noisy_oscillator_2d(meta);

	% Figure 11
	close all
	plot_pattern_formation_schematic(meta);

	% supplement
	close all
	experiment_ou_correlation_length()

	% supplement
	close all
	experiment_heterogeneity_variance_vs_spatial_resolution()
	example_exogenous_heterogeneity_artefacts();
	plot_heterogeneity_distribution();

	% supplement
	close all
	experiment_heterogeneity_artefacts()

	% supplement
	close all
	example_phase_noise_meander_tiles();

