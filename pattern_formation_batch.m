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
%%	dependencies_determine('dependencies.csv','mat/profile-3.mat',{'pattern_formation_patch','pdfprint'})

	% set to true to save fitures to files 
	pflag = false;

	meta = pattern_formation_metadata();
	meta.pflag = pflag;

	mkdir('mat/');
	mkdir('mat/may/');
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
		if (~license('test',toolbox_C{idx,1}))
			printf('%s is missing, execution will likely fail at a later point.\',toolbox_C{idx,2});
		end
	end
	url  = 'https://raw.githubusercontent.com/karlkastner/auxiliar/master/dependencies_fetch.m';
	dest = './lib/auxiliar/dependencies_fetch.m';

	urlwrite(url,dest);

	% this line needs only to be run when packing the code
	% dependencies_determine(meta.filename.dependencies,meta.filename.profile,{'pattern_analysis_batch','pdfprint'});

	% fetch library files
	% this line needs only to be run when the source code is downloaded
	% from the repository without dependencies
	% dependencies_fetch(meta.url,meta.filename.dependencies);

	addpath_recursive('lib/');

	minimum_working_example();

	% Figure 01 : natural regular patterns from aerial images
	close all;
	plot_pattern_observed_regular(meta);

	% Figure 01 : natural regular patterns from aerial images
	close all;
	plot_pattern_observed_irregular(meta);

	% Figure 02 : schematic filtering
	close all;
	plot_filter_schematic_2d(meta);

	% Figure 03 : fraction of ground covered by vegetation vs precipitation
	close all;
	rk_experiment_aridity_transition();
	
	% Figure  04 : heterogeneity map, spectrum and distribution
	close all;
	plot_exogenous_heterogeneity();

	% Figure 05 : irregular model
	% Figure SI 3, 4
	close all;
	grazing_model_experiment();

	% Figure 06 : rietkerk model generated patterns regular isotropic
	close all;
	rk_2d_heterogeneity_experiment(meta,0);

	% Figure 07a-d : pattern properties vs exogenous heterogeneity
	close all;
	clear tab
	rk_2d_heterogeneity_series_postprocess_isotropic();
	
	% Figure 07e : spectral coherence
	close all;
	rk_2d_heterogeneity_series_plot_isotropic_spectral_coherence();

	% figure 08: bandpass-like frequency response of the isotropci RK-model
	close all;
	rk_1d_frequency_response(meta);

	% figure 09 : bandpass bandpass generated patterns and density
	close all;
	plot_bandpass_2d(meta);

	% Figure 10 : rietkerk model generated patterns regular anisotropic
	rk_2d_heterogeneity_experiment(meta,1);

	% Figure 11 : ansisotropic pattern properties
	close all;
	clear tab
	rk_2d_heterogeneity_series_postprocess_anisotropic();

	% Figure 12a : phase-noise-integrating property of the anisotropic RK-model
	close all;
	rk_1d_experiment_phase_shift_single_bump(meta);

	% Figure 12b : oscillator with phase noise through integration of heterogeneity
	close all
	rk_1d_experiment_phase_noise_integration(meta);

	% Figure 13 : phase-noise-integration patterns and density
	close all;
	plot_noisy_oscillator_2d(meta);

	% Figure 14 : schematic pattern formation through filtering
	close all
	plot_pattern_formation_schematic(meta);

	% Figure SI1 - SI2 heterogeneity modelling
	close all
	experiment_heterogeneity_artefacts_2();
	close all
	experiment_heterogeneity_artefacts();
	close all
	experiment_heterogeneity_correlation_length_finite_domain_size();
	close all
	experiment_heterogeneity_correlation_finite_spatial_resolution();

	% Figure SI 6 : perturbation of periodic pattern by phase nois
	close all
	example_phase_noise_meander_tiles();


