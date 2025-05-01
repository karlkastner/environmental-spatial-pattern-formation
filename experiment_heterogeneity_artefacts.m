% 2025-04-02 10:29:37.827212507 +0200

%
% effect of a finite domain size when modelling the two dimensional spatial
% ornstein uhlenbeck process
% and techniques for minimizing artefacts:
% - windowing
% - oversampling
%

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% mean of the heterogeneity
mu_z=0.2;
% coefficient of variation of the heterogeneity
cv_z = 0.1;
% standard deviation of the heterogeneity
sd_z = cv_z*mu_z;

% spatial extent
L  = 256*[1,1];
% spatial resolution
dx = 1;
% number of grid points
nx=L/dx;
% relative size of fourier window
pw = [0,1,0,1];
% oversampling factor in the spectral domain 
m_spectral = [1,1,3,3];
% oversampling factor in the spatial domain
m_spatial = 5*[1,1,1,1]; % 11
theta_a = L(1)/4*[1,1,1,1];

% allocate memory
S_ = [];
Sr_ = [];
S_ = [];
% for each case
for idx=1:length(theta_a)
	disp(idx);
	% reset random number generator
	rng(0)
	% axes in the frequency domain
	fx = fourier_axis(nx(1),L(2));
	% generate heterogeneity
	[z, C, S] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_a(idx),L,nx,m_spatial(idx),m_spectral(idx),pw(idx));
	S = real(S);
	% radiam periodogram
	[Sr,fr]=periodogram_radial(real(S));

	splitfigure([2,2],[1,1],fflag);
	x = (0:nx(1)-1);
	plot(x,C(:,1)/C(1,1));
	hold on
	xlim([0,L(1)/2]);
	xlabel('Distance x');
	ylabel('Autocorrelation R');
	legend('plain','windowed','oversampled','combined');
	axis square

	splitfigure([2,2],[1,2],fflag);
	plot(fr,Sr.normalized)
	hold on

	splitfigure([2,2],[1,3],fflag);
	S__ = max(eps,S(1:end/2,1));
	loglog(fx(1:end/2),S__,'.-');
	hold on;
	ylim([1e-6,1]);
	ylabel('Density S(k_x,0)');
	xlabel('Frequency k_x/(2 \pi)');
	xlim([fx(1),max(fx)]);
	axis square;
	
	splitfigure([2,2],[1,4],fflag);
	plot(fr,1./Sr.normalized);
	hold on;
	Sr_(:,idx) = Sr.normalized;
	S_(:,idx) =real(S(:,1));
end

% goodness of fit
r2 = 1-rms(Sr_-Sr_(:,end)).^2./var(Sr_(:,end));
erel = max(S_-S_(:,end))./max(S_(:,end));
%max(Sr_-Sr_(:,end))./max(Sr_(:,end))

printf('r2 %g\n',r2);
printf('erel %g\n',erel);

if (pflag)
	ps = 2.5;
	pdfprint(11,'img/random-field-autocorrelation.pdf',ps);	
	pdfprint(13,'img/random-field-spectral-density.pdf',ps);	
end

