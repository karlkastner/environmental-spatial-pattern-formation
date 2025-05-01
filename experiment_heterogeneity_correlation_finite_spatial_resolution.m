% 2025-03-04 18:07:14.580252644 +0100

%
% demonstrate the effect of a finite spatial resolution on modelling the
% two dimensional spatial ornstein uhlenbeck process and reducing the artefacts
% by oversampling:
% - underestimation of variance 
% - underestimation of the correlation length
%

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% spatial extent
L     = 256*[1,1];
% spatial resolution
dx    = [1,2,4,8,16];

% mean of the spatial heterogeneity
mu_z = 1;
% coefficient of determination of the spatial heterogeneity
cv_z = 0.1;
% standard deviation of the spatial heterogeneity
sd_z = cv_z*mu_z;
% correlation length of the spatial heterogeneity
theta_z = theta;
% oversampling factor in the spatial domain
n_spatial = [1,3,5,7];
% oversampling factor in the spectral domain
m_spectral = 3;
% window in the spatial domain
pw = 1;

% declare variables
sd_dx = [];
theta_dx = [];
for idx=1:length(dx)
	disp(idx);
	for jdx=1:length(m_spectral)
		n  = L/dx(idx);
		% spatial axis
		x  = (0:n(1)-1)'*L(1)/n(1);
		% spectral axis
		fx = fourier_axis(n(1),L(2));
		% generate heterogeneity
		[z, C, S] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z,L,n,m_spectral(jdx),m_spectral,pw);
		S = real(S);
		C = real(C);
		% effection standard deviation
		sd_dx(idx,jdx) = sqrt(C(1,1));
		% radial correlation
	 	Rx = C(1:end/2,1)/C(1,1);
		% values are constant zeros for high lags and small theta
		Rx = make_monotonic(Rx,-1);
		% effecive correlation length
		theta_dx(idx,jdx) = interp1(Rx,x(1:end/2),exp(-1),'linear');
	end % for jdx
end % for idx

% 4 -> 4
splitfigure([2,2],[1,1],fflag);
semilogx(dx/theta,sd_dx/sd,'.-');
ylabel('\sigma_{\Delta{x}}/\sigma_0');
xlabel('\Delta{x}/\theta_0');
lh=legend(num2str(cvec(m_spectral)),'location','southwest');
title(lh,'n_i')
axis square
ylim([0.4,1.1]);

splitfigure([2,2],[1,2],fflag);
semilogx(dx/theta,theta_dx/theta,'.-');
xlabel('\Delta{x}/{\theta_0}');
ylabel('$\theta_{\Delta{x}} / \theta_0$','interpreter','latex');
axis square

if (pflag)
	ps = 3;
	pdfprint(11,'img/random-field-sigma-vs-dx',ps);
	pdfprint(12,'img/random-field-theta-vs-dx',ps);
end

