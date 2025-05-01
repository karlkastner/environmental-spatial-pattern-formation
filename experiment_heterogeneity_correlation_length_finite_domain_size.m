% Tue  4 Mar 11:56:22 CET 2025
%
% demonstrate the effect of a finite domain size on the two dimensional
% spatial ornstein uhlenbeck process with long correlation length and
% techniques reducing the artefacts such as
% - ripples in the tail of the spectral density
% - deviation of the correlation length
%

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;


% spatial extent
L  = 128*[1,1];
% spatial resoltion
dx = 1;
% relative width of window
pw = [0,1,0,1];
% oversampling factor in the spectral domain
m_spectral = [1,1,3,3];
% oversampling factor in the spatial domain
m_spatial = 5*[1,1,1,1];

% mean of the exogeneous spatial heterogeneity
mu_z=0.2;
% coefficient of variation of the exogenous spatial heterogeneity
cv_z = 0.1;
% standard deviation of the exogenous spatial heterogeneity
sd_z = cv_z*mu_z;
% correlation length of the exogenous spatial heterogeneity
theta_z = 2.^(0:0.25:7);

% number of grid points
n=L/dx;
% spatial axis
x=0:n(1)/2-1;
theta_effective = [];
Sr = [];
Rr = [];
for idx=1:length(theta_z)
	disp(idx);
	R_ = [];
	Sr = [];
	for jdx=1:length(pw)
		% generate heterogeneoity
		[z, Cxy, Sxy] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z(idx),L,n,m_spatial(jdx),m_spectral(jdx),pw(jdx));
		Sxy = real(Sxy);
		Cxy = real(Cxy);
		% radial density
		Sr(:,1) = real(S(:,1));
		%radial autocorrelation
		Rr(:,jdx) = Cxy(1:end/2,1)/Cxy(1,1);
		% values are constant zeros for high lags and small theta
		Rr(:,jdx) = make_monotonic(Rr(:,jdx),-1);
		% effective correlation length
		theta_effective(idx,jdx) = interp1(Rr(:,jdx),x,exp(-1),'linear');

		[d,ddx]      = max(diff(S(1:end/2,1)));
		[di,ddx]     = min(diff(1./S(1:end/2,1)));
		ddx_(jdx)    = ddx;
		d_(idx,jdx)  = d;
		di_(idx,jdx) = di;
	end
	figure(1);
	subplot(5,6,idx)
	cla
	plot(x,Rr);
	hold on
	plot(theta_effective(idx,:),exp(-1),'.');

	figure(3)
	subplot(5,6,idx)
	cla
	plot(Sr(1:end/2));
	hold on
% plot(t(idx,:),exp(-1),'.');
end % for idx

splitfigure([2,2],[4,1],fflag);
semilogx(theta_z/L(1),(theta_effective-theta_effective(:,1))./theta_effective(:,1),'.-','linewidth',1);
xlabel('Correlation length $\theta/L$','interpreter','latex');
ylabel('Deviation of correlation lentgh \Delta \theta / \theta');
xlim([1/L(1),0.5]);
set(gca,'xtick',2.^(-7:0),'xticklabel',rats(2.^(-7:0)'));
legend('plain','windowed','oversampled','combined','location','northwest')
axis square

splitfigure([2,2],[4,2],fflag);
semilogx(theta_z/L(1),max(d_,0),'.-','linewidth',1);
ylabel('Max Ripple in Desity S');
xlabel('Correlation length $\theta/L$','interpreter','latex');
xlim([1/L(1),1]);
set(gca,'xtick',2.^(-7:0),'xticklabel',rats(2.^(-7:0)'));
legend('plain','windowed','oversampled','combined','location','northwest')
axis square

if (pflag)
	ps = 2.5;
	pdfprint(41,'img/random-field-error-in-correlation-length.pdf',ps);
	pdfprint(42,'img/random-field-ripples-in-stop-band.pdf',ps);
end

