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

L     = 256*[1,1];
theta = 4;
dx    = [1,2,4,8,16]
sd    = 0.1;

pw = 0.2;
ni = [1,3,5,7];

mu_z = 1;
sd_z = 0.1*mu_z;
theta_z = theta;

pw = 1;
no = 2;

sd_dx = [];
theta_dx = [];
for idx=1:length(dx)
idx
 for jdx=1:length(ni)
	n  = L/dx(idx);
	x  = (0:n(1)-1)'*L(1)/n(1);
	fx = fourier_axis(n(1),L(2));
	tic();
	[z, C, S] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z,L,n,ni(jdx),no,pw);
	%S = real(S);
	%toc()
	%[Sr,fr]=periodogram_radial(real(S));
	
	sd_dx(idx,jdx) = sqrt(C(1,1));

 	C1 = C(1:end/2,1)/C(1,1);
	% values are constant zeros for high lags and small theta
	C1 = make_monotonic(C1,-1);
	theta_dx(idx,jdx) = interp1(C1,x(1:end/2),exp(-1),'linear');
end
end

% 4 -> 4
splitfigure([2,2],[1,1],fflag);
semilogx(dx/theta,sd_dx/sd,'.-');
ylabel('\sigma_{\Delta{x}}/\sigma_0');
xlabel('\Delta{x}/\theta_0');
lh=legend(num2str(cvec(ni)),'location','southwest');
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

