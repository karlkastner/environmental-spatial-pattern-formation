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

clf
Sr_ = [];
S_ = [];
if (1)
 L  = 256*[1,1];
 pw = [0,1,0,1];
 no = [1,1,3,3];
 ni = 5*[1,1,1,1]; % 11
 tt = L(1)/4*[1,1,1,1];
% tt = 512*[1,1,1,1]/2;
else
	L  = [128,128];
	ni = 1:11; %[1,3,5,7,9,11];
	pw = ones(length(ni),1); %[1,1,1,1,1,2];
	no = ones(length(ni),1); %[2,2,2,2,2,2];
	tt = ones(length(ni),1);% [1,1,1,1,1,1];
end

dx=1;
mu_z=0.2;
cv_z = 0.1;
sd_z = cv_z*mu_z;

S_ = [];
 for idx=1:length(tt)
disp(idx);
 theta_z = tt(idx);
 n=L/dx;
 fx = fourier_axis(n(1),L(2));
 tic();
rng(0)
 [z, C, S] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z,L,n,ni(idx),no(idx),pw(idx));
 S = real(S);
 toc()
 [Sr,fr]=periodogram_radial(real(S));

 splitfigure([2,2],[1,1],fflag);
 x = (0:n(1)-1);
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

r2 = 1-rms(Sr_-Sr_(:,end)).^2./var(Sr_(:,end))
max(Sr_-Sr_(:,end))./max(Sr_(:,end))
erel = max(S_-S_(:,end))./max(S_(:,end))

if (pflag)
	ps = 2.5;
	pdfprint(11,'img/random-field-autocorrelation.pdf',ps);	
	pdfprint(13,'img/random-field-spectral-density.pdf',ps);	
end

