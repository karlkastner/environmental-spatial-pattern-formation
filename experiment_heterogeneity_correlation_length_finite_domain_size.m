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

if (1)

Sr_ = [];
S_ = [];
 L  = 128*[1,1];

 pw = [0,1,0,1];
 no = [1,1,2,2];
 ni = 5*[1,1,1,1];
 
 tt = 2.^(0:0.25:7);

 dx=1;
 mu_z=0.2;
 sd_z = 0.1*mu_z;

 n=L/dx;
x=0:n(1)/2-1;
t = [];
 for idx=1:length(tt)
disp(idx);
C1_ = [];
 for jdx=1:length(pw)
 theta_z = tt(idx);
% tic
 [z, C, S] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z,L,n,ni(jdx),no(jdx),pw(jdx));
S = real(S); 
 S1_(:,1) = real(S(:,1));
 C1 = C(1:end/2,1)/C(1,1);
 % values are constant zeros for high lags and small theta
 C1 = make_monotonic(C1,-1);
 t(idx,jdx) = interp1(C1,x,exp(-1),'linear');

 C1_(:,jdx) = C1;

 [d,ddx]      = max(diff(S(1:end/2,1)));
 [di,ddx]     = min(diff(1./S(1:end/2,1)));
 ddx_(jdx)    = ddx;
 d_(idx,jdx)  = d;
 di_(idx,jdx) = di;
 % / S(ddx+1,1);
 if (4 == jdx)
	%d_(idx,:) = d_(idx,:) ./ rvec(S(ddx_,1))
 end

% toc()
if (0)
% [Sr,fr]=periodogram_radial(real(S));
 subplot(2,2,1)
 plot(C(:,1)/C(1,1));
 hold on
 subplot(2,2,2)
 plot(fr,Sr.normalized)
 hold on
 subplot(2,2,3)
 plot(fr,1./Sr.normalized);
 hold on;
 Sr_(:,idx) = Sr.normalized;
 S_(:,idx) =real(S(:,1));
end


 end
 figure(1);
 subplot(5,6,idx)
 cla
 plot(x,C1_);
 hold on
 plot(t(idx,:),exp(-1),'.');

figure(3)
 subplot(5,6,idx)
 cla
 plot(S1_(1:end/2));
 hold on
% plot(t(idx,:),exp(-1),'.');
end

end

splitfigure([2,2],[4,1],fflag);
semilogx(tt/L(1),(t-t(:,1))./t(:,1),'.-','linewidth',1);
xlabel('Correlation length $\theta/L$','interpreter','latex');
ylabel('Deviation of correlation lentgh \Delta \theta / \theta');
xlim([1/L(1),0.5]);
set(gca,'xtick',2.^(-7:0),'xticklabel',rats(2.^(-7:0)'));
legend('plain','windowed','oversampled','combined','location','northwest')
axis square

splitfigure([2,2],[4,2],fflag);
semilogx(tt/L(1),max(d_,0),'.-','linewidth',1);
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

