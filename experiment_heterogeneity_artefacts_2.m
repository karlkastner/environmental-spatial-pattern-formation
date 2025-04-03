% 2025-02-25 15:46:37.215456670 +0100

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

n = 1024*[1,1];
dx = 1;
L  = dx*n;


mu = 1;
sd = 0.1;
%[lmu,lsd] = lognpdf_moment2par(mu,sd)
%pause

[x,y] = fourier_axis_2d(n./L,n);
fx = fourier_axis(L(1),n(1));

theta_ = [1,L(1)/2];
ni_    = [1,11];
o_     = [1,1];
pw     = [0,0.5];
vv = [];
for tdx = 1:length(theta_)
theta = theta_(tdx);

for jdx=1:length(o_)

rng(0);

o = o_(jdx);
nint = ni_(jdx);
% oL = o*L;
% on = o*n;

hatS = 0;
n_ = 1;
for idx=1:n_
tic
 [v, C, S, lz, lC, lS,w] = geometric_ou_2d_grid_cell_averaged_generate(mu,sd,theta,L,n,nint,o,pw(jdx));
%v = real(v);
%C = real(C);
%C = max(C,0);
toc
 hatS = hatS + abs(fft2(v-mean(v(:)))).^2;
end
vv(:,:,jdx) = v;
if (0)
hatS = hatS/n_;
hatR = real(ifft2(hatS));
hatR=hatR/hatR(1);
[hatSr, fr] = periodogram_radial(hatS,L);
end
[Sr,fr] = periodogram_radial(S,L);


splitfigure([2,5],[tdx,1],fflag);
if (1==jdx) cla; end
plot(x,C(1,:)/sd^2,'.');
hold on
xlim([-0.5,min(5*theta,L(1)/2)+0.01]);
ylim([0,1.05])
axis square
xlabel('Lag distance r');
ylabel('Radial auto-covariance $C_{\hat z} / \sigma_z^2$','interpreter','latex');
if (1==tdx)
	legend('w artefacts','w/o artefacts');
end

splitfigure([2,5],[tdx,5+(jdx-1)*5],fflag);
if (1==jdx) cla; end
%plot([C(:,1)/C(1,1),lC(:,1)/lC(1,1)])
%imagesc(w)
%axis equal
%axis tight

splitfigure([2,5],[tdx,6],fflag);
if (1==jdx) cla; end
loglog(fr,Sr.normalized/Sr.normalized(2),'.-')
hold on
axis square
xlabel('Frequency k/(2 \pi)');
ylabel('Radial density $S_{\hat z}$','interpreter','latex');

splitfigure([2,5],[tdx,2+(jdx-1)*5],fflag);
if (1==jdx) cla; end
imagesc((v-mu)/sd)
axis equal
axis tight
ch=colorbar
title(ch,'$(\hat z - \mu_z)/\sigma_z$','interpreter','latex');
axis(min(5*theta,L(1)/2)*[0,1,0,1]+0.5);
caxis([-1,1]*4);
xlabel('Position $x$','interpreter','latex');
ylabel('Position $y$');
axis square

splitfigure([2,5],[tdx,3+(jdx-1)*5],fflag);
if (1==jdx) cla; end
imagesc(fftshift(log(C/C(1))))
axis equal
axis tight
colorbar
title('C')

splitfigure([2,5],[tdx,4+(jdx-1)*5],fflag);
if (1==jdx) cla; end
imagesc(real(log(fftshift(S))))
axis equal
axis tight
colorbar
title('S')

mu_(tdx,jdx) = mean(v,'all')
sd_(tdx,jdx) = std(v,[],'all')
cv_ = sd_./mu_

figure(3)
subplot(2,3,jdx)
histogram(flat(vv(:,:,jdx)))

end % for jdx
r2(tdx) = 1 - (rms(vv(:,:,1)-vv(:,:,2),'all')/sd)^2

if (pflag)
	ps = 3.5;
	base = sprintf('img/heterogeneity-theta-%d',theta_(tdx));
	pdfprint(tdx*10+1,[base,'-autocorrelation.pdf'],ps);
	pdfprint(tdx*10+6,[base,'-density.pdf'],ps);
	pdfprint(tdx*10+2,[base,'-with-artefacts.pdf'],ps);
	pdfprint(tdx*10+7,[base,'-without-artefacts.pdf'],ps);
end

end % for tdx
