% Tue 18 Feb 16:08:47 CET 2025
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% magnitude of perturbation
s = [0,1];
% number of tiles per row and column the patterns is split for analuysis
ntile = 2;


% number of pixels along unit tile
n1 = 8;
% side length of unit tile prototype
L1 = n1/2;
% oversampling factor of tiles
sdx = 8;
% number of tiles per row and column in the pattern
sL  = 8;
% radius of smoothing kernel for density estimates in pixels
nf = 0.5*sqrt(sL*n1);


ns = length(s);
nx = n1*sdx*sL*[1,1];

% axes in frequency domain
[fx,fy,frr] = fourier_axis_2d(L1*sL*[1,1],nx);
df=1./(sL*L1);

% plot limit
lim_ = 2;

aa = zeros(nx(1)/ntile,nx(2)/ntile,ntile,ntile);
bb = zeros(nx(1)/ntile,nx(2)/ntile,ntile,ntile);
barSS = zeros(nx(1)/ntile,nx(2)/ntile,ntile,ntile);
frSrc = zeros(ns,1);
frc = zeros(ns,1);
fcrr  = zeros(ns,1);

for idx=1:ns
% reset random number generator for reproducibility
rng(0)

% generate unit tile
z1 = meander_tile(n1);
% repeat unit tile to pattern and perturb
[b,x,ex,ey] = tesselate(z1,L1,sdx,sL,s(idx));

% Fourier transform
fb = (fft2(b-mean(b(:))));
% periodogram
hatS = abs(fb).^2;
% normalize periodogram
hatS = hatS/(sum(hatS,'all')*df^2);
hatS_ = ifftshift(gaussfilt2(fftshift(hatS),1));
% density estimage
barS  = ifftshift(gaussfilt2(fftshift(hatS),nf));
%barfb = ifftshift(gaussfilt2(fftshift(fb),nf));
% radial density
[Sr,fr_] = periodogram_radial(hatS,sL*L1*[1,1]); %,sdx*n1*[1,1]);
% characteristic frequency
[frSrc(idx,1),mdx] = max(cvec(fr_).*cvec(Sr.normalized));
frc(idx,1) = fr_(mdx); 
[Sc(idx,1),mdx] = max(barS(:));
fcrr(idx,1) = frr(mdx);
% noise spectrum
Sp = abs(fft2(ex-mean(ex(:)))).^2;
Spr = periodogram_radial(Sp,sL*L1*[1,1]);

% plot pattern
splitfigure([ns,5],[1,1+5*(idx-1)],fflag);
imagesc(frc(idx)*x,frc(idx)*x,b)
axis equal
axis tight
xlabel('Position x/\lambda_c')
ylabel('Position y/\lambda_c')
colormap(flipud(gray))

% plot periodogram
splitfigure([ns,5],[1,2+5*(idx-1)],fflag);
fx=fourier_axis(x);
imagesc(fftshift(fx)/frc(idx),fftshift(fx)/frc(idx),fftshift(hatS_));
axis(lim_*[-1,1,-1,1]);
axis square

% plot spectral density
splitfigure([ns,5],[1,4+5*(idx-1)],fflag);
imagesc(fftshift(fx)/frc(idx),fftshift(fx)/frc(idx),fftshift(barS));
axis(lim_*[-1,1,-1,1]);
axis square

% plot noise spectrum
splitfigure([ns,5],[1,5+5*(idx-1)],fflag);
imagesc(fftshift(fx/frc(idx)),fftshift(fx/frc(idx)),fftshift(Sp));
axis(2*[-1,1,-1,1]);
colormap(flipud(gray))
axis square
xlabel('Wavenumer k_x/k_c');
ylabel('Wavenumer k_y/k_c');

% split in quadrants and analyze separately
n   = size(b,1);
fx_ = fourier_axis(sL*L1/2,n(1)/2);
m   = ntile;
% row
for i1=1:m
	% column
	for i2=1:m
		i1_ = (1:n/m)+(i1-1)*n/m;
		i2_ = (1:n/m)+(i2-1)*n/m;
		% pattern of current quadrant
		b_ = b(i1_,i2_);
		b_ = b_-mean(b_,'all');
		bb(:,:,i1,i2) = b_;
		% window for smoothing transition at quadrant boundary
		w = tukeywin(size(b_,1));
		ww = cvec(w)*rvec(w);
		% fourier transform
		fb_ = fft2(ww.*b_);
		%bf__(:,:,i1,i2) = fb_;
		% angle
		a = angle(fb_);
		s_ = sin(a);
		c_ = cos(a);
		% smooth angle
		s_ = ifftshift(gaussfilt2(fftshift(s_),nf/sqrt(2)));
		c_ = ifftshift(gaussfilt2(fftshift(c_),nf/sqrt(2)));
		a_ = atan2(s_,c_);
		%ex_(i1,i2) = mean(ex(i1_,i2_),'all');
		%ey_(i1,i2) = mean(ey(i1_,i2_),'all');
		aa(:,:,i1,i2) = a_;
		% periodogram (not normmalized here)
		hatS_ = abs(fb_).^2;
		barS_ = ifftshift(gaussfilt2(fftshift(hatS_),nf/sqrt(2)));
		% density
		barSS(:,:,i1,i2) = barS_;

		% plot pattern
		splitfigure([m,m],[10*idx,(i1-1)*m+i2],fflag);
		imagesc(x(i1_),x(i1_),ww.*(b_));
		axis equal
		axis tight

		% plot angle
		splitfigure([m,m],[20*idx,(i1-1)*m+i2],fflag);
		imagesc(fftshift(fx_),fftshift(fx_),fftshift(a_));
		axis(lim_*[-1,1,-1,1]);
		axis square
		colormap(flipud(gray))		
		xlabel('Wavenumer k_x/k_c');
		ylabel('Wavenumer k_y/k_c');

		% plot density
		splitfigure([m,m],[30*idx,(i1-1)*m+i2],fflag);
		imagesc(fftshift(fx_)/frc(idx),fftshift(fx_)/frc(idx),fftshift(barS_));
		axis(lim_*[-1,1,-1,1]);
		axis square
		xlabel('Wavenumer k_x/k_c');
		ylabel('Wavenumer k_y/k_c');
		colormap(flipud(gray));	
	end
end

end

% phase deviation with increasing distance
e = zeros(length(x),1);
for idx=1:length(x)
	e(idx,1) = mean((ex-circshift(ex,idx-1,1)).^2 + (ey-circshift(ey,idx-1,1)).^2,'all');
end
figure(1e3);
plot(x,e)
xlabel('x');
ylabel('relative shift');

% correlation between tiles
C_b = zeros(ntile^2);
C_a = zeros(ntile^2);
C_S = zeros(ntile^2);
for idx=1:ntile^2
	for jdx=1:ntile^2
		C_b(idx,jdx)  = corr(flat((bb(:,:,idx))),flat((bb(:,:,jdx))));
		C_a(idx,jdx)  = corr_angle(flat(aa(:,:,idx)),flat(aa(:,:,jdx)));
		C_S(idx,jdx)  = corr(flat((barSS(:,:,idx))),flat((barSS(:,:,jdx))));
	end
end
fdx=tril(true(ntile^2),-1);
r2_b = C_b.^2
r2_b = mean(r2_b(fdx))
r2_S = C_S.^2
r2_S = mean(r2_S(fdx))
r2_a = C_a.^2
r2_a = mean(r2_a(fdx))

if (pflag)
	ps = 2.5;
	pdfprint(11,'img/tiled-meander-pattern.pdf',ps);
	pdfprint(15,'img/tiled-meander-pattern-phase-noise-spectrum.pdf',ps);
	ps = 3.5;
	for idx=1:4
		pdfprint(200+idx,sprintf('img/tiled-meander-phase-quadrant-%d.pdf',idx),ps);
		pdfprint(300+idx,sprintf('img/tiled-meander-density-quadrant-%d.pdf',idx),ps);
	end
end

