% Fri  4 Mar 14:09:21 CET 2022
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
%% schematic plot of filter densities, autocorrelation and patterns
%

function plot_filter_schematic(meta)

if (nargin()<1)
	meta = pattern_formation_metadata();
end
pflag = meta.pflag;
fflag = pflag;
ps = meta.plotscale;
fcmap = meta.fcmap;
s = 2.5;
f0 = 6/0.75;

type_C = {'iso','aniso'}

for tdx=1:length(type_C)
	type = type_C{tdx};

	% characteristic wavelength	
	lc = 1;
	% characteristic frequency
	fc = 1/lc;
	% scale factor
	m = 20;
	% spatial extend
	L = m*lc;
	% spatial resolution
	dx = lc/m;
	% number of grid points
	n = round(L/dx);
	p = 10;
	%if (isiso)
	switch (type)
		case {'iso'}
		pt = 0.8;
	otherwise
		pt = (0.5+0.5)/2;
	end
	x = linspace(-L/2,L/2,n)';
	[fx,fy,fr] = fourier_axis_2d(L*[1,1],n*[1,1]);
	axis tight
	
	% density
	% reset random number generator for reproducibility
	rng(0)
	% spectral density
	switch (type)
	case {'iso'}
		% radal density
		S0 = bandpass1d_continuous_pdf(fr,f0/10,1);
	otherwise
		% density perpendicular to stripes 
		%Sx = bandpass1d_continuous_pdf(fx,f0/10,1);
		Sx = phase_drift_pdf(fx,f0/10,0.5);
		% density parallel to stripes
		Sy = phase_drift_parallel_pdf(fx,2.5);
			%Sy = normpdf(fx,0,0.7*f0/10);
		%Sy = exppdf(abs(fx),15/f0);
		% 2D density
		S0 = cvec(Sy)*rvec(Sx);
	end

	% spatial noise	
	e = randn(n);
	
	% density raised to power (modulares regularity)
	S2d = S0.^p;
	% get maximum of the density
	[Sc,mdx]=max(S2d(1,:));
	fc=fx(mdx);
	lc = 1./fc;
	% autocorrelation
	R2d = ifft2(S2d);
	R2d = R2d/R2d(1,1);
	
	% transfer function
	T = sqrt(S2d);
	% spectral noise
	fe = fft2(e);
	% pattern
	b = real(ifft2(T.*fe));
	% periodogram
	Shat = abs(fft2(b)).^2;
	
	splitfigure([2,3],[tdx,1],fflag);
	imagesc(x/lc,x/lc,e);
	axis equal
	axis tight
	axis(s*[-1,1,-1,1]);
	axis off
	%colormap(flipud(gray))
	colormap(fcmap(256));
	caxis(2*[-1,1]);
	
	% acf
	splitfigure([2,3],[tdx,2],fflag);
	imagesc(x/lc,x/lc,ifftshift(R2d));
	axis equal
	axis tight
	axis(s*[-1,1,-1,1]);
	axis off
	colormap(fcmap(256));
	caxis([-1,1]/4)
	
	% pattern
	splitfigure([2,3],[tdx,3],fflag);
	imagesc(x/lc,x/lc,b>quantile(b(:),pt));
	axis equal
	axis tight
	axis off
	colormap(fcmap(256));
	
	splitfigure([2,3],[tdx,4],fflag);
	fehat = abs(fe).^2;
	imagesc(fftshift(fx)./fc,fftshift(fy)/fc,fehat)
	axis equal
	axis tight
	axis(s*[-1,1,-1,1]);
	axis off
	caxis([0,1]*quantile(fehat(:),0.975));
	colormap(fcmap(256));
	
	splitfigure([2,3],[tdx,5],fflag);
	imagesc(fftshift(fx)/fc,fftshift(fy)/fc,ifftshift(S2d));
	axis equal
	axis tight
	axis(s*[-1,1,-1,1]);
	axis off
	%colormap(flipud(gray))
	colormap(fcmap(256));
	
	splitfigure([2,3],[tdx,6],fflag);
	imagesc(fftshift(fx)/fc,fftshift(fy)/fc,fftshift(Shat));
	axis equal
	axis tight
	axis(s*[-1,1,-1,1]);
	axis off
	%colormap(flipud(gray))
	q = 1-0.25*(fc./L)^2
	%q = normcdf(3);
	caxis([0,1]*quantile(Shat(:),q));
	colormap(fcmap(256));
	
	figure(2);
	clf
	%splitfigure(2,2,1)
	df=1/L;S1d = S2d(1,:); S1d = 2*S1d/(sum(S1d)*df); plot(fx/fc,S1d*fc); xlim([0,2.5]);
	
	if (pflag)
			a = 1.3;
			pdfprint(10*tdx+1,['img/filter-2d-',type,'-heterogeneity.pdf'],ps,a);
			pdfprint(10*tdx+2,['img/filter-2d-',type,'-acf.pdf'],ps,a);
			pdfprint(10*tdx+3,['img/filter-2d-',type,'-pattern.pdf'],ps,a);
			pdfprint(10*tdx+4,['img/filter-2d-',type,'-heterogeneity-p.pdf'],ps,a);
			pdfprint(10*tdx+5,['img/filter-2d-',type,'-density.pdf'],ps,a);
			pdfprint(10*tdx+6,['img/filter-2d-',type,'-pattern-p.pdf'],ps,a);
	end
end % for tdx

