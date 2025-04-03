% Thu  5 Oct 17:10:16 CEST 2023
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
%% plot a satellite imag of an irregular pattern
%
function g2 = plot_pattern_observed_irregular(meta)
	if (nargin()<1)
		meta = pattern_formation_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	ps = meta.plotscale;
	fcmap = meta.fcmap;
	
	g = GeoImg();
	g.read('img/irregular-mont-ventoux.png');
	n  = 100;
	g2 = g.crop(1,400+n,935-n,1335-400-n);
	% imagesc(g2.img); axis equal
	b = g2.img;
	% coordinate axes
	x = g2.x; x = x-x(1);
	y = g2.y; y = y-y(1);
	L=[range(x),range(y)]
	%s = size(b);
	%x=0:s(1)-1;
	%y=0:s(2)-1;
	% spatial resolution
	dx = x(2)-x(1)
	%L = [x(end),y(end)]
	% spectral axis
	fx = fourier_axis(g2.x);
	fy = fourier_axis(g2.y);
	% spectral resolution
	df = fx(2)-fx(1);

	% pseudo biomass (quick and dirty grayscale conversion)
	b = mean(b,3);
	% periodogram
	Shat = abs(fft2(b-mean(b,'all'))).^2;
	% normalize
	Shat = Shat/(sum(Shat,'all')*df^2);
	% periodogram
	R    = ifft2(Shat);
	% normalize
	R = R/R(1,1);
	%size(b);
	% radial periodogram
	[Sr,fr]=periodogram_radial(Shat,L);
	Sr = Sr.normalized;
	% Sr(1)=0;
	nskip = 2;
	fdxs=1:nskip;
	w = fr;
	w(fdxs) = 0;
	% fit the spectral density of a lowpass-filter
	par=fit_spectral_density(fr,Sr,w,@lowpass2d_continuous_pdf,[1.5,1.5],'ls',0)
	Sr(:,2) = lowpass2d_continuous_pdf(fr,par(1),par(2));
	 fdx=find(Sr(nskip+1:end,2)<0.5*Sr(nskip+1,2),1,'first')+nskip;
	 fl=fr(fdx);
	 Sr=Sr./(sum(Sr)*(fr(2)-fr(1)));
	 Sr(fdxs,1)=NaN;
	
	splitfigure([2,3],[1,1],fflag);
	b = b-min(b,[],'all');
	b = b/max(b,[],'all');
	g = graythresh(b);
	imagesc(x*fl,y*fl,1-(b>g));
	axis equal;
	axis tight;
	xlabel('Position $x/\lambda_l$','interpreter','latex');
	ylabel('Position $y/\lambda_l$','interpreter','latex');
	set(gca,'xtick',0:5:20);
	set(gca,'ytick',0:5:20);
	axis square
	colormap(fcmap(256));
	
	splitfigure([2,3],[1,2],fflag);
	%Shat_ = gaussfilt2(Shat,0.7);
	imagesc(fftshift(fx)/fl,fftshift(fy)/fl,fftshift(Shat)*fl^2);
	%caxis([0,0.125*max(Shat*fl^2,[],'all')]);
	%caxis([0,0.125*max(Shat*fl^2,[],'all')]);
	caxis([0,0.4])
	xlabel('Wavenumber $k_x/k_l$','interpreter','latex');
	ylabel('Wavenumber $k_y/k_l$','interpreter','latex');
	axis square
	axis(3.5*[-1,1,-1,1]);
	cbh=colorbar
	title(cbh,'$\hat S/\lambda_l^2$','interpreter','latex');
	%colormap(flipud(gray(20)))
	colormap(fcmap(20));
	
	splitfigure([2,3],[1,3],fflag);
	cla
	imagesc(fl*(x-x(end)/2),fl*(y-y(end)/2),fftshift(R))
	cbh=colorbar
	title(cbh,'$\hat R$','interpreter','latex');
	xlim(2.5*[-1,1])
	ylim(2.5*[-1,1])
	axis square
	xlabel('Lag distance $x/\lambda_l$','interpreter','latex');
	ylabel('Lag distance $y/\lambda_l$','interpreter','latex');
	caxis([-0.11,1])
	colormap(fcmap(11))
	%	colormap(flipud(gray(10)))
	
	
	splitfigure([2,3],[1,4],fflag);
	plot(fr/fl,Sr*fl,'linewidth',1);
	xlim([0,4]);
	xlabel('Wavenumber $k_r/k_l$','interpreter','latex');
	ylabel('Radial Density $S_r/\lambda_l$','interpreter','latex');
	legend('empirical','LP-fit');
	axis square
	
	splitfigure([2,3],[1,5],fflag);
	cla
	S_ = trifilt2(Shat,7);
	[Sa,theta] = periodogram_angular(S_,L);
	%plot(theta,Sa)   
	plot(theta,Sa,'linewidth',1);
	xlabel('Angle $\theta$');
	ylabel('Density $S_\theta$','interpreter','latex');
	ylim([0,1])
	xlim([-1,1]*pi/2);
	hline(1/pi,'color','red','linewidth',1);
	axis square
	

	f_50 = fl;
	%f_50 = fc;
	dfr  = df(1)
	nf_test = round(0.25*f_50/dfr)
	%fmsk = (frr<4*fc);
	bmsk = [];
	% TODO use Spatial_Pattern/analyse_grid here
	fmsk = [];
	bmsk = [];
        [isp, pn, stati] = periodogram_test_periodicity_2d(...
					b, L,nf_test, bmsk, fmsk);
	r2 = wcorr(fr(nskip+1:end),Sr(nskip+1:end,1),fr(nskip+1:end),Sr(nskip+1:end,2)).^2;


	printf('Src/lc: %g\n',0); %sp(idx).stat.fc.radial.(field)*sp(idx).stat.Sc.radial.(field));
	printf('R2 %g\n',r2); %sp(idx).stat.fit.radial.bandpass.stat.goodness.r2);
	%printf('Sxc/lc: %g\n',sp(idx).stat.fc.x.(field)*sp(idx).stat.Sc.x.(field));
	%printf('R2 %g\n',sp(idx).stat.fit.x.phase_drift.stat.goodness.r2);
	printf('p-periodic %g\n',pn);
	
	%ylabel('Density S/\lambda_l');
	%legend('empirical','LP-fit');
	if (pflag)
		pdfprint(11,'img/pattern-irregular.pdf',ps)
		pdfprint(12,'img/pattern-irregular-Shat.pdf',ps)
		pdfprint(13,'img/pattern-irregular-correlogram.pdf',ps);
		pdfprint(14,'img/pattern-irregular-Sr.pdf',ps);
		pdfprint(15,'img/pattern-irregular-St.pdf',ps);
	end
end % function plot_pattern_observed_irregular
	
