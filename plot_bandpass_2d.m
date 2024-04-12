% 2022-06-12 11:37:43.189554667 +0200
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
%% plot bandpass generated patterns and their spectral density for with varying
%% degrees of regularity
%
function plot_bandpass_2d(meta)
	
	if (nargin()<1)
		meta=pattern_formation_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	plot_2d_only = pflag;
	fcmap = meta.fcmap;
	ps = meta.plotscale;

	normalize = true;
	lw = 1;
	% frequency of maximum of spectral density
	% n.b.: this still not works well when fc = 1/100, apparently accuracy is reduced because of suboptimal scaling in the hanke implementation
	fc = 1/10;
	% domain size
	L = 20/fc;
	% number of grid points
	n = 20*L*fc;
	% x-axis
	x = L*(0:n-1)'/n;
	% plot size
	Lp = 10;
	% regularity
	reg = [0.25,sqrt(2*0.25),2];
	reg = [1/sqrt(2)^3,1,sqrt(2)^3];
	reg = [0.5,1,2];
	%reg = [0.5,1,2]/sqrt(2);
	% regularity parameter of default density (reg=1)
	p = 2;
	
	b_thresh = 0.8;

	% initialize random number generator (for reproducibility)
	rng(1);
	
	df = 1/L;
	[fx, fy, fr] = fourier_axis_2d([L,L],[n,n]);

	L1     = L;
	n1     = n;
	x1 = L1*(0:n-1)'/n1;
	fr1    = fourier_axis(L1,n1);
	fdx = fr1>0;

	% normal noise
	e = randn(n,n);
	
	%e = single(e);
	%S2d = single(S2d);
	clf
	% create plots for vayring degrees of regularity
	for idx=1:length(reg)
		% reg = Sc/lc = Sc*fc -> Sc = reg/fc
		p2d = bandpass2d_continuous_pdf_mode2par(fc,reg(idx)/fc);
		printf('f0/fc %f fc %f Sc %f p_2d %f %f\n',p2d(1)/fc,fc,reg(idx)/fc,p2d);
		if (1)
			S2dp = bandpass2d_discrete_pdf([L,L],[n,n],1./p2d(1)*[1,1],p2d(2));
		%1./[fc,fc]/(2*pi/0.75),1);
		else
			S2d  = bandpass1d_continuous_pdf(fr,fc,p,normalize);
			S2d  = S2d/max(S2d(:));
			% transform density into arbitrary regularity
			S2dp = (S2d).^(p);
			% normalize
			S2dp  = 2*S2dp/sum(S2dp(:)*df*df);
		end

		% determine filter parameter yielding pattern with the desired regularity
		p = bandpass1d_continuous_pdf_max2par(fc,reg(idx)/fc);
		p = abs(p);
	
		% generate 1D (radial) density
		Sr1_     = bandpass1d_continuous_pdf(fr1,fc,p,normalize);

		% Transfer function
		T = sqrt(S2dp);
		% generate a pattern with the same specttral density by filtering
		b = real(ifft2(T.*e));
	
		% autocorrelation function
		R2d = ifft2(S2dp); %Sr1_);
		R2d = R2d/R2d(1);
		%R = R/R(1);
		Rr = R2d(1,:);

		% plot the pattern
		splitfigure([2,3],[1,idx], fflag);
		cla();
		b_bw = b>quantile(b(:),b_thresh);
		imagesc(x*fc,x*fc,b_bw);
		%q = quantile(abs(b(:)),0.01);
		%clim(q*[-1,1]);
		axis square
		axis equal
		axis tight
		axis xy
		axis([0,1,0,1]*Lp)
		colormap(gray)
		title(sprintf('$S_{rc}/\\lambda_c = %g$',roundn(reg(idx),2)),'interpreter','latex');
		xlabel('x/\lambda_c');
		ylabel('y/\lambda_c');
		colormap(fcmap(356));
		
		% plot th2 2D spectral density
		splitfigure([2,3],[1,3+idx], fflag);
		cla
		imagesc(fftshift(fx/fc),fftshift(fx/fc),fftshift(S2dp)*fc^2)
		colorbar
		axis(3.5*[-1,1,-1,1])

		%figure(2)
		% plot the radial (1D) spectral density
		splitfigure([2,2],[2,1], fflag);
		if (1==idx) cla; end
		%Src = max(Sr1_);
		[Sr_,fr_] = periodogram_radial(S2dp,[L,L]);
		if (~plot_2d_only)
			plot(fr1(fdx)/fc,(Sr1_(fdx))*fc,'linewidth',lw);
			hold on
			set(gca,'colororderindex',idx)
			plot(fr_/fc,Sr_.normalized*fc,'--');
		else
			plot(fr_/fc,Sr_.normalized*fc,'-','linewidth',1);
			
		end

%		xlim([0,2.5])
		axis square
		%axis equal
		axis tight
		xlim([0,3.5])
		hold on
		xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
		ylabel('Radial density $S_r/\lambda_c$','interpreter','latex');
		set(gca,'colororder',meta.colororder);
		leg_C = arrayfun(@(x) sprintf('%0.2f',x),reg,'uniformoutput',false);
		lh = legend(leg_C{:});%sprintf('%0.2f\n',cvec(reg)));
		%,num2str(roundn(cvec(reg),2)));
		title(lh,'Regularity $S_{rc}/\lambda_c$','interpreter','latex')
		if (idx==length(reg))	
			pos = lh.Position;
			pos(1)=0.84*pos(1);
			pos(2)=0.91*pos(2);
			set(lh,'Position',pos);
		end

		splitfigure([2,2],[2,2], fflag);
		if (1==idx) cla; end
		fdx = (x1*fc)>0.15;
		plot(x1(fdx)*fc,pi*sqrt(cvec(x1(fdx))*fc).*cvec(Rr(fdx)),'linewidth',lw)
		%plot(x1*fc,cvec(Rr),'linewidth',lw)
		hold on
		axis square
		axis equal
		axis tight
		xlim([0,2.5])
		ylim([-0.7,1])
		ylabel({'Rescaled Radial','Autocorrelation $\pi \sqrt{r/\lambda_c}{\cdot}R_r$'},'interpreter','latex');
		xlabel('Lag distance $r/\lambda_c$','interpreter','latex');
		set(gca,'colororder',meta.colororder);
		axis square	
	
		disp([sqrt(max(S2dp(:)))*fc, sqrt(max(Sr1_))*fc,reg(idx)])
	end
	
	if (pflag)
%		aspect = 1;
		aspect = [];
		for idx=1:length(reg)
			pdfprint(10+idx,sprintf('img/bandpass-2d-pattern-Sc-%0.2f.pdf',reg(idx)),ps,aspect);
%			pdfprint(12,'img/bandpass-2d-pattern-Sc-0.5.pdf',ps,aspect);
%			pdfprint(13,'img/bandpass-2d-pattern-Sc-1.0.pdf',ps,aspect);
		end
		pdfprint(21,'img/bandpass-2d-radial-density.pdf',ps,aspect);
		pdfprint(22,'img/bandpass-2d-radial-acf.pdf',ps,aspect);
	end % if pflag
end % plot_bandpass_2d

