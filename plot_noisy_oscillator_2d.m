% Sat 11 Jun 17:31:52 CEST 2022
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
% plot the spectral density, autocorrelation and a corresponding pattern
% for the Brownian phase filter
%
function plot_phase_noise_integration_2d(meta)
	if (nargin()<1)
		meta = pattern_formation_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	fcmap = meta.fcmap;
	col = meta.colormap;
	ps = meta.plotscale;
	linewidth = 1.5;

	mode = 's2d';
%	mode = 'exact';

	% length of domain
	L   = 20*[1,1];
	% spectral resolution
	df  = 1./L;

	% characteristic frequency fc = 1/\lambda_c
	fc=1;
	% number of grid points
	%n= (20^2)*[1,1];
	n= (256)*[1,1];
	% parameter influencing regularity
	sx0 = 0.1;

	% initialize the random number generator
	rng_seed = 5;

	% regularity in direction perbendicular to stripes
	%Sxpc_fc = 1.5*[1,2,3];
	%Sxpc_fc = [1,2,3];
	Sxpc_fc = [1,2,4];

	% regularity in direction parallel to stripes
	%Scy_fc = [1,2,3];
	Scy_fc = [1,2,4];

	% spatial axes
	x = linspace(0,L(1),n(1));
	y = linspace(0,L(2),n(2));
	% frequency axes
	fx = fourier_axis(L(1),(n(1)-1)+1);
	fxp = fx(1:n(2)/2);
	fy = fourier_axis(L(2),(n(2)-1)+1);
	fyp = fy(1:n(2)/2);
	
%	mSxpcfun = @(sx) max(spectral_density_brownian_phase(fx_,f0,sx,true)); 
	
	ns = [length(Sxpc_fc),length(Scy_fc)];
	
	% initialize arrays
	Sxp = zeros(n(1)/2,ns(1));
	Syp = zeros(n(2)/2,ns(2));
	Rx = zeros(n(1),ns(1));
	Ry = zeros(n(2),ns(2));
	Scy_ = zeros(ns(2),1);

	% compute densities and autocorrelation in the direction perpendicular to stripes
	for idx=1:ns(1)
		% sx(idx) = fzero(@(sx) f0*mSxpcfun(sx) - Sxpc_fc(idx),sx0);
		[f0(idx),sx(idx)] = phase_drift_pdf_mode2par(fc,Sxpc_fc(idx)/fc);
	end

	% compute densities and autocorrelation in the direction parallel to stripes
	for idx=1:ns(2)
		if (strcmp(mode,'s1d'))
			sc = 1;
		else
			sc = 4.5^2;
		end
		sy(idx)    = phase_drift_parallel_pdf_max2par(sc*Scy_fc(idx)/fc);
	end

	% 2D-plots
	% vary regularity in direction perpendicular to stripes
	for idx=1:ns(1)
		% vary regularity in direction parallel to stripes
		for jdx=1:ns(2)
			rng(rng_seed);
			try
			[b,xy,S,f,R] = anisotropic_pattern(L,n,f0(idx),[sx(idx),sy(jdx)],mode);
			% note, midpoint rule is more accurate for Syp as bin zero should be counted 1/2 time
			Sxp(:,idx)  = sum(S.S2d(:,1:n(1)/2),1)'*df(1);
			Syp(:,jdx)  = sum(S.S2d(1:n(2)/2,:),2)*df(2);
			Rx(:,idx)  = R(1,:)';
			Ry(:,jdx)  = R(:,1);
			hatS        = abs(fft2(b)).^2;
			hatSxp(:,idx) = mean(hatS,1)';
			hatSxp(:,idx) = 2*hatSxp(:,idx)/(sum(hatSxp(:,idx)*df(1)));
			hatSyp(:,jdx) = mean(hatS,2);
			hatSyp(:,jdx) = 2*hatSyp(:,jdx)/(sum(hatSyp(:,jdx)*df(1)));

			catch
			end

			R2d = fftshift(real(ifft2(S.S2d)));
						
			splitfigure([ns(1),ns(2)],[1,idx+(jdx-1)*ns(1)],fflag);
			cla();
			q = quantile(b,0.65);
			imagesc(x*fc,y*fc,b>q);
			xlabel('Position $x/\lambda_c$','interpreter','latex');
			ylabel('Position $y/\lambda_c$','interpreter','latex');
			axis square;
			axis equal
			axis tight
			colormap(fcmap(256));
			if (~pflag)
				title(num2str([Sxpc_fc(idx), Scy_fc(jdx)]))
			end
			%title(sprintf('S_c = %g, S_{cx} = %0.2g, s_x/s_y = %g', Sc(idx), Sxpc, sx/sy));
			
			splitfigure([ns(1),ns(2)],[2,idx+(jdx-1)*ns(1)],fflag);
			cla
			imagesc(fftshift(f.x)/fc,fftshift(f.y)/fc,fftshift(S.S2d));
			axis square;
			axis equal;
			axis tight;
			colormap gray;
			axis(2.5*[-1,1,-1,1])
			% e = rand(n,n);
			 
			splitfigure([ns(1),ns(2)],[3,idx+(jdx-1)*ns(1)],fflag);
			cla
			imagesc(R2d);
			axis equal
			axis square
			axis tight
			 
			 disp(sqrt(max(S.S2d(:)))*fc)
		end % for jdx
	end % for idx
	
	% 1D density perpendicular to stripes
	splitfigure([2,2],[4,1],fflag)
	cla
	plot(fxp/fc,Sxp*fc,'linewidth',linewidth)
if (0)
	hold on
	set(gca,'colororderindex',1)
	plot(fx/fc,hatSxp*fc,'*','linewidth',linewidth)
end
	xlim([0,2])
	xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
	ylabel('Density $S_{x^+}/\lambda_c$','interpreter','latex');
	if (~pflag)
		title('Perpendicular to stripes');
	end
	lh = legend(num2str(cvec(Sxpc_fc)),'location','northeast');
	%title(lh,'S_{x^+c}/\lambda_c');
	title(lh,'$S_{x^+c}/\lambda_c$','interpreter','latex');
	axis square
	pos = lh.Position;
	if (0) % if rectangular)
		lh.Position = [pos(1)+0.025,pos(2)+0.025, pos(3:4)];
	end
	set(gca,'colororder',col);

%	splitfigure([2,2],[5,1],fflag)
%	cla
%	plot(fx/fc,Sx*fc,'linewidth',linewidth)
	
	% 1D density parallel to stripes
	splitfigure([2,2],[4,2],fflag)
	cla()
	plot(fyp/fc,Syp*fc,'linewidth',linewidth)
%	fdx=(fy>=0);
%	df = fy(2)-fy(1)
%	sum(Syp)*df(2)
%	sum(Syp(fdx,:))*df(2)
%pause
if (0)
	hold on
	set(gca,'colororderindex',1)
	plot(fx/fc,hatSyp*fc,'*','linewidth',linewidth)
end
	xlim([0,2])
	xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
	ylabel('Density $S_{y^+}/\lambda_c$','interpreter','latex');
	if (~pflag)
		title('Parallel to stripes');
	end
	lh = legend(num2str(cvec(Scy_fc)));
	title(lh,'$S_{y^+c}/\lambda_c$','interpreter','latex');
	axis square
	set(gca,'colororder',col);
	
	% 1D autocorrelation perpendicular to stripes
	splitfigure([2,2],[4,3],fflag);
	cla
	plot(x*fc,Rx,'linewidth',linewidth);
	xlim([0,2]);
	xlabel('Lag distance $x/\lambda_c$','interpreter','latex');
	ylabel('Autocorrelation $R_x$','interpreter','latex');
	axis square
	set(gca,'colororder',col);
	
	% 1D autocorrelation parallel to stripes
	splitfigure([2,2],[4,4],fflag);
	cla
	plot(y*fc,Ry,'linewidth',linewidth);
	xlabel('Lag distance $y/\lambda_c$','interpreter','latex');
	ylabel('Autocorrelation $R_y$','interpreter','latex');
	xlim([0,10])
	axis square
	set(gca,'colororder',col);

	
	if (pflag)
		aspect = 1;
		s = [];
		% 2D patterns
		for idx=1:ns(1)
			for jdx=1:ns(2)
				pdfprint(10+idx+(jdx-1)*ns(1),['img/phase-noise-integration-2d-pattern-',mode,'-Sxpc-', ...
					 num2str(Sxpc_fc(idx)),'-Scy-',num2str(Scy_fc(jdx)),'.pdf'],ps,aspect,[],s);
			end % for jdx
		end % for idx
		s = [];
		aspect = 1
if (1)
		pdfprint(41,'img/phase-noise-integration-2d-density-Sx.pdf',ps,aspect);
		pdfprint(42,'img/phase-noise-integration-2d-density-Sy.pdf',ps,aspect);
		pdfprint(43,'img/phase-noise-integration-2d-autocorrelation-Rx.pdf',ps,aspect);
		pdfprint(44,'img/phase-noise-integration-2d-autocorrelation-Ry.pdf',ps,aspect);
end
	end % if pflag
end % plot_phase_noise_integration_2d

