% Tue 31 May 19:23:38 CEST 2022
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
%% plot patterns generated with the rietkerk model
function rk_2d_heterogeneity_plot(rk,sp,out,figid,aniso,nkc,meta,visible)
	ps      = meta.plotscale;
	fflag   = meta.pflag;
	fcmap   = meta.fcmap;
	figid  = 10*figid;
	lx_overlay = 6;
	col_overlay = [ 1,  1,1;
			1,  0,0;
			0,0.0,1; % 0 0.4 1
			0.5,0,0.5];
	x   = {};
	[x{1},x{2}]   = rk.x();
	lc = sp.lambda_c;	
	l1c = lc;
	fc = 1./lc;

	% plot 2D
		% biomass 2d
		 splitfigure([3,3],[figid+0,1],fflag,num2str(rk.pss.a),[],[],[],'Visible',visible);
		 cla();
		 sp.plot('b');
		 	%imagesc((x{2}-x{2}(1))/lc,(x{1}-x{1}(1))/lc,(b>p)');
		 	%xlabel('$x/\lambda_c$','interpreter','latex');
		 	%ylabel('$y/\lambda_c$','interpreter','latex');
		 colormap(fcmap(256)); %meta.colormap_b);
		 axis equal
		 axis tight
		 axis xy
		 if (~meta.pflag)
			title('Pattern Biomass');
		 end % if meta.pflag	

		 % correlogram 2d
		 splitfigure([3,3],[figid+0,4],fflag,[],[],[],[],'Visible',visible);
		 cla();
		 sp.plot('R.hat'); 
		 axis(nkc*[-1,1,-1,1]);
		 if (nkc == 3)
			%set(gca,'xtick',[-2,0,2]);
			%set(gca,'ytick',[-2,0,2]);
		 else
			%set(gca,'xtick',-2:2);
			%set(gca,'ytick',-2:2);
		 end % if nck
		 axis square
		 axis xy
		 %colormap(fcmap(256)); %meta.colormap_b);
		 caxis([-0.3,1]);
		 colormap(meta.fcmap(10))

		 if (~meta.pflag)
			title('Correlogram 2D');
		 end % ~meta.pflag

		 splitfigure([3,3],[figid+0,7],fflag,[],[],[],[],'Visible',visible);
		 cla();
		 sp.plot('S.hat'); 
		 axis(nkc*[-1,1,-1,1]);
		 if (nkc == 3)
			set(gca,'xtick',[-2,0,2]);
			set(gca,'ytick',[-2,0,2]);
		 else
			set(gca,'xtick',-2:2);
			set(gca,'ytick',-2:2);
		 end % if nck
		 axis square
		 colormap(fcmap(256)); %meta.colormap_b);
		 axis xy
		 if (~meta.pflag)
			title('Periodogram 2D');
		 end % ~meta.pflag

		%if (0)
		% splitfigure([2,4],[figid+0,4],fflag,[],[],[],[],'Visible',visible);
		% cla();
		% imagesc(fftshift(f.x{1})*l1c,fftshift(f.x{2})*l1c,fftshift(S2(:,:,idx))/l1c);
		% axis(nkc*[-1,1,-1,1]);
		% title('Density 2D');
		%end


		% 1d-quantities
		% for x and y or r and angle
		for jdx=1:2

		% biomass 1d
		 splitfigure([3,3],[figid+0,1+jdx],fflag,[],[],[],[],'Visible',visible);
		 cla();	
		 if (jdx == 1)
		 plot(x{jdx}/lc,sp.b(round(end/2),:));
		 xlabel('x/\lambda_c');
		 else
		 plot(x{jdx}/lc,sp.b(:,round(end/2)));
		 xlabel('y/\lambda_c');
		end
		 ylabel('b');
		 title('Pattern cross section');
		 axis square

		% 1d density
		 splitfigure([3,3],[figid+0,4+jdx],fflag,[],[],[],[],'Visible',visible);
		 cla();
		if (aniso)
		 if (1==jdx)
			 sp.plot('S.rot.x.hat');
			 hold on
			 %sp.plot('S.rot.x.hp');
			 sp.plot('S.rot.x.phase_drift');
			 %legend('RK-model','PNI-fit');
%		 	 legend('RK' ...
%		          , ['BM R^2 =', num2str(round(sp.stat.fit.x.phase_drift.stat.r2,3))] ...
%		          , ['BP R^2 =', num2str(round(sp.stat.fit.x.bandpass.stat.r2,3))] ...
%		         );
		 else
			 sp.plot('S.rot.y.hat');
			 hold on
			 sp.plot('S.rot.y.phase_drift_parallel');
			 %legend('RK-model','PNI-fit');
%		 	 legend('RK' ...
%		          , ['BM R^2 =', num2str(round(sp.stat.fit.y.phase_drift_parallel.stat.r2,3))] ...
%		         );
		 end % else of if jdx==1
		else

		if (1==jdx)
		 sp.plot('S.radial.hat');
		 hold on
		 sp.plot('S.radial.bandpass');
		 %legend('RK-model','BP-fit');
		 xlim([0,2.5]);
		else
		 sp.plot('S.rot.angular.hat');
		 hold on
		 hline(1/pi,'linewidth',1,'color','r');
		end
%		 legend('RK' ...
%		        , ['BM R^2 =', num2str(round(sp.stat.fit.radial.phase_drift.stat.r2,3))] ...
%		        , ['BP R^2 =', num2str(round(sp.stat.fit.radial.bandpass.stat.r2,3))] ...
%		       );

		end % if aniso
		hold on
		axis square
		set(gca,'colororder',meta.colororder);
	
		% 1D autocorrelation
		splitfigure([3,3],[figid+0,7+jdx],fflag,[],[],[],[],'Visible',visible);
		cla();
		if (aniso)
			if (1 == jdx)
				plot(sp.x/sp.lambda_c,fftshift(sp.R.rot.x.hat));
				xlabel('Lag Distance $x/\lambda_c$','interpreter','latex');
			else
				plot(sp.y/sp.lambda_c,sp.R.rot.y.hat);
				xlabel('Lag Distance $y/\lambda_c$','interpreter','latex');
			end
		else
			if (1==jdx)
				plot(sp.r/sp.lambda_c,sp.R.radial.hat);
				ylabel('$R_r$','interpreter','latex');
				xlabel('Lag Distance $r/\lambda_c$','interpreter','latex');
			end
		end % else of if band
		xlim([0,2.5]);
		axis square
	end % for jdx
	
	% plot overlay of rad and bandpass generated pattern
	splitfigure([3,3],[figid+1,1],fflag,[],[],[],[],'Visible',visible);
	imagesc(out.b_thresh);
	axis equal
	title('RK');
	drawnow

	splitfigure([3,3],[figid+1,2],fflag,[],[],[],[],'Visible',visible);
	imagesc(out.b_bp_thresh);
	axis equal
	title('BP');	
	drawnow

	splitfigure([3,3],[figid+1,3],fflag,[],[],[],[],'Visible',visible);
	cla();
	for pdx=1:4
		plot(NaN(1),NaN(1),'ko','markerfacecolor',col_overlay(pdx,:));
		hold on
	end			
	imagesc(x{1}/l1c,x{2}/l1c,out.b_overlay);
	colormap(col_overlay)
	if (~meta.pflag)
		title(num2str(out.cba))
	end
	%legend('bare','RK only','BP only','overlap')
	axis square;
	axis tight;
	xlim([0,lx_overlay]);
	ylim([0,lx_overlay]);
	xlabel('$x/\lambda_c$','interpreter','latex');
	ylabel('$y/\lambda_c$','interpreter','latex');
	drawnow

	% plot histogram
	splitfigure([3,3],[figid+1,4],fflag,[],[],[],[],'Visible',visible);
	cla();
	plot(out.patch_size.radius,out.patch_size.pr);
	ylabel('patch radius distribution')
	drawnow

	splitfigure([3,3],[figid+1,5],fflag,[],[],[],[],'Visible',visible);
	plot(out.bhist.x,out.bhist.p)
	title('histogram of log(biomass)');
	drawnow

	% plot coherence
	splitfigure([3,3],[figid+1,7],fflag,[],[],[],[],'Visible',visible);
	cla();
	fx = sp.f.x;
	fy = sp.f.y;
	imagesc(fftshift(fx/fc),fftshift(fy/fc),fftshift(out.coherence_2d));
	axis equal;
	axis tight;
	title('Spectral coherence C(b,a)');
	xlabel('k_x/k_c');
	xlabel('k_y/k_c');
	drawnow

	if (~isempty(out.coherence_r))
	splitfigure([3,3],[figid+1,8],fflag,[],[],[],[],'Visible',visible);
	cla();
	plot(sp.f.r/fc,out.coherence_r);
	vline(1);
	ylabel('C_r(a,b)');
	xlabel('k_x/k_c');
	drawnow
	end
	
	if (0)
	% evolution of regularity, Sc and lc over time

	fun = @(t,p) p(1)*(1-exp(-p(2)*t));
	hold on
	p  = lsqnonlin(@(p) double(fun(t,p) - Sc),double([Sc(end),1e-4]));

	%namedfigure(figid+3,'','Visible',visible);
	splitfigure([2,3],[figid+2,4],fflag,[],[],[],[],'Visible',visible);
	%subplot(2,2,1);
	cla();
	plot(t,Sc/lc(end));
	xlabel('t/day');
	ylabel('Sc/lc');
	%title(p);

	%subplot(2,2,2);
	splitfigure([2,3],[figid+2,5],fflag,[],[],[],[],'Visible',visible);
	cla();
	plot(t,lc);
	xlabel('t/day');
	ylabel('lc');

	%subplot(2,3,4);
	splitfigure([2,3],[figid+2,6],fflag,[],[],[],[],'Visible',visible);
	cla();
	r = rms(y(:,1:end/3)-y(end,1:end/3),2)./rms(y(end,1:end/3));
	r(:,2) = inner2outer(rms(diff(y(:,1:end/3)),2))./rms(y(end,1:end/3));
	semilogy(t,r);

	end

	if (meta.pflag)
		T = rk.T(end);

		str = sprintf('vxh-%g-eyh-%g-s_a-%0.2g-R-%g-L-%d-T-%1.0e-rng-%d',rk.pmu.vx(3),rk.pmu.ey(3),rk.pss.a/rk.pmu.a,rk.pmu.R,rk.L(1),T,rk.opt.rng);		
		pdfprint(figid*10+1,['img/rk-2d-pattern-',str,'.pdf'],ps);
		%namedfigure(figid*10+1,'','Visible',visible);
		%axis(10*[0,1,0,1]);
		figure(figid*10+1);
		drawnow();
		axis off
		pdfprint(figid*10+1,['img/rk-2d-pattern-',str,'.png'],ps,[],'png');
		axis on
		xlim(10*[0,1]);
		ylim(10*[0,1]);
		drawnow();
		pdfprint(figid*10+1,['img/rk-2d-pattern-Lp-10-',str,'.pdf'],ps);
		%namedfigure(figid*10+1,'','Visible',visible);
		axis off
		pdfprint(figid*10+1,['img/rk-2d-pattern-Lp-10-',str,'.png'],ps,[],'png');
		pdfprint(figid*10+4,['img/rk-2d-correlogram-',str,'.pdf'],ps);
		pdfprint(figid*10+7,['img/rk-2d-periodogram-',str,'.pdf'],ps);
%		namedfigure(NaN*figid*10+3,'','Visible',visible);
%		axis off
%		pdfprint(NaN*figid*10+3,['img/rk-2d-periodogram-',str,'.png'],ps,[],'png');
		aspect = [];
		if (aniso)
			pdfprint(figid*10+5,['img/rk-2d-density-Sx-',str,'.pdf'],ps,aspect);
			pdfprint(figid*10+6,['img/rk-2d-density-Sy-',str,'.pdf'],ps,aspect);
			pdfprint(figid*10+8,['img/rk-2d-acf-Rx-',str,'.pdf'],ps,aspect);
			pdfprint(figid*10+9,['img/rk-2d-acf-Ry-',str,'.pdf'],ps,aspect);
		else
			pdfprint(figid*10+5,['img/rk-2d-density-Sr-',str,'.pdf'],ps,aspect);
			pdfprint(figid*10+6,['img/rk-2d-density-Sa-',str,'.pdf'],ps,aspect);
			pdfprint(figid*10+8,['img/rk-2d-acf-Rr-',str,'.pdf'],ps,aspect);
			%pdfprint(figid*10+9,['img/rk-2d-acf-Ry-',str,'.pdf'],ps,aspect);
		end
		if (rk.pss.a>0)
			pdfprint((figid+1)*10+3,['img/rk-2d-pattern-vs-bandpass-overlap/rk-2d-pattern-vs-bandpass-overlap-',str],ps);
		end
	end % meta.pflag
end % rk_2d_heterogeneity_plot

