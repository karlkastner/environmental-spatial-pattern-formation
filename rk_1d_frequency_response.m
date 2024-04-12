% 2022-06-26 10:04:34.344120984 +0200 asinc.m
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
% demonstrate similarity of the rk-model to a bandpass filter
%
function [rk,bb,c] = rk_1d_frequency_response(meta)
	if (nargin<1)
		meta = pattern_formation_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	%ps    = 3.5;
	ps    = meta.plotscale;

	scale = 8;
	
	% length of domain
	L = 16000*scale;
	% spatial resolution
	dx = 1;
	% number of grid points
	nx = round(L./dx);
	a0 = 0.2;
	% relative magnitude of spatial variation of the bare soild infiltration 
	s_a = 0.05;
	% window size
	nwin = 2*round(800*sqrt(scale))+1;
	nwh = (nwin-1)/2;

	% final time
	To  = 4e2;
	% output time step
	dto = 20; 
	% time step
	dt  = 0.5;

	% frequencies at which pattern is plotted
	f0 = [0.5,1,2];

	% local frequency
	k = 20;
	r = 1-1/8;
	f  = 12*r.^(0:k/r)';
	f = f*2;
	F_ = scale*[0;cumsum(f)];
	F = mid(F_);

	% phi = omega*x <-> omega = unwrapp(phi) / x 
	
	% interpolate to grid
	n  = 0.5*nx+1;
	n_ = length(F);
	F  = interp1((1:n_)/(n_+1),F,(0:n-1)/(n-1),'spline')';
	df_dx = cdiff(F)./dx;
	% local bare soil infiltration
	ea  = sin(2*pi*(F-F(1)));

	% when we only increase the frequency, then there will be a jump at the end of the domain
	% so we ramp it up and down to make it symmetric
	ea = [ea;-flipud(ea(2:end-1))];
	a = a0*(1 + s_a*ea);  
	df_dx_ = [df_dx;-flipud(df_dx(2:end-1,:))];
	
	% model parameters
	param        = struct();
	% advection
	param.pmu.vx = [0,0,0];
	% diffusion
	param.pmu.ex = [1,1,100];
	% precipitation
	param.pmu.R  = 1.0;
	% boundary condition
	param.boundary_condition    = {'circular','ciruclar'};
	% reload values of intermediate time steps
	param.opt.loadfinal = false;
	param.L = L;
	param.nx = nx;
	param.T = To;
	param.opt.dt = dt;
	param.opt.dto = dto;
	param.pss.a  = 0;
	param.opt.compute_class = @single;
	param.opt.path_str = 'mat/';
	param.pmu.a = a;

	% initial condition
	rk = Rietkerk(param);
	[b0,w0,h0] = rk.homogeneous_state(rk.pmu);
	rng(0);
	ic = [rand(nx,1);
              w0*ones(nx,1);
	      mean(h0)*ones(nx,1)];
	%param.initial_condition = 'obj.random_state()';
	param.initial_condition = ic;

	% run model
	rk = Rietkerk(param);
	[t,y]	  = rk.run(); %param);
	y = single(y);
	[bb,ww,hh] = rk.extract1(y);
	b = bb(end,:)';
	
	%ax = plotyy(rk.x,[bb(end,:)'/(2*rms(bb(end,:))),0.5*(1+(rk.pmu.a-a0)/(s_a*a0))],rk.x,ff)
	
	% periodically extend, as the first and last window exceed the range
	b = [b;b;b];
	a = [a;a;a];

	% goodness of fit (correlation)
	c=NaN(rk.nx,2);
	w=triwin(1:nwin)'; % triwin
	wh=hanwin(1:nwin)'; % triwin
	for idx=1:rk.nx
		% rectwin
		c(idx,1) = corr(b(rk.nx+idx+(-nwh:nwh)),a(rk.nx+idx+(-nwh:nwh)));
		% triwin
		c(idx,2) = wcorr(w,b(rk.nx+idx+(-nwh:nwh)),w,a(rk.nx+idx+(-nwh:nwh)));
		c(idx,3) = wcorr(wh,b(rk.nx+idx+(-nwh:nwh)),wh,a(rk.nx+idx+(-nwh:nwh)));
	end
	nbartlett = round(sqrt(rk.nx));
	S = periodogram_bartlett(bb(end,:)-mean(bb(end,:)), rk.L,nbartlett,rk.nx);
	fx = fourier_axis(rk.x);
	[Sc,mdx] = max(S);
	fc = fx(mdx);
	%Sa = periodogram_bartlett(rk.pmu.a-mean(rk.pmu.a), rk.L, 40*scale,rk.nx);
	%[Sca,mdxa] = max(Sa);
	%fca = fx(mdxa);

	% plot pattern and frequency
	splitfigure([2,3],[1,1],fflag);
	yyaxis left
	cla
	plot(rk.x,[bb(end,:)'/(2*rms(bb(end,:))),0.5*(1+ea)]);
	%set(gca,'colororder',{'k','r'})
	set(gca,'colororder',[0,0,0;1 0 0])  
	yyaxis right
	cla
	plot(rk.x,df_dx_);
	hline(fc)

	% plot goodness of fit (correlation^2)
	splitfigure([2,3],[1,2],fflag);
	if (~pflag)
		plot(df_dx_(1:nx)/fc,(c(1:nx,:)).^2,'linewidth',1)
	else
		plot(df_dx_(1:nx)/fc,(c(1:nx,3)).^2,'linewidth',1)
	end
	%plot(df_dx(1:nx/2)/fc,(c(1:nx/2,:)).^2,'linewidth',1)
	ylabel('Goodness of fit $R^2(a,b)$','interpreter','latex');
	xlabel('Wavenumber $k_a/k_c$','interpreter','latex');
	axis tight
	ylim([0,1])
	xlim([0,3.5])
	if (pflag)
		axis square
	end

	% plot spectrum	
	splitfigure([2,3],[1,3],fflag);
	plot(fx/(fc+sqrt(eps)),S*(fc+sqrt(eps)));  
	xlim([0,3.5]); %+sqrt(eps))]);
	xlabel('Wavenumber $k_a/k_c$','interpreter','latex');
	ylabel('S/\lambda_c');
	title(num2str(fc))	
	drawnow

	% plot parts of pattern
	d = [0,0,0];
	x = rk.x;
	dt = [0, 0, 0];
	for idx=1:length(f0)
		[mf,mdx] = min(abs(df_dx-f0(idx)*fc));
		splitfigure([2,3],[1,3+idx],fflag);
		drawnow
		%cla
		%[ax,h1,h2] = plotyy( (rk.x-x(mdx))*fc+d(idx),bb(end,:)'/(rms(bb(end,:))), ...
		%	(rk.x-x(mdx)+d(idx))*fc,rk.pmu.a/mean(rk.pmu.a));
		yyaxis left
		drawnow
		ax = gca;
		cla(ax)
		plot( (rk.x-x(mdx))*fc+d(idx),bb(end,:)'/(rms(bb(end,:))),'linewidth',1);
		xlim(ax(1),2*[-1,1]);
		ylim(ax(1),[0,2.3])
		set(ax(1),'xtick',[-2:2]);
		set(ax(1),'ytick',[0,1,2]);
		xlabel(ax(1),'$x/\lambda_c$','interpreter','latex');
		if (idx == 1)
			ylabel(ax(1),'Biomass $b/\mathrm{rms}(b)$','interpreter','latex')
		else
			set(ax(1),'yticklabel',[]);
		end
		text(ax(1),-2+dt(idx) + 0.1,2.2 - 0.025,sprintf('$k_a = %g \\cdot k_c$',f0(idx)),'interpreter','latex'); 
		ax(1).YColor = [0,0,0.7];
		set(ax(1),'colororder',[0,0,0.7;0.8,0,0]) 
		axis(ax(1),'square')
		drawnow		

		yyaxis right
		drawnow
		ax(2) = gca;
		cla(ax(2))
		plot( (rk.x-x(mdx)+d(idx))*fc,rk.pmu.a/mean(rk.pmu.a),'linewidth',1);
		hold on
		yr = 0.95+s_a*bb(end,:)'/(rms(bb(end,:)));
		plot( (rk.x-x(mdx)+d(idx))*fc,yr,'--','color',[0,0,0.8],'linewidth',1);
		linkaxes(ax,'off');
%		yyaxis left
%		hold on
		%hold(ax(1),'on')
%		col = colormap();
%		plot((rk.x-x(mdx))*fc+d(idx),bb(end,:)'/(rms(bb(end,:))),'--',	'color',[0,1,0]); %col(1,:));
		%h1.LineWidth = 1;
		%h2.LineWidth = 1;
		%set(gca, 'SortMethod', 'depth')

		%(min(rk.pmu.a) + [-sqrt(eps),1.25*range(rk.pmu.a)])/mean(rk.pmu.a));
		%ylim(ax(2), (min(rk.pmu.a) + [-sqrt(eps),1.25*range(rk.pmu.a)])/mean(rk.pmu.a));

		xlim(ax(2),2*[-1,1]);
		ylim(ax(2), 1+[-(s_a+sqrt(eps)),+1.3*s_a]);
		set(ax(2),'ytick',[0.95,1,1.05])
		set(ax(2),'colororder',[0.8,0,0;0,0,0.7]) 
		ax(1).YColor = [0.8,0,0];
		%ylim(ax(2), 1+[-(s_a+sqrt(eps)),+1.1*s_a]);
		if (idx==3)
			ylabel(ax(2),'Infiltration Coefficient $a/\bar a$','interpreter','latex')
		else
			ylabel(ax(2),[]);
			set(ax(2),'yticklabel',[]);
		end
		axis(ax(2),'square')
		drawnow
	end
	
	if (pflag)
		pdfprint(12,'img/rk-frequency-response-corr-a-b.pdf',ps);
		for idx=1:3
			pdfprint(13+idx,['img/rk-frequency-response-fa-',num2str(f0(idx)),'.pdf'],ps);
		end	
	end
end % rk_bandpass_similarity

