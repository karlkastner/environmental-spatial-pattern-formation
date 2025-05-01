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
%
%
% we can just take one long pattern where that we split
% note: negative corr might be due to ambiguity of initial phase shift
%       if the first patch is off by 180 deg, the sign inverts
function [stat, bb, c] = rk_1d_experiment_phase_noise_integration(meta)
	if (nargin<1)
		meta = pattern_formation_metadata();
	end
	dflag = 0;
	pflag = meta.pflag;
	fflag = pflag;
	ps    = meta.plotscale;

	scale = 1;
	% length of domain
	L   = scale*2000;
	x0  = 0;
	% spatial resolution
	dx  = 4; %scale;
	% time step
	dt  = 0.5; %scale;
	% number of grid points
	nx = round(L./dx);
	% final time
	To  = scale*3.2e4;
	% output time step
	dto = 1000;

	% mean infiltration coefficient
	a0 = 0.2;

	% relative magnitude of spatial variation of the bare soild infiltration
	s_a = [0,0.001,0.01];
	s_a = [0,0.1];
	
	% seeds for random number generator (for repro
	ns = 10;
	seed = ones(length(s_a),1)*(1:ns);


	fc0 = 1/100;

	% model parameters
	param        = struct();
	% advection
	param.pmu.vx = [0,0,-10];
	% diffusion
	param.pmu.ex = [1,1, 0.0];
	% precipitation
	param.pmu.R  = 1.0;
	% boundary condition
	param.boundary_condition    = {'circular','ciruclar'};
	% reload values of intermediate time steps
	param.opt.loadfinal = true;
	param.L = L;
	param.nx = nx;
	% length of simulation
	param.T = To;
	% time step
	param.opt.dt = dt;
	% output time step
	param.opt.dto = dto;
	% meand of infiltration coefficient
	param.pmu.a = a0;
	% standard deviation of infiltration coefficient
	param.pss.a  = 0;
	param.opt.compute_class = @single;
	param.opt.path_str = 'mat/rk-experiment-phase-noise-integration/';
	mkdir(param.opt.path_str);
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_spectral';
	param.opt.adapt_time_step = 0;

	param.initial_condition = 'obj.random_state()';

	rk = Rietkerk(param);
	rk.hashfield_C{end+1} = 'p.a';
	[b0,w0,h0] = rk.homogeneous_state();
	ic = struct();
	ic.mu = [15,3,6];
	ic.sd = [];
	ic.s  = [0,0,0];
	ic.c  = [15,0,0];
	ic.fc = fc0;

	param.initial_condition = ic;
	corr_ = NaN(length(s_a),1);
	corr__ = NaN(length(s_a),1);

	stat = struct();
	for adx=1:length(s_a)
	    for sdx=1:size(seed,2)
		disp([sdx adx]);
		if (s_a(adx)==0)
		rng(1);
		else
		rng(seed(adx,sdx));
		end
	
		% varying half
		sd_a = a0*s_a(adx);
	
		ne = sum(rk.x>=x0);

		%[lmu,lsd] = logn_moment2par(a0,sd_a);
		%a = lognrnd(lmu,lsd,ne,1);
		[p1,p2] = uniform_moment2par(a0,sd_a);
		a = uniformrnd(p1,p2,ne,1);

		% ensure mean and sd
		a = a0 + (sd_a/std(a))*(a-mean(a));
		% append non-varying half
		a = [repmat(a0,nx-ne,1);a];
		% make a double bridge
		%a(nx/2+1:end) = (a(nx/2+1:end) - mean(a(nx/2+1:end))+a0);
		param.p.a = a;
	
		% run model
		rk = Rietkerk(param);
		rk.hashfield_C{end+1} = 'p.a';
		[t,y]	  = rk.run();
		y = single(y);
		[bb,ww,hh] = rk.extract1(y);
		size(bb)
		b = bb(end,:)';
		b_(:,adx,sdx) = b;
		%rms(bb)
	
		fx = fourier_axis(rk.x);
		S = abs(fft2(b-mean(b)));
		[Sc,mdx] = max(S);
		%fc_(adx,sdx) = fx(mdx);
		[Sc_, fc(adx,sdx)] = extreme3(abs(fx),S,mdx);
	
		% plot pattern and frequency
		bold = [];
		for idx_ = 1
			ca=cumsum(rk.p.a-a0)*dx;
			ca(end)
			ca_ = [ca;ca;ca];
			nl = round(3/(fc(adx,sdx)*dx));
			%win = cvec(triwin(1:nl));
			win = rectwin(1:nl);
			ca_ = conv(ca_,win,'same');
			ca_ = ca_(nx+1:2*nx);
			% periodically extend
			b3 = [b;b;b];
			nx = rk.nx;
			x = cvec(rk.x);
			nl = round(5/(fc(adx,sdx)*dx));
			xl = x(1:nl);
			swin = sin(2*pi*fc(adx,sdx)*xl);
			cwin = cos(2*pi*fc(adx,sdx)*xl);
			%twin = triwin(1:nl)';
			twin = 1;
			s = conv(b3,swin.*twin,'same');
			c = conv(b3,cwin.*twin,'same');
			s=s(nx+1:2*nx);
			c=c(nx+1:2*nx);
			% phase
			phi = atan2(s,c);
			phi = unwrap(phi);
			phi = phi-phi(1);
			%fdx = find((x>3300/2) & (x<3700/2));
			fdx = 1:length(x);
		
			%fc(adx) = 1/(2*pi)*(phi(end)-phi(1))/(x(end)-x(1))
		end % for idx_
		%dphidx = (phi(fdx(end))-phi(fdx(1)))/(x(fdx(end))-x(fdx(1)));
		%dphi = phi - dphidx*x;
		if (adx == 1)
			phi0 = phi;
		end % if adx == 1
		dphi         = phi - 2*pi*fc(adx,sdx)*x;
		phi_(:,adx,sdx)  = phi; %/fc(adx,sdx);
		ddphi(:,adx,sdx) = dphi; %/fc(adx,sdx);
	
		corr_(adx,sdx) = corr(dphi,ca);
		corr__(adx,sdx) = corr(dphi,ca_);
		A = [ones(nx,1),ca_];
		cc(:,adx) = A \ dphi;
		dphi_ = A*cc(:,adx);
		r2(adx,sdx) = 1 - mean((dphi_-dphi).^2)/var(dphi);
	
		if (dflag)
		splitfigure([2,3],[1e3*adx+sdx,1],fflag);
		cla
		imagesc(bb)

		splitfigure([2,3],[1e3*adx+sdx,2],fflag);
		cla
		plot(fc(adx)*(rk.x-x0),b);
		if(1==adx)
			b0 = b;
		end
		hold on
		plot(fc(adx)*(rk.x-x0),[b0]);
		printf('R2 %f\n',r2(adx,sdx))
		title(s_a(adx))
	
		splitfigure([2,3],[1e3*adx+sdx,3],fflag);
		cla();
		%fdx_ = round((x0/L)*nx+1)-nl
		fdx_ = 1;
		%plot((rk.x-x0)/lc,[ca-ca(1),(ca_-ca_(fdx_))]/(lc*a0));
		plot((rk.x-x0)*fc(adx),[(ca_-ca_(fdx_))]/(a0/fc(adx)));
		ylabel({'Heterogeneity-Integral', ...
			'$(\int (a-\bar a) \mathrm{d} x)/(\lambda_c \bar a)$'},'interpreter','latex');
		%if (adx==2) ylim([-6,2]*1e-4*s_a(adx)/s_a(1)); end
		yl = ylim;
		ylim(max(abs(ylim))*[-1,1]);
		yyaxis right;
		plot((rk.x-x0)*fc(adx),(dphi-dphi(fdx_))/(2*pi));
		ylabel({'Phase Deviation','$(\phi - k_c x)/(2 \, \pi)$'},'interpreter','latex');
		xlabel('Position $x/\lambda_c$','interpreter','latex');
		%$ylim([-6,2]*1e-3*s_a(adx)/s_a(1));
		%if (adx==2) ylim([-6,2]*1e-3*s_a(adx)/s_a(1)); end
		%yl = ylim; ylim(max(abs(ylim)*[-1,1]));
		yl = ylim; ylim(max(abs(ylim))*[-1,1]);
		xlim([0,0.5*(L-x0)]*fc(adx));
		if (pflag && sdx == 5 && adx == 1)
			ps   = 3.5;
			name = sprintf('img/phase-random-walk-and-int-ea-dx-sa-%g-seed-%g.pdf',s_a(adx),seed(adx,sdx));
			pdfprint(10053,name,ps);
		else
			title(corr_(adx,sdx))
		end % if pflag
	if (1)
		%hist(log(b),100)
		% determine centroids of stripes
		isveg = 0;
		bveg  = 2.0; % quantile(b,0.2);
		bbare = 1.0; % quantile(b,0.4);
		xc = [];
		b = cvec(b);
		x = cvec(rk.x);
		for idx=1:nx
			if (~isveg)
				if (b(idx)>bveg)
					isveg = true;
					i1 = idx;
				end
			else
				if (b(idx)<bbare)
					isveg = false;
					i2 = idx;
					% determine centroid of current stripe
					xc(end+1,1) = sum(x(i1:i2).*b(i1:i2))./sum(b(i1:i2));
				end
			end % if ~isveg
		end % for idx
		% interpolate patch number to x
		np = interp1(xc,1:length(xc),x,'linear','extrap');
	end % if 1
		dnp_dx = (np(fdx(end))-np(fdx(1)))/(x(fdx(end))-x(fdx(1)));
		dnp_dx = (np(end)-np(end/2))/(x(end)-x(end/2));
		dnp = np-dnp_dx*x;
		if (1 == adx)
		np0 = np;
		else
			dnp = np - np0;
		end % else of is adx == 1
		corr__(adx,sdx) = corr(dnp,ca);
		%[cvec(s_a), corr_(:,1),corr__]
		nanmedian([cvec(s_a), corr_(:,1),corr__])
	
		splitfigure([2,3],[1e3*adx+sdx,4],fflag);
		plot(rk.x,[phi,2*pi*fc(adx)*rk.x']/(2*pi));
	
		splitfigure([2,3],[1e3*adx+sdx,5],fflag);
		plot(fx,S);
		vline(fc(adx),'color',[1,0,0])
	
		splitfigure([2,3],[1e3*adx+sdx,6],fflag);
		plot(rk.x,dnp);
		yyaxis right
		plot(rk.x,ca);
		end % if dflag
		bold = b;
	    end % for sdx 
	end % for adx

	figure(1e4);

	subplot(2,2,1);
	s_a
	corr_
	a = atanh(corr_);
 	ma = mean(a,2);
	sa=std(a,[],2);
	a=tanh(ma+sa*[-1,0,1]);
	%errorbar(s_a/a0,a(:,2),a(:,2)-a(:,1),a(:,3)-a(:,2));
	bar(s_a,a(:,2).^2);
	xlim([0,s_a(end)*1.01])
	xlabel('Exogenous heterogeneity $s_a$','interpreter','latex');
	ylabel('Goodness of fit $R^2$','interpreter','latex');
	%(\int (a-\bar a) \dx,phi)

	subplot(2,2,2)
	hist(corr_)

	subplot(2,2,3)
	hist(r2);
	stat.corr_ = corr_;
	stat.corr__ = corr__;
	stat.r2 = r2;
	stat.phi = phi_;
	stat.dphi = ddphi;
	stat.fc = fc;
	stat.s_a = s_a;
	stat.b = b_;
	stat.x = x;


	figure(1e5);
	bar(median(stat.r2(2:end,:),2),'facecolor',[0,0,0.6]);
	ylabel('Goodness of fit R^2');
	xlabel('Degree of heterogeneity $s_a$','interpreter','latex');
	set(gca,'xticklabel',stat.s_a(2:end));
	axis square;
	if (0) %pflag)
		pdfprint(1e5,'img/phase-prediction-r2.pdf',4);
	end

	fflag = 1;
	for idx=1:length(stat.s_a)
		d_=squeeze(stat.dphi(:,idx,:));
		fc=mean(stat.fc(idx,:));
		dphi =rms(d_,2);
	
		v=sqrt(x.*(x(end)-x));
		v=v/max(v);
	
		id = min(48,ns/2);
	
		splitfigure([2,3],[1e3,idx],fflag);
		cla
		if (pflag)
		errorarea2(x*fc,dphi*[-1,0,1]/(2*pi*stat.s_a(idx)),[1/2,1/2,1])
		else
		errorarea2(x*fc,dphi*[-1,0,1]/(2*pi*stat.s_a(idx)),[0,0,0])
		end
		hold on;
		plot(x*fc,[d_(:,id)]/(2*pi*stat.s_a(idx)),'color',[0,0,0.7],'linewidth',1);
		xlim([0,fc*x(end)/4]);
		ylabel('Phase Error $\Delta \phi / (2 \pi s_a)$','interpreter','latex');
		xlabel('Position $x/\lambda_c$','interpreter','latex');
		xlim([0,10])
	
		splitfigure([2,3],[1e3,idx+3],fflag);
		cla();
		dx=[-0.05,-0.25];
		xx=flat(stat.x+[-1,0,1]*stat.x(end));
	 	bb=repmat(stat.b,3,1);
		for jdx=1:2;
			barb(jdx) = mean(stat.b(:,jdx,id));
			plot(xx*stat.fc(jdx,id)-dx(jdx),bb(:,jdx,id)/barb(jdx),'linewidth',1);
			hold on;
		end;
		 xlim([0,10])
		lh=legend(num2str(cvec(s_a)));
		title(lh,'s_a');
		xlabel('Position $x/\lambda_c$','interpreter','latex');
		%ylabel('Biomass $b$ / (g/m$^2$)','interpreter','latex');
		ylabel('Biomass $b / \bar b$','interpreter','latex');
	
		% fit the curve
		fdx = (x*fc < 10);
		A   = sqrt(fc*x);
		s   = A(fdx) \ dphi(fdx);
		%plot(x*fc,A*s);
		xi = 0.32 + (0:9)';
		q0 = normcdf([-1,1]);
		qi = xi + norminv(q0,0,s*sqrt(xi))/(2*pi);
		%vline(qi)
		y = 2.75*ones(size(qi));
		cm = colormap('lines');
		lh = line(qi',y','color',cm(2,:),'linewidth',2);
		lh(end+1) = plot(mean(qi,2),mean(y,2),'k.','markersize',6);
		for jdx=1:length(lh)
			lh(jdx).HandleVisibility = 'off';
		end
		hold off
		ylim([0,3.5]);
	end

	figure(1e3);
	clf;
%	plot(stat.x*stat.fc(1),(stat.phi(:,1)+10*stat.dphi)/(2*pi));
	xlim([0,20]);
	lh=legend(num2str(cvec(stat.s_a)),'location','southeast');
	title(lh,'s_a');
	xlabel('Distance x/\lambda_c');
	ylabel('Phase \phi/(2 \pi)');

if (pflag)
	pdfprint(1e4+2,'img/rk-phase-error.pdf',3.5);
	pdfprint(1e4+5,'img/rk-phase-error-pattern.pdf',3.5);
end

end % rk_1d_experiment_phase_noise_integration

