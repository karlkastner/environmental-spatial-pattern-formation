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
function [stat, bb, c] = rk_1d_pni(meta)
	if (nargin<1)
		meta = pattern_formation_metadata();
	end
	dflag = 0;
	pflag = meta.pflag;
	fflag = pflag;
	ps    = meta.plotscale;
	bold = [];
	scale = 2;
	% length of domain
	L   = scale*4000;
	x0  = 0;
	% spatial resolution
	dx  = scale*2;
	% time step
	dt  = scale;
	% number of grid points
	nx = round(L./dx);

	a0 = 0.2;
	% relative magnitude of spatial variation of the bare soild infiltration 
%	s_a   = [1e-4,2e-4,5e-4,0.001:0.001:0.01,0.02,0.05,0.1]'*a0;
%	s_a = [1e-2:1e-2:1e-1];
%	s_a  = [0;logspace(-3,-1,11)'];
	%seed = round(s_a*1e4); %+1e4*(0:3);
	%s_a = [0, 0.01:0.01:0.1]';
	%s_a = [0, 0.1*ones(1,10)]';
	s_a = [0,0.001,0.01,0.1];
	s_a = 0.01;
%	seed(:) = 25;
%	seed(1) = 0;
	seed = ones(length(s_a),1)*(1:10);
%	seed   = round(s_a/a0*1000)+[0,1e3,2e3];
%	seed(1,1:3) = 1e4+[0,1e3,2e3];
%	seed(2,1:3) = 2e4+[0,1e3,2e3];
%	seed(3,1:3) = 3e4+[0,1e3,2e3];

	%% window size
	%nwin = 2*round(800*sqrt(scale))+1;
	%nwh = (nwin-1)/2;

	% final time
	To  = scale*3.2e4;
	% output time step
	dto = 100; 

	% frequencies at which pattern is plotted
	%f0 = [0.5,1,2];

	% local frequency
	%k = 20;
	%r = 1-1/8;
	%f  = 12*r.^(0:k/r)';
	%f = f*2;
	%F_ = scale*[0;cumsum(f)];
	%F = mid(F_);

	% phi = omega*x <-> omega = unwrapp(phi) / x 
	
	% interpolate to grid
	%n  = 0.5*nx+1;
	%n_ = length(F);
	%F  = interp1((1:n_)/(n_+1),F,(0:n-1)/(n-1),'spline')';
	%df_dx = cdiff(F)./dx;
	% local bare soil infiltration
	%ea  = sin(2*pi*(F-F(1)));

	% when we only increase the frequency, then there will be a jump at the end of the domain
	% so we ramp it up and down to make it symmetric
	%ea = [ea;-flipud(ea(2:end-1))];
	%a = a0*(1 + s_a*ea);  
	%df_dx_ = [df_dx;-flipud(df_dx(2:end-1,:))];

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
	param.opt.loadfinal = false;
	param.L = L;
	param.nx = nx;
	param.T = To;
	param.opt.dt = dt;
	param.opt.dto = dto;
	param.pss.a  = 0;
	param.opt.compute_class = @single;
	param.opt.path_str = 'mat/';
	param.pmu.a = a0;
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_implicit_q_fft';
%	param.opt.inner_solver = 'step_advect_diffuse_trapezoidal';
%	param.opt.inner_solver = 'step_advect_diffuse_spectral';
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
	% 'obj.initial_condition_periodic()';
	corr_ = NaN(length(s_a),1);
	corr__ = NaN(length(s_a),1);
	for sdx=1:size(seed,2)
	for adx=1:length(s_a)
	rng(seed(adx,sdx));

	% vayring half
	sd_a = a0*s_a(adx);

	ne = sum(rk.x>=x0);
if (0)
	[lmu,lsd] = logn_moment2par(a0,sd_a);
	a = lognrnd(lmu,lsd,ne,1);
else
	[p1,p2] = uniform_moment2par(a0,sd_a);
	a = uniformrnd(p1,p2,ne,1);
end
	% ensure mean and sd
	a = a0 + (sd_a/std(a))*(a-mean(a));
	% append non-varying half
	a = [repmat(a0,nx-ne,1);a];
	% make a double bridge
	%a(nx/2+1:end) = (a(nx/2+1:end) - mean(a(nx/2+1:end))+a0);
	param.p.a = a;



	% initial condition
%	rk = Rietkerk(param);
%	[b0,w0,h0] = rk.homogeneous_state(rk.pmu);
%	rng(0);
%	ic = [rand(nx,1);
%             w0*ones(nx,1);
%	      mean(h0)*ones(nx,1)];
	%param.initial_condition = 'obj.random_state()';
	%param.initial_condition = ic;

	% run model
	rk = Rietkerk(param);
	rk.hashfield_C{end+1} = 'p.a';
	[t,y]	  = rk.run(); %param);
	y = single(y);
	[bb,ww,hh] = rk.extract1(y);
	size(bb)
	b = bb(end,:)';
	%rms(bb)	

	fx = fourier_axis(rk.x);
	S = abs(fft2(b-mean(b)));
	[Sc,mdx] = max(S)
	%fc_(adx,sdx) = fx(mdx);
	[Sc_, fc(adx,sdx)] = extreme3(abs(fx),S,mdx);
%pause
	%fc(adx,1) = abs(fx(mdx))

	%ax = plotyy(rk.x,[bb(end,:)'/(2*rms(bb(end,:))),0.5*(1+(rk.pmu.a-a0)/(s_a*a0))],rk.x,ff)
	
%	% periodically extend, as the first and last window exceed the range
%	b = [b;b;b];
%	a = [a;a;a];
%
%	% goodness of fit (correlation)
%	c=NaN(rk.nx,2);
%	w=triwin(1:nwin)'; % triwin
%	wh=hanwin(1:nwin)'; % triwin
%	for idx=1:rk.nx
%		% rectwin
%		c(idx,1) = corr(b(rk.nx+idx+(-nwh:nwh)),a(rk.nx+idx+(-nwh:nwh)));
%		% triwin
%		c(idx,2) = wcorr(w,b(rk.nx+idx+(-nwh:nwh)),w,a(rk.nx+idx+(-nwh:nwh)));
%		c(idx,3) = wcorr(wh,b(rk.nx+idx+(-nwh:nwh)),wh,a(rk.nx+idx+(-nwh:nwh)));
%	end
%	nbartlett = round(sqrt(rk.nx));
%	S = periodogram_bartlett(bb(end,:)-mean(bb(end,:)), rk.L,nbartlett,rk.nx);
%	fx = fourier_axis(rk.x);
%	[Sc,mdx] = max(S);
%	fc = fx(mdx);
%	%Sa = periodogram_bartlett(rk.pmu.a-mean(rk.pmu.a), rk.L, 40*scale,rk.nx);
%	%[Sca,mdxa] = max(Sa);
%	%fca = fx(mdxa);

	% plot pattern and frequency
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
	xl = x(1:nl);
	s = conv(b3,sin(2*pi*fc(adx,sdx)*xl),'same');
	c = conv(b3,cos(2*pi*fc(adx,sdx)*xl),'same');
	s=s(nx+1:2*nx);
	c=c(nx+1:2*nx);
	phi = atan2(s,c);
	phi = unwrap(phi);
	phi = phi-phi(1); 
	fdx = find((x>3300) & (x<3700));

	%fc(adx) = 1/(2*pi)*(phi(end)-phi(1))/(x(end)-x(1))
	end
	%dphidx = (phi(fdx(end))-phi(fdx(1)))/(x(fdx(end))-x(fdx(1)));
	%dphi = phi - dphidx*x;
	if (adx == 1)
		phi0 = phi;
	%	dphi = 0.*phi;
	%else
	%	dphi = phi-phi0;
	end
	dphi = phi - 2*pi*fc(adx,sdx)*x;
	%ddphi(:,adx) = dphi;

	corr_(adx,sdx) = corr(dphi,ca);
	corr__(adx,sdx) = corr(dphi,ca_);
	A = [ones(nx,1),ca_];
	cc(:,adx) = A \ dphi;
	dphi_ = A*cc(:,adx);
	r2(adx,sdx) = 1 - mean((dphi_-dphi).^2)/var(dphi);
fc
r2

	if (1) %sdx == 1)
	%figure(adx*10);
	%clf
	splitfigure([2,3],[1e3*adx+sdx,1],fflag);
	%yyaxis left
	cla
	imagesc(bb)
%	plot(rk.x,[bb(end,:)'/(2*rms(bb(end,:)))]); %,0.5*(1+ea)]);
	%plot(rk.x,b);
	%set(gca,'colororder',{'k','r'})
	%set(gca,'colororder',[0,0,0;1 0 0])  
	%yyaxis right
	%cla
	%plot(rk.x,df_dx_);
	%hline(fc)
	splitfigure([2,3],[1e3*adx+sdx,2],fflag);
	cla
		%plot(rk.x,bb(1,:));
	%hold on
	plot(fc(adx)*(rk.x-x0),b);
	if(1==adx)
		b0 = b;
	end
	hold on
	plot(fc(adx)*(rk.x-x0),[b0]);
	%plot((rk.x-x0)/l,[b0,bold]);
	%bb(round(end-10),:));
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
	xlabel('Distance $x/\lambda_c$','interpreter','latex');
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
	end
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
		end
	end % for idx
	% interpolate patch number to x
	np = interp1(xc,1:length(xc),x,'linear','extrap');
end
	dnp_dx = (np(fdx(end))-np(fdx(1)))/(x(fdx(end))-x(fdx(1)));
	dnp_dx = (np(end)-np(end/2))/(x(end)-x(end/2));
	dnp = np-dnp_dx*x;
	if (1 == adx)
	np0 = np;
	else
		dnp = np - np0; 
	end
	corr__(adx,sdx) = corr(dnp,ca);
	%[cvec(s_a), corr_(:,1),corr__]
	nanmedian([cvec(s_a), corr_(:,1),corr__])
	
% phase shift with respect to the s_a = 0 pattern
	%for idx=1:nx-nl+1
	%	fdx = idx+(1:nl);
	%	xc = xcorr(b0(fdx),b(fdx))
	%end
	

%	corr_(adx,2) = corr(a_(nx/2+1:end),ca_(nx/2+1:end))
	% tan = sin/cos
 %x=rk.x';
 %n=round(l/dx(1));
 %m=floor(length(b)/n);
 %x=x(1:n);
 %b = reshape(b(1:n*m),n,m);
 %c=[ones(n,1),sin(2*pi*x/100),cos(2*pi*x/100)]\b;

	splitfigure([2,3],[1e3*adx+sdx,4],fflag);
	plot(rk.x,[phi,2*pi*fc(adx)*rk.x']/(2*pi));
%ww(end,:))

	splitfigure([2,3],[1e3*adx+sdx,5],fflag);
	%plot(rk.x,hh(end,:))
	plot(fx,S);
	vline(fc(adx),'color',[1,0,0])
	%vline(fc_(adx),'color','b')
	%rk.x,rk.p.a/a0)

	splitfigure([2,3],[1e3*adx+sdx,6],fflag);
	plot(rk.x,dnp);
	yyaxis right
	plot(rk.x,ca);
	%rk.p.a)
	end% if sdx ==1
	bold = b;
end
end

	figure(1e4);
	subplot(2,2,1);
	a = atanh(corr_);
 	ma = mean(a,2);
	sa=std(a,[],2);
	a=tanh(ma+sa*[-1,0,1]);
	%errorbar(s_a/a0,a(:,2),a(:,2)-a(:,1),a(:,3)-a(:,2));
	plot(s_a,a(:,2).^2);
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
	%stat.dphi = ddphi;
	stat.fc = fc;
	stat.s_a = s_a;


	figure(1e5);
 bar(median(stat.r2(2:end,:),2),'facecolor',[0,0,0.6]);
 ylabel('Goodness of fit R^2');
 xlabel('Degree of heterogeneity $s_a$','interpreter','latex');
 set(gca,'xticklabel',stat.s_a(2:end));
 axis square;
if (0) %pflag)
 pdfprint(1e5,'img/phase-prediction-r2.pdf',4);
end

end % rk_1d_pni


