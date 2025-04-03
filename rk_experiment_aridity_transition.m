% Sat  9 Nov 15:59:40 CET 2024

	if (~exist('pflag','var'))
		pflag = 0;
	end
	fflag = pflag;

	rng(0);

	s = 1;

	dim = [1,1];
%	dim = 1;

	% length of domain
	if (length(dim)>1)
		L = 2*128*dim;
	else
		L = 1000;
	end

	% spatial resolution
	dx = 1*dim;

	% number of grid points
	nx = round(L./dx);

	a0 = 0.2;
	% relative magnitude of spatial variation of the bare soild infiltration 
	%s_a = 0.000;

	%a = a0;
	%a = a0+a0*s_a*randn(nx);

	% final time
	To  = s*1e4;

	% output time step
	dto = s*100;
 
	% time step
	dt  = 0.5;

	% model parameters
	param        = struct();

	% advection
	param.pmu.vx = [0,0,0];

	% diffusion
	param.pmu.ex = [0.1,0.1,100];

	% precipitation
	%param.pmu.R  = 1;
	Rfun = sprintf('@(t) 1.15*(1 - t/%g)+0.35',To)                                           
	param.pmu.R = eval(Rfun);
	%t0 = 0.2*To;
	%param.p.rw = @(t) 0.1*(1 + 4*max(0,(t-t0)/(To-t0)));
	%param.p.db = @(t) 2/3*0.25*(1 + t/To);


	% boundary condition
	param.boundary_condition    = {'circular','ciruclar'};

	a0 = 0.2;
	sa = 0.1;

	% reload values of intermediate time steps
	param.opt.loadfinal = false;
	param.L  = L;
	param.nx = nx;
	param.T  = To;
	param.opt.dt  = dt;
	param.opt.dto = dto;
	param.pss.a   = sa*a0;
	param.psdist.a = 'lognormal';
	param.pstrel.b  = 0;
	param.pstrel.w  = 0;
	param.pstrel.h  = 0;
	param.opt.compute_class = @single;
	param.opt.path_str = 'mat/';
	%param.pmu.a = 0.2;
	param.opt.solver = 'solve_split';
	param.opt.inner_solver = 'step_advect_diffuse_spectral';
	param.opt.adapt_time_step = false;

	% initial condition
	rk = Rietkerk(param);
	[b0,w0,h0] = rk.homogeneous_state(0,rk.pmu);
	ic = [rand(prod(nx),1);
             w0*ones(prod(nx),1);
	      mean(h0)*ones(prod(nx),1)];
	param.initial_condition = ic;
	%param.initial_condition = 'obj.random_state()';

	% run model
	rk = Rietkerk(param);
	rk.hashfield_C{end+1} = 'pstrel';
	[t,y]	  = rk.run();
	y = single(y);
	[bb,ww,hh] = rk.extract1(y);
	
	splitfigure([2,3],[1,1],fflag)
	imagesc(bb)

	p = [];
	thresh = 5;
	p = mean(bb>thresh,2);

	ystat.mean = mean(bb,2);
	ystat.std = std(bb,[],2);
	ystat.sk = skewness(bb,[],2);
	ystat.ku = kurtosis(bb,[],2);
	R = param.pmu.R(t);
	b = y(:,1:end/3);
	for idx=1:length(t)
		%yt(idx,1) = graythresh_scaled(y(idx,:));
		bthresh(idx,1) = graythresh(b(idx,:));
		bhigh(idx,1) = mean(b(idx,b(idx,:)>bthresh(idx,1)));
		blow(idx,1)  = mean(b(idx,b(idx,:)<bthresh(idx,1)));
		%yp(idx,1) = mean(b(idx,:)>yt(idx));
		%ystat.p(idx,1) = mean(b(idx,:)>0.5);
		ystat.p(idx,1) = mean(b(idx,:)>bthresh(idx));
	end
	ypm = flipud(make_monotonic(flipud(ystat.p),1));
	tdx = interp1(ypm,1:length(t),0.25,'nearest');
	b = rk.extract2(y(tdx,:));%round((0.25*tdx1+0.75*tdx0)),:));

	t1 = interp1(ypm,t,0.999);
%	tdx1 = find(yp<=0.99 & cvec(t) > 100,1,'first')
%	t1 = t(tdx1);
	tdx0 = find(ystat.p<=0.01 & cvec(t) > 100,1,'first')
	t0 = t(tdx0);
	R1 = interp1(t,R,t1)
	R0 = interp1(t,R,t0)
	%bmin = 0;
	%b1   = interp1(t,ystat.mean,t1);
	fdx  = cvec(t >= t1 & t<=t0);
	for idx=1:2
	if (1 == idx)
		yyaxis left
		rhs = ystat.mean;
	else
		yyaxis right
		rhs = ystat.p;
	end
	y1(idx) = interp1(t,rhs,t1)
	slope = y1(idx)/(R1-R0)
	A    = vander_1d(R(fdx),1);
	cb   = A \ rhs(fdx)
	Rmin = roots([cb(2),cb(1)]);
	yp_ = A*cb;
	yp_(:,2) = y1(idx)*(R(fdx)-R0)/(R1-R0);
	rmse = rms(rhs(fdx)-yp_)
	R2 = 1 - rmse.^2/var(rhs(fdx))	

	ypp(:,idx) = yp_(:,2);
%	hold on
%	if (1 == idx)
%	plot(R(fdx),yp_(:,2));
%	else
%		plot(R(fdx),yp_(:,2))
%		%(R(fdx)-R0)/(R1-R0));
%	end

	end
	printf('Rmax %g\n',R1);
	printf('Rmin %g\n',R0);
	printf('bmax %g\n',y1);


	splitfigure([2,3],[1,1],fflag);
	imagesc(bb);

	splitfigure([2,3],[1,2],fflag);
	imagesc(b);
	axis equal;
	axis tight;

	splitfigure([2,3],[1,3],fflag);
	plot(t,[ystat.mean,ystat.std,ystat.sk,ystat.ku,p*10]);
	legend('mean','std','sk','ku','p');

	splitfigure([2,3],[1,5],fflag)
	yyaxis left
	plot(R,ystat.mean,'linewidth',1);
	hold on
	plot(R(fdx),ypp(:,1),'k--');
%	plot(R,[bhigh,ystat.mean./ystat.p])
	yyaxis right
	plot(R,ystat.p,'linewidth',1);
	hold on
	plot(R(fdx),ypp(:,2),'r--');
	yyaxis left
	xlim([0.3,1.4]);
	ylim([0,16]);
	ylabel('Mean biomass $\bar b$ / (g/m$^2$)','interpreter','latex');
	xlabel('Precipitation rate $R$ / (mm/day)','interpreter','latex')
	vline(R0,'color','k','linestyle','--','linewidth',1);	
	vline(R1,'color','k','linestyle','--','linewidth',1);
	yyaxis right;
	ylabel('Vegetation cover $p$','interpreter','latex');
	ylim([0.00,1.1]);
	axis square

	splitfigure([2,3],[1,6],fflag);
	plot(ystat.p(t>100),ystat.mean(t>100));

	splitfigure([2,3],[1,4],fflag);
	plot(R,[blow,bthresh,bhigh])

	if (pflag)
		ps = 3.5;
		pdfprint(15,'img/rietkerk-biomass-vs-precipitation.pdf',ps);
	end

