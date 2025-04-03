% Sun 12 Nov 21:48:28 CET 2023
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
%% generate irregular patterns with the spatially extended grazing model
%% of May1977
function out = grazing_model_experiment(meta)
	
	if (nargin()<1)
		meta = pattern_formation_metadata();
	end
	ps    = meta.plotscale;
	pflag = meta.pflag;
	fflag = pflag;
	fcmap = meta.fcmap;
	
	optimize = false;
	% spatial extent
	Lx = 1024;
	% spatial extent of plot
	Lp = Lx;
	% spatial resolution
	dx = 1;
	%dt = 0.125;
%	dt = 0.0625;


% 5 0.5
% 6 0.25
	% spectral resolution
	df = 1./Lx;
	% numper of grid points
	nx = Lx/dx;
	% end time
	T = 1e3;

	% model parameters, these have been calibrated to ensure the same
	% mean biomass irrespectively of the regularity
	p0   = 0.33;
	g    = 2;
	% mean carrying capacity of the May
	k    = 8;
	mu0lr = [1.9,2.6];

	% plot colors
	cm=[1 1 1;
	 1 0 0;
	 0 0 1;
	 0.5 0 0.5];
	
	% declare arrays
	tab      = table();
	cv_k = [0.01   0.02      0.05      0.1000    0.2000    0.3000    0.4000    0.5000    0.6,      0.7,      0.8,      0.9,      1.0];
	tab.cv_k(1:length(cv_k)) = cv_k;
	tab.mu_b0(:) = [2.5288 2.5297    2.5203    2.4724    2.3106    2.1415    2.0336    1.9897    1.9492    1.9447    1.9394    1.9745    2.0200];
	nsa      = length(tab.cv_k);

	% for each experiment with varying degrees of heterogeneity
	for adx=1:nsa
	    adx
	    for idx=2
	    
	    % determine initial biomass concentration
	    if (idx==3)
	    	c = vander_1d(mu0_a(adx,1:2),1) \ p_a(adx,1:2)';
	    	mu0_a(adx,idx) = (p0-c(1))/c(2)
	    %else
	    %	mu0_ = mu0_a(idx);
	    end
	    if (optimize)
	    	opt.tolfun = 0.02;
	    	mu0_a(adx) = fminbnd(@(mu0) abs(fobj(mu0,sdk_a(adx)) - p0),mu0lr(1),mu0lr(2),opt);
	    else
	    	mu0_a(adx) = tab.mu_b0(adx);
	    end
	    % run model, extract values, or just load the result if the script has run before
	    [prec,rad] = precompute(mu0_a(adx),k*tab.cv_k(adx));
	    tab.runtime(adx) = sum(prec.runtime);
	    tab.filter_order(adx) = prec.filter_order; %_par(2); % par_a(adx,:) = prec.par;
	    tab.p_coverage(adx) = prec.p_coverage;
	    tab.fl(adx) = prec.fl;
	    %fl_a(adx,2) = prec.fl_lp;
	    tab.r2_bt(adx) = prec.r2_bt;
	    tab.r2_Sr(adx) = prec.r2_Sr;
	    tab.p_disagreement(adx) = prec.p_disagreement;
	    %pd(adx)     = prec.pd;
	    %gg(adx)     = prec.g;
	    tab.coherence(adx) = prec.coherence;
	    tab.mu_b(adx) = prec.mu_b;	    

	    % pattern, unscaled
	    splitfigure([3,ceil(nsa/3)],[1,adx],fflag,'',1000);
	    cla();
	    imagesc(prec.x,prec.x,prec.b);
	    xlabel('$x$','interpreter','latex');
	    ylabel('$y$','interpreter','latex');
	    clim(quantile(prec.b(:),[0.1,0.9]));
	    axis square;
	    colormap(fcmap(256));
	    
	    % pattern scaled
	    splitfigure([3,ceil(nsa/3)],[2,adx],fflag,'',1000);
	    cla();
	    imagesc(prec.x*prec.fl_lp,prec.x*prec.fl_lp,prec.b);
	    xlabel('Position $x/\lambda_l$','interpreter','latex');
	    ylabel('Position $y/\lambda_l$','interpreter','latex');
	    axis(10*[0,1,0,1]);
	    axis square;
	    clim(quantile(prec.b(:),[0.1,0.9]));
	    colormap(fcmap(256));
	    
	    % value histogram
	    splitfigure([3,ceil(nsa/3)],[3,adx],fflag,'',1000);
	    cla();
	    %xh = linspace(0.1,quantile(p.b(:),0.999),nx);
	    xh = logspace(log10(quantile(prec.b(:),0.001)),log10(quantile(prec.b(:),0.999)),nx);
	    %histogram(p.b(:),xh,'Normalization','pdf','DisplayStyle','stairs');
	    histogram(log(prec.b(:)),log(xh),'Normalization','pdf','DisplayStyle','stairs');
	    %vline(g);
	    xlabel('Biomass $\ln(b)$','interpreter','latex');
	    ylabel('Probability Density $P(\ln(b))$','interpreter','latex');
	    ylim([0,0.9]);
	    xhc = mid(xh);
	    xlim([log(xhc(1)),log(xhc(end))]+0.1*[-1,+1])
	    rad_ = Grazing_May1977();
	    rad.p = rad.pmu;
	    b0 = rad.homogeneous_states(0);
	    vline(log(b0));
	    vline(log(prec.g),'linestyle','-','color','b');
	    %log(xhc(1))
	    %pause
	    
	    % periodogram 2d
	    splitfigure([3,ceil(nsa/3)],[4,adx],fflag,'',1000);
	    cla();
	    imagesc(fftshift(prec.fx)/prec.fl_lp,fftshift(prec.fy)/prec.fl_lp,fftshift(prec.hatS)*prec.fl_lp^2);
	    s = 4;
	    xlim([-1,1]*s);
	    ylim([-1,1]*s);
	    S2d = trifilt2(fftshift(prec.hatS),3);
	    caxis([0,0.5*max(S2d*prec.fl_lp^2,[],'all')]);
	    axis square
	    xlabel('Wavenumber $k_x/k_l$','interpreter','latex');
	    ylabel('Wavenumber $k_y/k_l$','interpreter','latex');
	    cbh=colorbar();
	    title(cbh,'$\hat S_{xy}/\lambda_l^2$','interpreter','latex');
	    colormap(fcmap(256));
	    
	    % correlogram 2d
	    splitfigure([3,ceil(nsa/3)],[5,adx],fflag,'',1000);
	    cla();
	    x = real_axis_1d(Lx,nx); 
	    imagesc(fftshift(x)*prec.fl_lp,fftshift(x)*prec.fl_lp,fftshift(prec.Rhat));
	    s = 1.0;
	    xlim([-1,1]*s);
	    ylim([-1,1]*s);
	    %S2d = trifilt2(fftshift(prec.hatS),3);
	    %caxis([0,0.5*max(S2d,[],'all')]);
	    %caxis([-0.5,1]);
	    axis square
	    xlabel('Lag distance $x/\lambda_l$','interpreter','latex');
	    ylabel('Lag distance $y/\lambda_l$','interpreter','latex');
	    cbh=colorbar();
	    title(cbh,'$\hat R$','interpreter','latex');
	    caxis([-0.25,1])
	    colormap(fcmap(10));
	    
	    % radial density 
	    splitfigure([3,ceil(nsa/3)],[6,adx],fflag,'',1000);
	    yyaxis left
	    cla
	    plot(prec.fr/prec.fl_lp,prec.fl_lp*prec.Sr.normalized,'.','linewidth',1);
	    hold on;
	    plot(prec.fr/prec.fl_lp,prec.fl_lp*prec.Srlp,'-','linewidth',1);
	    s=4.0;
	    xlim([0,s]);
	    xlabel('Radial wavenumber $k_r/k_l$','interpreter','latex');
	    ylabel('Radial density $S_r/\lambda_l$','interpreter','latex');
	    hold off;
	    yyaxis right
	    cla;
	    plot(prec.fr/prec.fl_lp,prec.S.coherence.radial,'linewidth',1);
	    ylim([0,1]);
	    ylabel('Spectral coherence $c$','interpreter','latex');
	    legend('GM-Model','LP-fit','Coherence');
	    axis square;
	    
	    % angular density
	    splitfigure([3,ceil(nsa/3)],[7,adx],fflag,'',1000);
	    cla();
	    plot(prec.theta,prec.St,'linewidth',1)
	    hold on;
	    xlabel('Angle $\theta$','interpreter','latex');
	    ylabel('Angular density $S_\theta/\lambda_l$');
	    xlim([-1,1]*pi/2);
	    ylim([0,0.5]);
	    axis square
	    hline(1/pi,'color','r','linewidth',1); 
	    
	    % thresholded pattern compared to lowpass
	    % unscaled
	    splitfigure([3,ceil(nsa/3)],[8,adx],fflag,'',1000);
	    cla();
	    imagesc(prec.x,prec.x,prec.b_overlay);
	    axis(Lp*[0,1,0,1]);
	    axis equal;
	    xlabel('$x$','interpreter','latex');
	    ylabel('$y$','interpreter','latex');
	    colormap(cm);
	    axis square
	    
	    % thresholded pattern, scaled
	    splitfigure([3,ceil(nsa/3)],[9,adx],fflag,'',1000);
	    cla();
	    imagesc(prec.x*prec.fl_lp,prec.x*prec.fl_lp,prec.b_overlay);
	    axis equal;
	    xlim([0,10]);
	    ylim([0,10]);
	    xlabel('Position $x/\lambda_l$','interpreter','latex');
	    ylabel('Position $y/\lambda_l$','interpreter','latex');
	    colormap(cm);
	    axis square
	    
	    % mean and sd over time
	    splitfigure([3,ceil(nsa/3)],[10,adx],fflag,'',1000);
	    yyaxis left
	    cla();
	    plot(prec.t,[prec.mu,prec.sd]);
	    yyaxis right
	    cla();
	    plot(mid(prec.t),prec.d);
	    xlabel('time')
	    legend('\bar b','std(b)','1/||b||db/dt','dbarb/dt','||b(t)-b(T)||/||b(T)||');
	    end % for idx
	end % for adx
	
	splitfigure([2,3],[8,1],fflag);
	%plot(tab.cv_k,par_a(:,2),'.-','linewidth',1);
	plot(tab.cv_k,tab.filter_order,'.-','linewidth',1);
	xlim([0,1.05*max(tab.cv_k)]); 
	ylim([min(tab.filter_order)*0.9,1.1*max(tab.filter_order)]); 
	xlabel('Exogenous heterogeneity $cv(k)/k$','interpreter','latex');
	ylabel('Filter order','interpreter','latex');
	axis square;
	
	splitfigure([2,3],[8,2],fflag);
	plot(tab.cv_k,1./tab.fl,'linewidth',1);
	ylabel('Cutoff Wavelength $\lambda_l$','interpreter','latex');
	xlabel('Heterogeneity $cv(k)/k$','interpreter','latex');
	axis square;
	
	splitfigure([2,3],[8,3],fflag);
	plot(tab.cv_k,mu0_a,'linewidth',1);
	ylabel('Initial condition $b(0)$','interpreter','latex');
	xlabel('Heterogeneity $cv(k)/k$','interpreter','latex');
	axis square;
	
	splitfigure([2,4],[8,5],fflag);
	cla
	plot(tab.cv_k,tab.p_coverage,'linewidth',1);
	ylim([0,1]);
	ylabel('$Fraction of ground covered by vegetation$','interpreter','latex');
	xlabel('Heterogeneity $cv(k)/k$','interpreter','latex');
	axis square;
	
	splitfigure([2,4],[8,6],fflag);
	yyaxis left
	cla();
	markersize = 4;
	plot(tab.cv_k,tab.r2_bt,'^','markerfacecolor','k','color','none','markersize',markersize);
	hold on
	plot(tab.cv_k,tab.r2_Sr,'v','markerfacecolor','k','color','none','markersize',markersize);
	xlabel('Heterogeneity $CV(K)$','interpreter','latex');
	ylabel('Goodness of fit $R^2$','interpreter','latex');
	ylim([0,1.05]);
	yyaxis right
	cla();
	plot(tab.cv_k,tab.coherence,'o','color','none','markerfacecolor','r','markersize',markersize);
	ylabel('Spectral coherence $\bar c$','interpreter','latex');
	legend('Pattern $b$','Density $S_r$','Coherence $\bar c$','interpreter','latex','location','southeast');
	ylim([0,1.05]);
	axis square;
	if (0)
		yyaxis right;
		cla();
		plot(tab.cv_k,100*(1-pd));
		ylabel('Agreement / %');
	end
	
	splitfigure([2,4],[8,7],fflag);
	plot(tab.cv_k,tab.runtime)
	ylabel('Runtime / seconds');
	xlabel('Heterogeneity cv(k)');
	axis square;
	
	splitfigure([2,4],[8,8],fflag);
	plot(tab.cv_k,tab.p_disagreement)
	ylabel('Area of disagreement');
	xlabel('Heterogeneity cv(k)');
	axis square;

	%splitfigure([2,2],[10,1],fflag);
	%plot(tab.cv_k,coherence);

	tabw          = [tab(:,'cv_k') ...
			 , tab(:,'mu_b0') ...
			 , tab(:,'mu_b') ...
			 , tab(:,'p_coverage') ...
			 , tab(:,'fl') ...
			 , tab(:,'filter_order') ...
			 , tab(:,'coherence') ...
			 , tab(:,'r2_Sr') ...
			 , tab(:,'r2_bt') ...
			];
	tabw = round(tabw,3,'significant');
	writetable(tabw,'output/grazing-model-experiment.csv');

	if (pflag)
		for adx=1:length(tab.cv_k);
			base = sprintf('img/may/may1977-cvk_%g-L_%g-dx-%g',cv_k(adx),Lx,dx);
			pdfprint(1000+adx,[base,'-pattern.pdf'],ps);
			pdfprint(2000+adx,[base,'-pattern-scaled.pdf'],ps);
			pdfprint(3000+adx,[base,'-probability-density.pdf'],ps);
			pdfprint(4000+adx,[base,'-periodogram-2d.pdf'],ps);
			pdfprint(5000+adx,[base,'-correlogram-2d.pdf'],ps);
			pdfprint(6000+adx,[base,'-density-radial.pdf'],ps);
			pdfprint(7000+adx,[base,'-density-angular.pdf'],ps);
			pdfprint(8000+adx,[base,'-pattern-rad-model-vs-lowpass.pdf'],ps);
			pdfprint(9000+adx,[base,'-pattern-rad-model-vs-lowpass-scaled.pdf'],ps);
			
		end
		base = sprintf('img/may/may1977-L_%g-dx-%g',Lx,dx);
		pdfprint(80+1,[base,'-filter-order-vs-cvk.pdf'],ps);
		pdfprint(80+2,[base,'-lc-vs-sk.pdf'],ps);
		pdfprint(80+3,[base,'-b0-vs-sk.pdf'],ps);
		pdfprint(80+5,[base,'-p-dense-vs-cvk.pdf'],ps);
		pdfprint(80+6,[base,'-goodness-of-fit-vs-cvk.pdf'],ps);
		pdfprint(80+7,[base,'-runtime-vs-cvk.pdf'],ps);
		pdfprint(80+8,[base,'-area-fraction-of-disagreement-vs-cvk.pdf'],ps);
	end
	
	function val = fobj(mu_b0,sd_k)
		prec = precompute(mu_b0,sd_k);
		%,prec.p_coverage]);
		val = prec.p_coverage;
	end
	
	function [prec,rad] = precompute(mu_b0,sd_k)
		disp([mu_b0,sd_k]); 
	
		param         = struct();
		% param.opt.dt  = dt; % 1; % 2
		param.opt.adapt_time_step=1;
		param.opt.dt           = 1/400;
		param.opt.dt_min       = 1/400;
		param.opt.dt_max       = 0.5;
		param.opt.outer_abstol = 1e-4;
		param.opt.outer_reltol = 1e-2;
		param.opt.dt_max_scale_up   = sqrt(2);
		param.opt.dt_min_scale_down = 0;
		param.opt.rms_delta_zo_rel_max = 0.2;

		param.opt.dto = inf;
		param.T       = T;
		param.opt.path_str = 'mat/may/';
		param.nx      = nx*[1,1];
		param.L       = Lx*[1,1];
		
		param.initial_condition.mu = mu_b0;
		param.initial_condition.sd = 0;
		param.initial_condition.dist = {'normal'};
		param.pmu.k       = k;
		param.pss.k       = sd_k; 
		%param.psdist.k    = 'gamma';
		param.psdist.k = 'geometric-ornstein-uhlenbeck';
		%param.psl.k       = [min(256,param.L(1)/4)];
		param.psl.k       = 1;

		param.opt.solver  = 'solve_split';
		param.opt.inner_solver = 'step_advect_diffuse_spectral';
		%param.opt.solver  = @ode23;
		
	        rad    = Grazing_May1977(param);
		f_str0 = rad.filename();
		f_str  = [f_str0(1:end-4),'-analyzed.mat']
		
		if (exist(f_str,'file'))
			load(f_str);
			load(f_str0,'rad','out');
		else
			[t,bb,out]=rad.run();
			prec.runtime = out.runtime;
			bb  = single(bb);
			b   = rad.extract2(bb(end,:));
			ae  = reshape(rad.p.k,nx,nx);
			prec.x = rad.x;

			% TODO make this is afunction analyze anisotropic in spatial_pattern
			sp = Spatial_Pattern
			sp.L = rad.L
			sp.b = b;
			sp.source = ae;
			sp.opt.suppress_low_frequency_components = 0;
			sp.analyze_grid();
			
			%tab.mu(adx,idx) = mean(b,'all');
			prec.mu_b = mean(b,'all');
			prec.sd_b = std(b,[],'all');
			prec.me_b = median(b,'all');
			g = graythresh_scaled(b);
			%bmin = min(b(:));
			%bmax = max(b(:));
			%g = bmin+graythresh((b-bmin)/(bmax-bmin))*(bmax-bmin);
			
			sd = std(bb,[],2);
			mu = mean(bb,2);
			
			hatS = abs(fft2(b-mean(b,'all'))).^2;
	    		hatS = hatS/(sum(hatS,'all')*df*df);
			Rhat = ifft2(hatS);
			Rhat = Rhat/Rhat(1);
			
			[Sr,fr] = periodogram_radial(hatS,[Lx,Lx]);
			w = fr;
			%try
			[par,Srlp,fitstat]=fit_spectral_density(fr,Sr.normalized,fr,@lowpass1dpdf,[1.5,15],'hellinger',0);
			prec.filter_fc = par(1);
			prec.filter_order = par(2);
			
			if (0)
				Cr=cumsum((cvec(fr).*cvec(Sr.normalized)));
				fl=interp1(Cr/Cr(end),fr,0.5);
			else
				% note, this is not a good estimate of fl,
				% as the value for the mean is set to 0
				% setting it Sr(1) to max(Sr) fixes this somewhat
				%Sr_ = Sr.normalized;
				%Sr_(1) = max(Sr_);
				Ss = sort(Sr.normalized,'descend');
				fl = interp1(Ss,fr,0.5*Ss(1));
				Ss = sort(Sr.normalized,'descend');
				fl_lp = interp1(Srlp,fr,0.5*Srlp(1));
				%fl_a(adx,idx) = fl_lp;
			end
			% for plot
			Sr.normalized(1) = NaN;
			
			bb = single(bb);
			%fl = fl_lp;
			
			[St,theta] = periodogram_angular(hatS,[Lx,Lx]);
			
			[fx,fy,frr] = fourier_axis_2d([Lx,Lx],[nx,nx]);
			%Srlp = lowpass2d_pdf(fr,par(1),par(2),1);
			Slp2d = interp1(fr,Srlp,frr,'linear',0);
			blp = ifft2(sqrt(Slp2d).*fft2(reshape(ae,nx,nx)));
			me_b   = median(b,'all');
			me_blp = median(blp,'all');
			p_coverage     = mean(b>g,'all');
			q_     = quantile(blp(:),1-p_coverage);
			b_thresh = (b>g);
			blp_thresh = (blp>q_);
			b_overlay   = b_thresh + 2*(blp_thresh);
			d      = rms(diff(bb),2)./diff(t)./mean(bb(:));
			d(:,2) = diff(mean(bb,2))'./diff(t)./mean(bb(:));
			d(:,3) = mid(rms(bb - bb(end,:),2))/rms(bb(end,:));
			
			prec.d = d;
			prec.p_coverage = p_coverage;
			prec.t = t;
			prec.fl = fl;
			prec.fl_lp = fl_lp;
			prec.fx = fx;
			prec.fy = fy;
			prec.fr = fr;
			prec.b = b;
			prec.fx = fx;
			prec.b_overlay = b_overlay;
			prec.b_thresh = b_thresh;
			prec.blp_thresh = blp_thresh;
			prec.sd = sd;
			prec.mu = mu;
			prec.Sr = Sr;
			prec.Srlp = Srlp;
			prec.hatS = hatS;
			prec.Rhat = Rhat;
			prec.theta = theta;
			prec.St = St;
			prec.p_disagreement = mean(b_overlay == 1,'all') + mean(b_overlay == 2,'all');
			prec.g = g;
		
			prec.r2_b = corr(flat(blp),flat(b));
			% corr_a(adx,2) = corr(flat(Sr.normalized),flat(Srlp))
			%prec.corr(2) = wcorr(cvec(fr),sqrt(flat(Sr.normalized)),cvec(fr),sqrt(flat(Srlp))); 
			Sr.normalized(1) = 0;
			c = corr(flat(blp_thresh),flat(b_thresh));
			prec.r2_bt = sign_to_pearson(c.^2);

			prec.r2_Sr = fitstat.goodness.r2;
			%df = 1./param.L(1);
			% hellinger_distance(flat(Sr.normalized),flat(Srlp),df,fr);
			prec.coherence = sp.stat.coherence.radial;
			prec.S.coherence.radial = sp.S.coherence.radial;

			save(f_str,'prec','sp');
		end % if not yet precomputed
	end % function precompute
end % function grazing_model_experiment

