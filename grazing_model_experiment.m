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
	Lp = 400;
	% spatial extent
	Lx = 400;
	% spatial resolution
	dx = 1;
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
	k    = 8;
	s_a = [0.05 0.1000    0.2000    0.3000    0.4000    0.5000 0.6, 0.7, 0.8, 0.9, 1.0];
	%s_a = [0.05 0.1000    0.2000    0.3000      0.5000 0.6, 0.7, 0.8, 0.9, 1.0];
	mu0lr = [1.9,2.6];
	%mu0o = [     2.5203    2.4724    2.3106    2.1415    2.0336    1.9897    1.9492 ];
	mu0o = [    2.5203    2.4724    2.3106    2.1415    2.0336    1.9897    1.9492    1.9447    1.9394    1.9745    2.0200];
	

	% plot colors
	cm=[1 1 1;
	 1 0 0;
	 0 0 1;
	 0.5 0 0.5];
	
	nsa = length(s_a);
	% declare arrays
	mu0_ = [];
	corr_a = [];
	fc_a = [];
	p_a = [];
	S0_a = [];
	par_a = [];
	sdk_a = k*s_a;
	par_a = [];
	rt_a = [];
	% for each experiment with varying degrees of heterogeneity
	for adx=1:length(sdk_a)
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
	    	mu0_a(adx) = mu0o(adx);
	    end

	    % run model, extract values, or just load the result if the script has run before
	    [p,rad] = precompute(mu0_a(adx),sdk_a(adx));
	    rt_a(adx) = sum(p.runtime);
	    par_a(adx,:) = p.par;
	    p_a(adx)    = p.pp;
	    fc_a(adx,1) = p.fc;
	    fc_a(adx,2) = p.fc_lp;
	    corr_a(adx,:) = p.corr;
	    r2_S(adx,1) = p.r2_S;
	    pd(adx)     = p.pd;
	    gg(adx)     = p.g;
	    
	    % pattern, unscaled
	    splitfigure([3,ceil(nsa/3)],[1,adx],fflag,'',1000);
	    cla();
	    imagesc(p.x,p.x,p.b);
	    xlabel('$x$','interpreter','latex');
	    ylabel('$y$','interpreter','latex');
	    clim(quantile(p.b(:),[0.1,0.9]));
	    axis square;
	    colormap(fcmap(256));
	    
	    % pattern scaled
	    splitfigure([3,ceil(nsa/3)],[2,adx],fflag,'',1000);
	    cla();
	    imagesc(p.x*p.fc_lp,p.x*p.fc_lp,p.b);
	    xlabel('Position $x/\lambda_l$','interpreter','latex');
	    ylabel('Position $y/\lambda_l$','interpreter','latex');
	    axis(10*[0,1,0,1]);
	    axis square;
	    clim(quantile(p.b(:),[0.1,0.9]));
	    colormap(fcmap(256));
	    
	    % value histogram
	    splitfigure([3,ceil(nsa/3)],[3,adx],fflag,'',1000);
	    cla();
	    %xh = linspace(0.1,quantile(p.b(:),0.999),nx);
	    xh = logspace(log10(quantile(p.b(:),0.001)),log10(quantile(p.b(:),0.999)),nx);
	    %histogram(p.b(:),xh,'Normalization','pdf','DisplayStyle','stairs');
	    histogram(log(p.b(:)),log(xh),'Normalization','pdf','DisplayStyle','stairs');
	    %vline(g);
	    xlabel('Biomass $\ln(b)$','interpreter','latex');
	    ylabel('Probability Density $P(\ln(b))$','interpreter','latex');
	    ylim([0,0.9]);
	    xhc = mid(xh);
	    xlim([log(xhc(1)),log(xhc(end))]+0.1*[-1,+1])
	    rad_ = May1977();
	    rad.p = rad.pmu;
	    b0 = rad.homogeneous_states(0);
	    vline(log(b0));
	    vline(log(p.g),'linestyle','-','color','b');
	    %log(xhc(1))
	    %pause
	    
	    % periodogram 2d
	    splitfigure([3,ceil(nsa/3)],[4,adx],fflag,'',1000);
	    cla();
	    Shat = p.Shat;
	    Shat = 2*Shat/(sum(Shat,'all')*df*df);
	    imagesc(fftshift(p.fx)/p.fc_lp,fftshift(p.fy)/p.fc_lp,fftshift(Shat)*p.fc_lp^2);
	    s = 4;
	    xlim([-1,1]*s);
	    ylim([-1,1]*s);
	    S2d = trifilt2(fftshift(Shat),3);
	    caxis([0,0.5*max(S2d*p.fc_lp^2,[],'all')]);
	    axis square
	    xlabel('Wavenumber $k_x/k_l$','interpreter','latex');
	    ylabel('Wavenumber $k_y/k_l$','interpreter','latex');
	    cbh=colorbar();
	    title(cbh,'$\hat S/\lambda_l^2$','interpreter','latex');
	    colormap(fcmap(256));
	    
	    % correlogram 2d
	    splitfigure([3,ceil(nsa/3)],[5,adx],fflag,'',1000);
	    cla();
	    x = real_axis_1d(Lx,nx); 
	    imagesc(fftshift(x)*p.fc_lp,fftshift(x)*p.fc_lp,fftshift(p.Rhat));
	    s = 1.0;
	    xlim([-1,1]*s);
	    ylim([-1,1]*s);
	    %S2d = trifilt2(fftshift(p.Shat),3);
	    %caxis([0,0.5*max(S2d,[],'all')]);
	    %caxis([-0.5,1]);
	    axis square
	    xlabel('Lag distance $x/\lambda_l$','interpreter','latex');
	    ylabel('Lag distance $y/\lambda_l$','interpreter','latex');
	    cbh=colorbar()
	    title(cbh,'$\hat R$','interpreter','latex');
	    caxis([-0.25,1])
	    colormap(fcmap(10));
	    
	    % radial density 
	    splitfigure([3,ceil(nsa/3)],[6,adx],fflag,'',1000);
	    cla;
	    plot(p.fr/p.fc_lp,p.fc_lp*p.Sr.normalized,'linewidth',1);
	    hold on;
	    plot(p.fr/p.fc_lp,p.fc_lp*p.Srlp,'linewidth',1);
	    xlabel('Radial wavenumber $k_r/k_l$','interpreter','latex');
	    ylabel('Radial density $S/\lambda_l$','interpreter','latex');
	    s=4.0;
	    xlim([0,s]);
	    legend('GM-Model','LP-fit');
	    axis square;
	    
	    % angular density
	    splitfigure([3,ceil(nsa/3)],[7,adx],fflag,'',1000);
	    cla();
	    plot(p.theta,p.St,'linewidth',1)
	    hold on;
	    xlabel('Angle $\theta$','interpreter','latex');
	    ylabel('Angular Density $S/\lambda_l$');
	    xlim([-1,1]*pi/2);
	    ylim([0,0.5]);
	    axis square
	    hline(1/pi,'color','r','linewidth',1); 
	    
	    % thresholded pattern compared to lowpass
	    % unscaled
	    splitfigure([3,ceil(nsa/3)],[8,adx],fflag,'',1000);
	    cla();
	    imagesc(p.x,p.x,p.b_overlay);
	    axis(Lp*[0,1,0,1]);
	    axis equal;
	    xlabel('$x$','interpreter','latex');
	    ylabel('$y$','interpreter','latex');
	    colormap(cm);
	    axis square
	    
	    % thresholded pattern, scaled
	    splitfigure([3,ceil(nsa/3)],[9,adx],fflag,'',1000);
	    cla();
	    imagesc(p.x*p.fc_lp,p.x*p.fc_lp,p.b_overlay);
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
	    plot(p.t,[p.mu,p.sd]);
	    yyaxis right
	    cla();
	    plot(mid(p.t),p.d);
	    xlabel('time')
	    legend('\bar b','std(b)','db/dt','||b(t)-b(T)||');
	    end % for idx
	end % for adx
	
	splitfigure([2,3],[8,1],fflag);
	plot(s_a,par_a(:,2),'.-','linewidth',1);
	xlim([0,1.05*max(s_a)]); 
	ylim([min(par_a(:,2))*0.9,1.1*max(par_a(:,2))]); 
	xlabel('Exogenous heterogeneity $s_k/k$','interpreter','latex');
	ylabel('Filter order','interpreter','latex');
	axis square;
	
	splitfigure([2,3],[8,2],fflag);
	plot(s_a,1./fc_a(:,2),'linewidth',1);
	ylabel('Cutoff Wavelength $\lambda_l$','interpreter','latex');
	xlabel('Heterogeneity $s_k/k$','interpreter','latex');
	axis square;
	
	splitfigure([2,3],[8,3],fflag);
	plot(s_a,mu0_a,'linewidth',1);
	ylabel('Initial condition $b(0)$','interpreter','latex');
	xlabel('Heterogeneity $s_k/k$','interpreter','latex');
	axis square;
	
	splitfigure([2,4],[8,5],fflag);
	cla
	plot(s_a,p_a,'linewidth',1);
	ylim([0,1]);
	ylabel('$P(b=\mathrm{high})$','interpreter','latex');
	xlabel('Heterogeneity $s_k/k$','interpreter','latex');
	axis square;
	
	splitfigure([2,4],[8,6],fflag);
	%yyaxis left
	cla();
	plot(s_a,[corr_a(:,1).^2,r2_S],'linewidth',1);
	legend('Pattern $b$','Density $S$','interpreter','latex','location','southeast');
	xlabel('Heterogeneity $s_K/K$','interpreter','latex');
	ylabel('Goodness of fit $R^2$','interpreter','latex');
	ylim([0,1.01]);
	axis square;
	if (0)
		yyaxis right;
		cla();
		plot(s_a,100*(1-pd));
		ylabel('Agreement / %');
	end
	
	splitfigure([2,4],[8,7],fflag);
	plot(s_a,rt_a)
	ylabel('Runtime / seconds');
	xlabel('Heterogeneity s_k');
	axis square;
	
	splitfigure([2,4],[8,8],fflag);
	plot(s_a,pd)
	ylabel('Area of disagreement');
	xlabel('Heterogeneity s_k');
	axis square;

	if (pflag)
		for adx=1:length(s_a);
			base = sprintf('img/may1977-sk_%g-L_%g-dx-%g',s_a(adx),Lx,dx);
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
	
		base = sprintf('img/may1977-L_g-dx-%g',Lx,dx);
		pdfprint(80+1,[base,'-filter-order-vs-sk.pdf'],ps);
		pdfprint(80+2,[base,'-lc-vs-sk.pdf'],ps);
		pdfprint(80+3,[base,'-b0-vs-sk.pdf'],ps);
		pdfprint(80+5,[base,'-p-dense-vs-sk.pdf'],ps);
		pdfprint(80+6,[base,'-goodness-of-fit-vs-sk.pdf'],ps);
		pdfprint(80+7,[base,'-runtime-vs-sk.pdf'],ps);
		pdfprint(80+8,[base,'-area-fraction-of-disagreement-vs-sk.pdf'],ps);
	end
	
	function val = fobj(mu0,s)
		prec = precompute(mu0,s);
		[mu0,s,prec.pp]
		val = prec.pp;
	end
	
	function [prec,rad] = precompute(mu0,s)
	
		param         = struct();
		param.opt.dt  = 1; % 2
		param.opt.dto = 10;
		param.T       = T;
		param.opt.path_str = 'mat/may/';
		param.nx      = nx*[1,1];
		param.L       = Lx*[1,1];
		
		param.initial_condition.mu = mu0;
		param.initial_condition.sd = 0;
		param.initial_condition.dist = {'normal'};
		param.pmu.k    = k;
		param.pss.k    = s; 
		param.psdist.k = 'gamma';
		param.opt.solver  = 'solve_split';
		param.opt.solver  = @ode23;
		
		rad    = May1977(param);
		f_str0 = rad.filename();
		f_str  = [f_str0(1:end-4),'-analyzed.mat']
		
		if (exist(f_str,'file'))
			load(f_str);
			load(f_str0,'rad','out');
		else
			[t,bb,out]=rad.run();
			prec.runtime = out.runtime;
			bb = single(bb);
			b = rad.extract2(bb(end,:));
			ae = reshape(rad.p.k,nx,nx);
			prec.x = rad.x;
			
			mu_a(adx,idx) = mean(b,'all');
			sd_a(adx,idx) = std(b,[],'all');
			me_a(adx,idx) = median(b,'all');
			bmin = min(b(:));
			bmax = max(b(:));
			g = bmin+graythresh((b-bmin)/(bmax-bmin))*(bmax-bmin);
			
			sd = std(bb,[],2);
			mu = mean(bb,2);
			
			Shat = abs(fft2(b-mean(b,'all'))).^2;
			Rhat    = ifft2(Shat);
			Rhat = Rhat/Rhat(1);
			
			[Sr,fr] = periodogram_radial(Shat,[Lx,Lx]);
			w = fr;
			%try
			[par,Srlp,fitstat]=fit_spectral_density(fr,Sr.normalized,fr,@lowpass2d_continuous_pdf,[1.5,15],'hellinger',0);
			prec.par = par;
			
			if (0)
				Cr=cumsum((cvec(fr).*cvec(Sr.normalized)));
				fc=interp1(Cr/Cr(end),fr,0.5);
			else
				% note, this is not a good estimate of fc, as the value near 0 is spoiled
				% setting it Sr(1) to max(Sr) fixes this somewhat
				Sr_ = Sr.normalized;
				Sr_(1) = max(Sr_);
				Ss = sort(Sr.normalized,'descend');
				fc = interp1(Ss,fr,0.5*Ss(1));
				Ss = sort(Sr.normalized,'descend');
				fc_lp = interp1(Srlp,fr,0.5*Srlp(1));
				%fc_a(adx,idx) = fc_lp;
			end
			% for plot
			Sr.normalized(1) = NaN;
			
			bb = single(bb);
			%fc = fc_lp;
			
			[St,theta] = periodogram_angular(Shat,[Lx,Lx]);
			
			[fx,fy,frr] = fourier_axis_2d([Lx,Lx],[nx,nx]);
			%Srlp = lowpass2d_pdf(fr,par(1),par(2),1);
			Slp2d = interp1(fr,Srlp,frr,'linear',0);
			blp = ifft2(sqrt(Slp2d).*fft2(reshape(ae,nx,nx)));
			me_b   = median(b,'all');
			me_blp = median(blp,'all');
			pp     = mean(flat(b)>g);
			q_     = quantile(blp(:),1-pp);
			b_thresh = (b>g);
			blp_thresh = (blp>q_);
			b_overlay   = b_thresh + 2*(blp_thresh);
			d      = rms(diff(bb),2)./(t(2)-t(1))./mean(bb(:));
			d(:,2) = diff(mean(bb,2))'./(t(2)-t(1))./mean(bb(:));
			d(:,3) = mid(rms(bb - bb(end,:),2));
			
			prec.d = d;
			prec.pp = pp;
			prec.t = t;
			prec.fc = fc;
			prec.fc_lp = fc_lp;
			prec.fx = fx;
			prec.fy = fy;
			prec.fr = fr;
			prec.p = pp;
			prec.b = b;
			prec.fx = fx;
			prec.b_overlay = b_overlay;
			prec.b_thresh = b_thresh;
			prec.blp_thresh = blp_thresh;
			prec.sd = sd;
			prec.mu = mu;
			prec.Sr = Sr;
			prec.Srlp = Srlp;
			prec.Shat = Shat;
			prec.Rhat = Rhat;
			prec.theta = theta;
			prec.St = St;
			prec.pd = mean(b_overlay == 1,'all') + mean(b_overlay == 2,'all');
			prec.g = g;	
		
			prec.corr(1) = corr(flat(blp),flat(b));
			% corr_a(adx,2) = corr(flat(Sr.normalized),flat(Srlp))
			%prec.corr(2) = wcorr(cvec(fr),sqrt(flat(Sr.normalized)),cvec(fr),sqrt(flat(Srlp))); 
			Sr.normalized(1) = 0;
			prec.corr(2) = corr(flat(blp_thresh),flat(b_thresh));

			prec.r2_S = fitstat.goodness.r2;
			%df = 1./param.L(1);
			% hellinger_distance(flat(Sr.normalized),flat(Srlp),df,fr);

			save(f_str,'prec');
		end % if not yet precomputed
	end % function precompute
end % function grazing_model_experiment

