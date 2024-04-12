% Fri  4 Mar 14:19:02 CET 2022
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
% demonstrate the phase shift a banded pattern experiences when it is
% perturbed, here the perturbation is a locally deviating bare soil
% infiltration coefficient
%
% TODO reduce width of perturbation
% TODO initial wavelength of 250 m to identify
%
function [x,bb,w,tp,rk] = rk_1d_phase_shift_experiment_single_bump(meta)
	if (nargin()<1)
		meta = pattern_formation_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	lw = 1;
	
	% domain size
	L  = 2*meta.example.L;
	% spatial discretisation
	dx = 2*meta.example.dx;
	% perturbation of the bare soil infiltration coefficient
	sd_a = 0;
	wwin   = 3; % if larger, left transition influences right transition
	L0 = 200;
	da = [0.2] %0.199 0.2 0.201];
	ap = 1-da;
	sig = -1;
	vh   = sig*meta.example.vh;
	eh   = meta.mopt.eh;
	lambda0 = 250;
	db = 0;
	fc0 = 1./lambda0;
	fc0 = [];
	% initialize random number generator

	rkmap = Rietkerk_Map('path_str','mat/');
	rkmap.init();
	rkmap.opt.hashvectors = true;

	for ldx = 1:length(L0)
	for adx = 1:length(da)

	param = struct();
	param.L = L;
	param.n = round(L/dx);
	param.T = meta.mopt.To;
	param.dt = meta.mopt.dt;
	param.dto = meta.mopt.dto;
	param.pss.a  = sd_a;
	param.pmu.eh = eh;
	param.pmu.vh = vh;
	param.pst.db = db;
	param.initial_condition = 'random';
	if (isempty(fc0))
		[t,y,rk]   = rkmap.run(param);
		y0 = rk.initial_condition_from_central_frequency(y(end,:));
	else
		rk = Rietkerk(param); 
		y0 = rk.initial_condition_periodic(fc);
	end

	param.initial_condition = y0;
	% run model
	[t,y,rk] = rkmap.run(param);
	y = double(y);
	
	[bb,ww,hh] = rk.extract2(y);

	[b,w,h] = rk.extract1(y(end,:)');

	a = rk.pmu.a;

	n = length(b);
	%y0 = [];
	%y0(:,1) = 0.1*rand(n,1)*mean(b);
	%y0(:,2) = mean(w)*ones(n,1);
	%y0(:,3) = mean(h)*ones(n,1);
	y0 = y0(:);
	T = 20*meta.mopt.To; % 20
	dt = meta.mopt.dt;
	x0 = [0.1,0.6]*L;
	param.pmu.a = a.*(1 + ...
		(ap(adx)-1)*( ...
				  exp(-0.5*(cvec((rk.x-x0(1)))/L0).^2) ...
				+ exp(-0.5*(cvec((rk.x-x0(2)))/L0).^2) ...
	    )); 
	param.pss.a = 0;
	param.T =  T;
	param.dt = dt;
	param.initial_condition = y0;
	[t,y,rk]   = rkmap.run(param);
	[b,w,h] = rk.extract1(y(end,:)');
	[bb,ww,hh] = rk.extract1(y);

	Li = logspace(log10(50),log10(350),2*200);
	b_cwt = cwt_man2(b,Li,dx,wwin);
	[b_c,mdx] = max(abs(b_cwt),[],2);
	Lc = Li(mdx);


	phi = zeros(length(mdx),1);
	for idx=1:length(mdx)
		phi(idx,1) = angle(b_cwt(idx,mdx(idx)));
	end
	phi   = cvec(phi);
	phi   = unwrap(phi);
	fdx   = rk.x > 0.2*L & rk.x < 0.45*L;
	dphi  = median(diff(phi(fdx)));
	phi  = phi - dphi*(0:length(phi)-1)';

	% mean
	lambda_c = mean(Lc(fdx));
	printf('lambda_c %g\n',lambda_c);
	b_bar    = mean(b(fdx));
	h_bar    = mean(h(fdx));
	a_bar = mean(rk.p.a(fdx));
	b_c_bar = mean(b_c(fdx));

	xl = x0(2) + sort(sig*lambda_c*[-10,5]);
	%xl = x0(2) + sort(sig*lambda_c*[-12.5,7.5]);
	fdx_ =find(rk.x>=xl(1) & rk.x<=xl(2));

	% filter whiggles
	nf = round(2*lambda_c/dx);
	b_c_bar = trifilt1(b_c_bar,nf);
	phi  = trifilt1(phi,nf);
	% must come after filtering
	phi  = wrapToPi(phi);
	Lc = trifilt1(Lc,nf);
	b_c = trifilt1(b_c,nf);
	
	% plot bare soil infiltration
	splitfigure([2,3],[(adx-1)*10+1,1],fflag,'',100);
	plot((rk.x-x0(2))/lambda_c,rk.p.a/a_bar,'linewidth',lw);
	xlim((xl-x0(2))/lambda_c);
	xlabel('$x/\bar\lambda_c$','interpreter','latex');
	ylabel('$\displaystyle\frac{a}{\bar a}$','interpreter','latex','rot',0);
	title('Bare soil infiltration');
	ylim([0.75,1.05])

	% plot biomass
	splitfigure([2,3],[(adx-1)*10+1,2],fflag,'',100);
	plot((rk.x-x0(2))/lambda_c,b/b_bar,'linewidth',0.5)
	xlim((xl-x0(2))/lambda_c);
	xlabel('$x/\bar\lambda_c$','interpreter','latex')
	ylabel('$\displaystyle\frac{b}{\bar b}$','interpreter','latex','rot',0)
	title('Biomass');
	vline((x0-x0(2))/lambda_c,'color','r');

	% plot wavelength
	splitfigure([2,3],[(adx-1)*10+1,3],fflag,'',100);
	yyaxis left
	cla
	plot((rk.x-x0(2))/lambda_c,Lc/lambda_c,'linewidth',lw)
	xlim((xl-x0(2))/lambda_c);
	xlabel('Distance $x/\bar\lambda_c$','interpreter','latex')
	ylim([0,1.05]);
	hold on
	plot((rk.x-x0(2))/lambda_c,b_c/b_c_bar,'-r','linewidth',lw);
	xlim((xl-x0(2))/lambda_c);
	ylim([0,1.05]);
	legend('Wavelength $\lambda_c/\bar \lambda_c$','Amplitude $b_c/\bar b_c$','interpreter','latex','location','southeast');
	ylabel('$\displaystyle\frac{\lambda_c}{\bar\lambda_c}$','interpreter','latex','rot',0)
	title('Wavelength of central frequency component');
	yyaxis right;
	cla();
	phi = wrapToPi(phi-phi(fdx_(1)));
	plot((rk.x-x0(2))/lambda_c,phi,'color',[0,0,0.8],'linewidth',0.5);
	ylim(1.1*pi*[-1.2,1]);
	xlim((xl-x0(2))/lambda_c);
	xlabel('Distance $x/\bar\lambda_c$','interpreter','latex');
	set(gca,'ycolor',[0,0,0.8]);
	set(gca,'ytick',[-pi/2,0,pi/2],'yticklabel',{'$\displaystyle-\frac{\pi}{2}$','0','$\displaystyle+\frac{\pi}{2}$'},'yticklabelrot',0,'TickLabelInterpreter','latex'); 
	hline(0,'color',[0.5,0.5,0.5],'linestyle','-');
	legend('Wavelength $\lambda_c/\bar \lambda_c$','Amplitude $b_c/\bar b_c$','Phase shift $\varphi_c-\bar \varphi_c$','interpreter','latex','location','southeast');
	yyaxis left;
	ylim([0,2]);
	yyaxis right;
	ylim([-1,1]*pi);

	% plot amplitude
	splitfigure([2,3],[(adx-1)*10+1,4],fflag,'',100);
	plot((rk.x-x0(2))/lambda_c,b_c/b_c_bar,'linewidth',lw);
	xlim((xl-x0(2))/lambda_c);
	xlabel('Distance $x/\bar\lambda_c$','interpreter','latex');
	ylabel('Biomass $\displaystyle\frac{b_c}{\bar b_c}$','interpreter','latex','rot',90);
	title([da(adx),lambda_c]);

	% plot phase
	splitfigure([2,3],[(adx-1)*10+1,5],fflag,'',100);
	cla();
	plot((rk.x-x0(2))/lambda_c,phi,'linewidth',0.5);
	xlim((xl-x0(2))/lambda_c);
	xlabel('Distance $x/\bar\lambda_c$','interpreter','latex');
	ylabel('Phase shift $\varphi_c-\bar \varphi_c$','interpreter','latex','rot',90);
	hline(0,'color',[0.5,0.5,0.5]);

	% plot biomass
	splitfigure([2,3],[(adx-1)*10+1,6],fflag,'',100);
	yyaxis right
	cla
	yyaxis left
	cla
	plot((rk.x-x0(2))/lambda_c,b/b_bar,'linewidth',0.5)
	xlim((xl-x0(2))/lambda_c);
	xlabel('Distance $x/\bar\lambda_c$','interpreter','latex')
	ylim([0,7]);
	ylabel('Biomass $b/\bar b$','interpreter','latex');
	yyaxis right
	plot((rk.x-x0(2))/lambda_c,rk.p.a/a_bar,'linewidth',1,'linewidth',lw);
	xlim((xl-x0(2))/lambda_c);
	ylim([0.75,1.15]);
	legend('Biomass $b/\bar b$','Bare soil infiltration $a/\bar a$','interpreter','latex');
	ylabel('Bare soil infiltration $a/\bar a$','interpreter','latex');

	%figure(3)
	splitfigure([2,2],[(adx-1)*10+2,1],fflag,'',100);
	cla();
	%subplot(2,2,1)
	imagesc(rk.x,Li,abs(b_cwt)')

	%subplot(2,2,2);
	splitfigure([2,2],[(adx-1)*10+2,2],fflag,'',100);
	cla
	imagesc(rk.x,Li,abs(b_cwt)')

	%subplot(2,2,3)
	splitfigure([2,2],[(adx-1)*10+2,3],fflag,'',100);
	cla
	imagesc((rk.x-x0(2))/lambda_c,t,bb);
	xlabel('$x/\bar\lambda_c$','interpreter','latex')



	splitfigure([2,2],[(adx-1)*10+2,4],fflag,'',100);
	cla
	%figure(100)
	%clf
	x=rk.x;
	xi = xl(1):(0.1*(x(2)-x(1))):xl(2);
	bi = interp1(x(fdx_),trifilt1(b(fdx_),3),xi,'spline');
	hi = interp1(x(fdx_),trifilt1(h(fdx_),3),xi,'spline');
	%wi = interp1(x(fdx_),trifilt1(w(fdx_),3),xi,'spline');
	plot(hi/h_bar,bi/b_bar);
	hold on
	gdx = (xi-xi(1))<lambda_c;
	plot(hi(gdx)/h_bar,bi(gdx)/b_bar,'linewidth',2);
	gdx = (xi(end)-xi)<lambda_c;
	h = plot(hi(gdx)/h_bar,bi(gdx)/b_bar,'b--','linewidth',2);
	h.HandleVisibility = 'off';
	plot(NaN,NaN,'b-','linewidth',2)
	xlabel('$h/\bar h$','interpreter','latex');
	ylabel('$b/\bar b$','interpreter','latex');
	legend('at perturbation','uphill','downhill','location','northwest');

	if (pflag)
		ps = 4.5;
		base = sprintf('img/phase-shift-experiment-local-perturbation-dl_%g-da_%0.2f',L0(ldx),da(adx));
		pdfprint(101,[base,'-a.pdf'],ps);
		pdfprint(102,[base,'-b.pdf'],ps);
		ps_ = 3.5;
		pdfprint(103,[base,'-lambda-and-amplitude.pdf'],ps_);
		pdfprint(104,[base,'-amplitude.pdf'],ps);
		pdfprint(105,[base,'-phase.pdf'],ps);
		pdfprint(106,[base,'-a-and-b.pdf'],ps_);
	end
	end % for adx
	end % for ldx
end  % plot_rietkerk_dzdt_1d_initial


