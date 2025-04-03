% Karl Kästner, Berlin
% Tue 31 May 19:23:38 CEST 2022
% 2023-01-21 22:10
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
% generated 2D patterns with the rietkerk model for varying degrees of
% spatial heterogeneity of the bare soil infiltration
%
%% generate series of patterns with the Riektkerk model for vaying
%% degrees of spatial heterogeneity
%
function [tab,r2] = rk_2d_heterogeneity_experiment(meta,aniso)
	maxNumCompThreads(1)

	if (nargin()<1||isempty(meta))
		meta = pattern_formation_metadata();
	end
	if (nargin()<2)
		% 0 for spotted and 1 for striped patterns
		aniso = 0; 
	end
	fflag = meta.pflag;
	fflag = 1;
	ps    = 3.5;

	visible = meta.visible;
	scalefield = 'hat';

	%[param, vp, p_noise, sap,nkc] = rk_2d_heterogeneity_setup_intermediate(aniso);
	[param, vp, p_noise, sap,nkc] = rk_2d_heterogeneity_setup_complete(aniso);
	cva = vp.sa_times_a / param.pmu.a;

		sp_a = Spatial_Pattern(); 

		relstd = [];		
		Sc  = [];
		p_periodic  = [];
		cba = [];
		r2  = [];
		lc  = [];
		fi  = [];
		Si = [];
		xi = [];
		Ri = [];
		celerity = [];

		% variable parameters	
		name_C = {};
		val_C  = {};
		f_C = fieldnames(vp);
if (0)
		for idx=1:length(f_C)
			name_C{end+1} = f_C{idx};
			if (isnumeric(vp.(f_C{idx})))
				val_C{end+1} = num2cell(vp.(f_C{idx}));
			else
				val_C{end+1} = vp.(f_C{idx});
			end
		end
else
		name_C{end+1}  = 'opt.rng'; 
                val_C{end+1}   = num2cell(vp.seed);
		name_C{end+1}  = 'pss.a';
		val_C{end+1}   = num2cell(vp.sa_times_a);
		name_C{end+1}  = 'pmu.R';
		val_C{end+1}   = num2cell(vp.pmu.R);
		name_C{end+1}  = 'pmu.ey';
		val_C{end+1}   = vp.pmu.ey;
		name_C{end+1}  = 'psl.a';
		val_C{end+1}   = num2cell(vp.psl.a); 
end

		% initial condition (uniform / random)
		% heterogeneity spectrum (near white / near pink)
		% R = 
		% T/4
		% L/4
		 


		k  = 0;
		nn = prod(cellfun(@length,val_C));
		function run_(param)
			k = k+1;
                        param.opt.path_str = sprintf( ...
				'mat/server/vxh-%g-eyh-%g-R-%g-psl-%g-L-%g-%g-T-%g-seed-%g/' ...
				 , param.pmu.vx(3) ...
				 , param.pmu.ey(3) ...
				 , param.pmu.R ...	
				 , param.psl.a ...
				 , param.L...
				 , param.T ...
				 , param.opt.rng ...
				 );
			disp(param.opt.path_str)
			mkdir(param.opt.path_str);

			% run time estimate
			printf('Iteration k %d/%d CV(a) = %f\n',k,nn,param.pss.a./param.pmu.a);
			rk = Rietkerk(param);
		        [t, y, out] = rk.run();
			% continue model run
			if (length(param.T)>1)
				rk.T  = T(2);
				% quick fix
				rk.opt.path_str = rkmap.path_str;
				rk.opt.base_str = 'rietkerk-';
				[oname,oname_final] = rk.filename();
				if (exist(oname_final,'file'))
					% loading
					load(oname_final)
				else
					disp(['Continuing ', num2str(rk.hash)]);
					% create empty mat file as semaphore
					% for parallel computation
					empty = struct();
					save(oname_final,'-struct','empty');
					[t,y,out] = rk.continue_solve(t,y,T(2));
					rk.save(t,y,out);
				end
			end
			% write data as struct to be readable without class code
			filename = rk.filename;
			srad = struct(rk);
			%save(filename,'-append','srad');

			y = single(y);

			if (meta.analyze || meta.dflag)
				[sp, out] = rk_2d_heterogeneity_analyze(t,y,rk,aniso,p_noise,scalefield);
				Sc(k,1)    = sp.Sc;
				lc(k,1)    = sp.lambda_c;
				cba(k,:) = out.cba;


				p_periodic(k,1) = sp.stat.p_periodic;
				relstd(k,1) = out.relstd;
				if (0)
				if (aniso)
					fi      = out.fi.x;
					xi      = out.xi;
					Si(:,k) = out.Si.x.pdf.(scalefield);
					Ri(:,k) = out.Ri.x.(scalefield);
					r2(k,1) = sp.stat.fit.x.phase_drift.stat.goodness.r2;
				else
					fi      = out.fi.radial;
					xi      = out.xi;
					Si(:,k) = out.Si.radial.pdf.(scalefield);
					Ri(:,k) = out.Ri.radial.(scalefield);
					r2(k,1) = sp.stat.fit.x.bandpass.stat.goodness.r2;
				end
				end				
			end
			if (aniso)
if (0)
				rk.pmu.vy = rk.pmu.vx;
				rk.pmu.vx(:) = 0;
				rk.opt.dto = 1;
				rk.T   = 10;
				rk.opt.output_class = @single;
				rk.opt.solver = 'solve_split';
				rk.opt.inner_solver = 'step_advect_diffuse_spectral';

					rk.initial_condition = y(end,:);
				
				[t,y_] = rk.run();
end
				y_ = y;
				celerity(k,1:2,:) = rk.celerity(y_(end,:),true);
			end
			if (0)
			squeeze(celerity(:,1,:))
			squeeze(celerity(:,2,:))
			clf
			b1 = rk.extract2(y_(1,:));
			b2 = rk.extract2(y_(end,:));
			subplot(2,2,1)
			imagesc(b1);
			subplot(2,2,2)
			imagesc(b2);
			subplot(2,2,3)
			%plot([b1(1,:)',b2(1,:)']);
			plot([b1(:,1),b2(:,1)]);
			end
		%y_(1,:)',y_(end,:)']);
%			pause
			if (0) % meta.dflag || ismember(round(100*param.pss.a./param.pmu.a),round(100*sap)))
				rk_2d_heterogeneity_plot(rk,sp,out,k,aniso,nkc,meta,visible);
			%[Sc_,lc(k),cf,cS(:,k),cfi,cSi(:,k), ...
			%	pt(k),pt_,cba(k,:),r2(k),sp_a(k),RR(:,k)] = ...
			%	 rk_2d_heterogeneity_plot(t,y,rk,sp,k,aniso,nkc,meta,visible);
			%Sc(k) = Sc_;
			end % if dflag
		end % run

		% iterate, changing parameters accordingly
		iterate_cell_struct(@run_,param,name_C,val_C);

if (meta.analyze)
		regularity   = Sc ./ lc;
if (length(vp.seed)>1)
		% TODO, this does not work for iterator
			    regularity_p(idx,:) = lognfit(regularity(idx,:));
			    regularity_q(idx,:) = logninv([0.25,0.5,0.75],regularity_p(idx,1),regularity_p(idx,2));
			    z  = atanh(squeeze(cba));
			    mu = mean(z,3);
			    sd = std(z,[],3);
			    regularity_q_cba1 = tanh(norminv([0.25,0.5,0.75],mu(:,1),sd(:,1)));
			    regularity_q_cba2 = tanh(norminv([0.25,0.5,0.75],mu(:,2),sd(:,2)));
			    regularity_q_cba1(1,:) = 0;
			    regularity_q_cba2(1,:) = 0;
else
			req_q = regularity;
end % else of if length
end % if analyze


if (meta.analyze)
if (length(vp.seed)>1)
	for idx=1:size(regularity_p,1)
	   regularity_p(idx,:) = lognfit([regularity(idx,:)]);
	   regularity_q(idx,:) = logninv([0.16,0.5,0.84],regularity_p(idx,1),regularity_p(idx,2));
	   lc_p(idx,:) = lognfit([lc(idx,:)]);
	   lc_q(idx,:) = logninv([0.16,0.5,0.84],lc_p(idx,1),lc_p(idx,2));
	   Sc_p(idx,:) = lognfit([Sc(idx,:)]);
	   Sc_q(idx,:) = logninv([0.16,0.5,0.84],Sc_p(idx,1),Sc_p(idx,2));
	end % for idx
end % if length
	load(meta.filename.observed_patterns,'sp');

	% plot regularity vs heterogeneity
	splitfigure([2,2],[1e3,1],fflag,[],[],[],[],'Visible',visible);
	yyaxis left
	cla();
	if (length(vp.seed)>1)
		errorbar(cva,regularity_q(:,2),regularity_q(:,2)-regularity_q(:,1),regularity_q(:,3)-regularity_q(:,2),'o','markerfacecolor','k','linewidth',1,'markersize',3)
	else
		plot(cva,regularity,'ko','markerfacecolor','k','markersize',3);
	end
	hold on
	sa_ = mid(cva);
	p_periodic_ = median(p_periodic,2);
	%p_periodic_ = medfilt1([p_periodic_(1);p_periodic;p_periodic_(end)],3);
	%p_periodic_ = medfilt1([p_periodic_(1);p_periodic_(1);p_periodic_;p_periodic_(end);p_periodic_(end)],5);
	%p_periodic_ = p_periodic_(3:end-2);
	p_periodic_ = p_periodic;
	fdx = find(diff(p_periodic_>0.05)>0,1,'first');
%	set(gca,'xtick',0:dsap:1);
	xlim([0,max(cva)+sqrt(eps)]);
	ylim([0, 1.05*max(req_q(:))]);
	set(gca,'xtick',0:0.1:1);
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	axis square;
	if (~isempty(fdx))
		vline(sa_(fdx),'linestyle','-','color','red','linewidth',1);
	end
	if (aniso)
		ylabel('Regularity $S_{xc}/\lambda_c$','interpreter','latex');
	else
		ylabel('Regularity $S_{rc}/\lambda_c$','interpreter','latex');
	end
	yyaxis right
	cla
	plot(cva,relstd(:,1),'linewidth',1,'color',[0,0,0.7]);
	set(gca,'ycolor',[0,0,0.7])
	ylabel('Fraction of exogenous heterogeneity      ');
%Heterogeneity ratio $\displaystyle \frac{\mathrm{std}(a) \bar a_v}{\mathrm{std}(a_v) \bar a}$','interpreter','latex');
	drawnow();

	% plot max of density Sc vs heterogeneity
	splitfigure([2,2],[1e3,2],fflag,[],[],[],[],'Visible',visible);
	cla();
	if (length(vp.seed)>1)
		errorbar(cva,Sc_q(:,2),Sc_q(:,2)-Sc_q(:,1),Sc_q(:,3)-Sc_q(:,2),'o','markerfacecolor','k','linewidth',1,'markersize',3)
	else
		plot(cva,Sc,'ko','markerfacecolor','k','markersize',3);
	end
%	set(gca,'xtick',0:dsap:1);
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Density maximum $S_c$','interpreter','latex');	
	%ylim([0,625]);
	xlim([0,max(cva)+sqrt(eps)]);
	drawnow();
	axis square
	drawnow();

	% plot characteristic wavelength lambda_c
	splitfigure([2,2],[1e3,3],fflag,[],[],[],[],'Visible',visible);
	cla();
	if (length(vp.seed)>1)
		errorbar(cva,lc_q(:,2),lc_q(:,2)-lc_q(:,1),lc_q(:,3)-lc_q(:,2),'o','markerfacecolor','k','linewidth',1,'markersize',3)
		ylim([0,1.05*max(lc_q,[],'all')]);
	else
		plot(cva,lc,'ko','markerfacecolor','k','markersize',3);
		ylim([0,1.05*max(lc,[],'all')]);
	end
%	set(gca,'xtick',0:dsap:1);
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Wavelength $\lambda_c$ / m','interpreter','latex')
	axis square
	xlim([0,max(cva)+sqrt(eps)]);
	drawnow();

if (aniso)
	r2_ = r2; %arrayfun(@(x) x.phase_drift(end),r2);
	%plot(cva,r2.stoch(:,1),'-','linewidth',1);
else
	r2_ = r2; %arrayfun(@(x) x.bandpass(end),r2);
end
	% plot correlation between rad-pattern and bandpass pattern
	splitfigure([2,2],[1e3,4],fflag,[],[],[],[],'Visible',visible);
	cla();
	sai = linspace(cva(1),cva(end))';
	if (~aniso)
	plot(NaN,NaN,'k.-','linewidth',1);
	hold on
	plot(NaN,NaN,'k.--','linewidth',1);
	plot(NaN,NaN,'r.-','linewidth',1);
	else
	plot(NaN,NaN,'r-','linewidth',1);
	hold on
	end
if (length(vp.seed)>1)
	cba_ = mean(cba,3); 
	cbai = interp1(cvec(cva),cba_,sai,'pchip');
	plot(sai,cbai(:,1),'k-','linewidth',1);
	plot(sai,cbai(:,2),'k--','linewidth',1);
	errorbar(cva,regularity_q_cba1(:,2),regularity_q_cba1(:,2)-regularity_q_cba1(:,1),regularity_q_cba1(:,3)-regularity_q_cba1(:,2),'ko','markerfacecolor','k','linewidth',1,'markersize',3)
	errorbar(cva,regularity_q_cba2(:,2),regularity_q_cba2(:,2)-regularity_q_cba2(:,1),regularity_q_cba2(:,3)-regularity_q_cba2(:,2),'ko','markerfacecolor','k','linewidth',1,'markersize',3)
else
	if (~aniso)
		plot(cva,cba(:,1).^2,'k-','linewidth',1);
		plot(cva,cba(:,2).^2,'k--','linewidth',1);
		plot(cva,r2_,'r-','linewidth',1);
	else
		plot(cva,r2_,'r-','linewidth',1);
	end
end
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Goodness of fit','interpreter','latex');
	if (aniso)
		legend('$R^2(S_b,S_{PNI})$', ... %'$R^2(b,a)$','$R^2(S_b,S_{PNI})$', ...
				'location','southeast','interpreter','latex');
		%legend('$R^2(b,PNI(a))$','$R^2(b,a)$','$R^2(S_b,S_{PNI})$', ...
		%		'location','southeast','interpreter','latex');
	else
		legend('$R^2(b,BP(a))$','$R^2(b,a)$','$R^2(S_b,S_{BP})$', ...
				'location','southeast','interpreter','latex');

	end
	xlim([0,max(cva)+sqrt(eps)]);
	ylim([0,1]);
	set(gca,'xtick',0:0.1:1);
%	set(gca,'xtick',0:dsap:1);
	ylim([min(0,min(cba(:))),1]);
	axis square
	set(gca,'ytick',0:0.2:1);
if (0)
	yyaxis right
	plot(cva,r2_,'-','linewidth',1);
	ylim([0,1]);
	ylabel('R$^2$ of density fit');
end
	drawnow()

	% p-value of periodicity test
	splitfigure([2,2],[2e3,1],fflag,[],[],[],[],'Visible',visible);
	cla
	plot(cva,p_periodic);
	hold on
	plot(cva,p_periodic_);
	xlim([0,max(cva)+sqrt(eps)]);
	ylabel('P-test');
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	axis square
	drawnow()

	% plot spectral density
	if (~isempty(sap))
	splitfigure([2,2],[2e3,2],fflag,[],[],[],[],'Visible',visible);
	cla();
	fdx = ismember(round(100*cva),round(100*sap));
	plot(fi,mean(Si(:,fdx,:),3),'-','linewidth',1);
	xlabel('Wave number $k_x/k_c$','interpreter','latex');
	if (aniso)
	ylabel('Density $S_x / \lambda_c$','interpreter','latex')
	else
	ylabel('Density $S_r / \lambda_c$','interpreter','latex')
	end
	%axis square
	set(gca,'colororder',meta.colormap);

	% plot empirical
	hold on
if (0)
	if (1 == aniso)
		fdx_ = sp(1).f.x>0;
		s = sum(sp(1).S.rot.x.(field)(fdx_))*diff(sp(1).f.x(1:2));
		plot(sp(1).f.x(fdx_)/sp(1).stat.fc.x.(field),1./s*sp(1).S.rot.x.(field)(fdx_)*sp(1).stat.fc.x.(field),':','linewidth',1.5)
		regularity_natural = sp(1).stat.Sc.x.(field).*sp(1).stat.fc.x.(field);
	else
		s = sum(sp(2).S.radial.(field))*diff(sp(2).f.r(1:2));
		plot(sp(2).f.r/sp(2).stat.fc.radial.(field),1./s*sp(2).S.radial.(field)*sp(2).stat.fc.radial.(field),':','linewidth',1.5)
		regularity_natural = sp(2).stat.Sc.radial.(field).*sp(1).stat.fc.radial.(field);
	end
end 
	xlim([0,2.5])
	ylim([0,7]);
	axis square
	leg = arrayfun(@(x) num2str(x,'RK CV(a) = %0.2f'),cva(fdx),'uniformoutput',false);
	leg{end+1} = 'natural';
	lh = legend(leg{:});
	drawnow()

if (0)
	splitfigure([2,2],[1,4],fflag,[],[],[],[],'Visible',visible);
	hline(regularity_natural,'linestyle',':','color','k','linewidth',1.5);
end

	% plot autocorrelation
	splitfigure([2,2],[2e3,3],fflag,[],[],[],[],'Visible',visible);
	cla()
	fdx = find(ismember(round(100*cva),round(100*sap)));
	fdx
	for idx=rvec(fdx)
		%if (0 == aniso)
		%	x = sp_a(idx,1).r;
		%	fc=sp_a(idx,1).stat.fc.radial.(field);
		%else
		%	x = sp_a(idx,1).x;
		%	x = x-x(1);
		%	fc = sp_a(idx,1).stat.fc.x.(field);
		%end
		if (aniso)
			plot(xi,Ri(:,idx),'linewidth',1);
		else
			x0 = 0.1497;
			fvalid = xi>=x0;
			plot(xi(fvalid),pi*sqrt(cvec(xi(fvalid))).*Ri(fvalid,idx),'linewidth',1);
		end
		%mean(RR(:,idx,:),3),
		hold on
	end
	xlim([0,2.5]);
	axis square
	if (0)
	% TODO real pattern
	if (0 == aniso)
		plot(sp(1).r.*sp(1).stat.fc.radial.bar,sp(1).R.radial.bar,'k:','linewidth',1);
	else
		x = sp(2).x;
		x = x-x(1);
		plot(x.*sp(2).stat.fc.x.bar,sp(2).R.rot.x.bar,'k:','linewidth',1);
		%plot(x.*sp(2).stat.fc.x.bar,sp(2).R.rot.x.bar,'k--');
	end
	end
	if (aniso)
		xlabel('Lag Distance $x/\lambda_c$','interpreter','latex');
		ylabel('Autocorrelation $R_x$','interpreter','latex');
	else
		xlabel('Lag Distance $r/\lambda_c$','interpreter','latex');
		ylabel('Autocorrelation $\pi \sqrt{r/\lambda_c} R_r$','interpreter','latex');
	end
	%lh=legend(num2str(cvec(sap)));
	%title(lh,'Heterogeneity $cva$','interpreter','latex');
	set(gca,'colororder',meta.colormap);
	drawnow()
%catch e
%	disp(e);
%end
end % if sap

	splitfigure([2,2],[1e3,4],fflag,[],[],[],[],'Visible',visible); 
	cla();
if (0)
	r2.stoch = median(r2.stoch,2);
	r2.logn  = median(r2.logn,2);
	r2.gamma = median(r2.gamma,2);
%	r2.white = median(r2.white,2);
	r2.periodic = median(r2.periodic,2);
	printf('median(r2): %f\n',median(r2.stoch));
	printf('median(r2_logn): %f\n',median(r2.logn));
	printf('median(r2_gamma): %f\n',median(r2.gamma));
end
	%plot(cva,[r2.stoch,r2.logn,r2.gamma,r2.white,r2.periodic],'.-','linewidth',1);
%	r2_ = [arrayfun(@(x) x.bandpass(2),r2), arrayfun(@(x) x.phase_drift(2),r2), arrayfun(@(x) x.logn(2),r2),arrayfun(@(x) x.gamma(2),r2),arrayfun(@(x) x.white(2),r2),arrayfun(@(x) x.periodic(2),r2)];
%	median(r2_)
if (0)
%	plot(cva,r2_,'linewidth',1);
end
%	ylim([0,1]);
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	if (aniso)
		plot(cva,r2,'.','color',[0,0,0.6]);
		ylabel('Goodness of density fit $R^2_{S_x+}$','interpreter','latex');
		%ylabel('R$^2$ of $S_x$ density fit','interpreter','latex');
		%legend({'PNI','Log-normal','Gamma'},'location','southeast');
		tab.r2  = r2;
	else
		ylabel('R$^2$ of $S_r$ density fit','interpreter','latex');
		legend({'BP','Log-normal','Gamma'},'location','southeast');
	end
	ylim([0,1]);
	set(gca,'colororder',meta.colormap);
	axis square
	drawnow()
	str = sprintf('img/rk-2d-sl-%g-vxh-%g-eyh-%g-R-%g-L-%d-T-%1.0e',param.psl.a,param.pmu.vx(3),vp.pmu.ey{1}(3),vp.pmu.R(1),param.L(1),param.T(end));
        xlabel('Exogenous Heterogeneity $CV(a)$','interpreter','latex');          
	ps = 4;
%	pdfprint(1e4+4,[str,'-r2-density-fit'],ps);

	tab = table();
	tab.cva = cvec(cva);	
	tab.Sc = cvec(Sc);
	tab.wavelength_c = cvec(lc);
	tab.regularity = cvec(regularity);
	tab.p_periodic = cvec(p_periodic);
	tab.correlation_ba = cba;
%	tab.r2  = r2;

	fflag = true();

	splitfigure([2,2],[3e3,1],fflag);
	%celerity
	celerity = squeeze(celerity(:,2,1));
	celerity = celerity*sign(celerity(1));
	plot(cva,celerity,'.','color',[0,0,0.6]);
	xlabel('Exogenous Heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Migration celerity $c$ / (m/d)','interpreter','latex');
	xlim([0,cva(end)]);
	ylim([-sqrt(eps),1.05*max(celerity)]);
	axis square
%	pdfprint(3e3*10+1,[str,'-migration-celerity.pdf'],ps);

	if (0) %meta.pflag)
		aspect = [];
	
		pdfprint(11,[str,'-regularity-vs-sa.pdf'],ps,aspect);
		pdfprint(12,[str,'-Sc-vs-sa.pdf'],ps,aspect);
		pdfprint(13,[str,'-wavelength-vs-sa.pdf'],ps,aspect);
		pdfprint(14,[str,'-correlation-a-b-vs-sa.pdf'],ps);

		pdfprint(21,[str,'-periodicity-test.pdf'],ps);
		if (aniso)
		%pdfprint(22,[str,'-density-Sx.pdf'],ps);
		%pdfprint(23,[str,'-acf-Rx.pdf'],ps);
		else
		%pdfprint(22,[str,'-density-Sr.pdf'],ps);
		%pdfprint(23,[str,'-acf-Rr.pdf'],ps);
		end
	%	pdfprint(2e5+3,[str,'-autocorrelation.pdf'],ps);
		pdfprint(2e5+4,[str,'-r2-density-fit'],ps);
		%pdfprint(100001,sprintf('img/rietkerk-2d-regularity-vs-sa-vh-%g-R-%g-L-%d-T-%1.0e.pdf',rk.pmu.vh(2),rk.pmu.R,L(1),T),ps,aspect);
		%pdfprint(100002,sprintf('img/rietkerk-2d-Sc-vs-sa-vh-%g-R-%g-L-%d-T-%1.0e.pdf',rk.pmu.vh(2),rk.pmu.R,L(1),T),ps,aspect);
		%pdfprint(100003,sprintf('img/rietkerk-2d-wavelength-vs-sa-vh-%g-R-%g-L-%d-T-%1.0e.pdf',rk.pmu.vh(2),rk.pmu.R,L(1),T),ps,aspect);
		%pdfprint(100004,sprintf('img/pattern-2d-correlation-a-b-vh-%g',(param.pmu.vh(2))),ps);
	end
end % if analyze
end % rk_2d_heterogeneity_experiment

