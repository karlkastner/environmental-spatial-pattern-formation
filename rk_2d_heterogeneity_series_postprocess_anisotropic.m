% 2025-03-21 20:46:59.824292501 +0100

mat_filename = 'mat/pattern-formation-statistic-isotropic.mat';

if (exist(mat_filename,'file'))
	load(mat_filename);
else
	folder = 'mat/server/vxh-10-eyh-20-R-0.9-L-1024-1024-T-500000-seed-0/';
	f_a=dir([folder,'/*final.mat']);
	clear sp_a
	tab = table();
	for idx=1:length(f_a)
		disp(idx);
		f =[f_a(idx).folder,'/',f_a(idx).name];
		f_analyzed = [f(1:end-4),'-analyzed.mat'];
		clear t y rad cp c
		load(f);
		if (exist(f_analyzed,'file'))
			load(f_analyzed);
		else
		[b,h] = rad.extract2(y(:,2));
		a = rad.p.a;
		sp = Spatial_Pattern();
		if (1 ~= numel(a))
			sp.source = reshape(a,rad.nx);
		end
		% pattern
		sp.b = b;
		% spatial extent
		sp.L = rad.L;
		% postprocess pattern
		sp.opt.suppress_low_frequency_components = 0;
		sp.opt.angle_deg = 0;
		sp.analyze_grid();
		sp.fit_parametric_densities();

		a = rad.p.a;
		% decompose variance of infiltration infiltration
		if (~isscalar(a))
		%cva(idx,1) = rad.pss.a/rad.pmu.a;
			bb(:,idx) = flat(b);
			hh(:,idx) = flat(h);
			aa(:,idx) = flat(rad.p.a);
			ie_ = rad.infiltration_enhancement(b);
			ie(:,idx) = flat(ie_);
			I(:,idx) = aa(:,idx).*ie(:,idx).*hh(:,idx);
			% note that the mean is subtracted, so no column for the mean is required
			A = [   (aa(:,idx) - mean(aa(:,idx)))/std(aa(:,idx)), ...
				(ie(:,idx)-mean(ie(:,idx)))/std(ie(:,idx)),...
				(hh(:,idx)-hh(:,idx))/std(hh(:,idx))];
			c = A\((I(:,idx)-mean(I(:,idx)))/std(I(:,idx)));
		 else
			c= [NaN,NaN,NaN];
		 end
		cc(idx,:)=c;

		save(f_analyzed,'sp','c');
		end
		sp_a(idx,1)   =sp;
		%$tab.pssa(idx) = rad.pss.a;
		tab.cva(idx) = rad.pss.a/rad.pmu.a;
		tab.ca(idx)   = c(2);
	end % for idx
	sp_a = cvec(sp_a);

	[pssa_,sdx]=sort(tab.cva);
	f_C = {'con','bandpass','logn','gamma','normalmirrored','phase_drift'};

	Sc = [];
	fc = [];
	r2 = [];
	p_periodic = [];

	% extract parameters
	for idx=1:length(f_C)
		fc=[fc,arrayfun(@(x) x.stat.fc.xp.(f_C{idx}),sp_a)];
		Sc=[Sc,arrayfun(@(x) x.stat.Sc.xp.(f_C{idx}),sp_a)];
		try
		r2(:,idx) = [arrayfun(@(x) x.stat.fit.xp.(f_C{idx}).stat.goodness.r2,sp_a)];
		catch
		end
	end

	reg = Sc.*fc;

	% TODO migration rate
	tab.Scxp  = Sc(:,end);
	tab.fcxp  = fc(:,end);
	tab.regx  = reg(:,end);
	tab.p_periodic = arrayfun(@(x) x.stat.p_periodic,sp_a);
	tab.r2Sxp = r2(:,end);

else
	save(mat_filename,'tab');
end % if not exist

figure(1);
clf
subplot(2,3,1)
plot(tab.cva(sdx),tab.Sc(sdx,:),'.-')
legend(f_C)
subplot(2,3,2)
plot(tab.cva(sdx),1./tab.fc(sdx),'.')
ylim([0,250]);
subplot(2,3,3)
plot(tab.cva(sdx),tab.Sc(sdx,:).*tab.fc(sdx),'.-')
legend(f_C)

subplot(2,3,4)
plot(tab.cva(sdx),tab.r2Sxp(sdx,:),'.-');
legend(f_C)
% plot(pssa,Sc.*fc,'.');

subplot(2,3,5)
plot(tab.cva(sdx),tab.p_periodic(sdx),'.-');

figure(2);
clf;
plot(NaN)
yyaxis right
plot(tab.cva,tab.ca,'o','markersize',2,'color','none','markerfacecolor',[0,0,0.8])
set(gca,'ycolor',[0,0,0.8])
ylim([0,1]);
set(gca,'ytick',0:0.2:1);
xlabel('Exogeneous heterogeneity $CV(a)$','interpreter','latex');
ylabel('Fraction of exogenous heterogeneity $s_a$','interpreter','latex');
axis square;
xlim([0,0.5])

tabs = sortrows(tab);
tabs = round(tabs,3,'significant');
writetable(tabs,[folder,'/rietkerk-striped.csv']);
writetable(tabs,'output/rietkerk-striped.csv');

if (pflag)
	ps = 4;
	pdfprint(1,'img/anisotropic-fraction-of-exogenous-heterogeneity.pdf',ps);
end

