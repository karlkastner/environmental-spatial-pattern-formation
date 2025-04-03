% 2025-03-21 20:46:59.824292501 +0100

if (~exist('tab','var'))
folder = '/home/pia/large/cottbus/fourier-vegegation/mat/server/vxh-10-eyh-20-R-0.9-L-1024-1024-T-500000-seed-0/';
f_a=dir([folder,'/*final.mat']);
clear sp_a 
tab = table();
for idx=1:length(f_a)
	idx
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
	sp.b = b;
	sp.L = rad.L;
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
		cc(idx,:)=c;
	 catch
		c(idx,1) = [NaN,NaN,NaN];
	 end



	save(f_analyzed,'sp','c');
	end
	sp_a(idx,1)   =sp;
	tab.pssa(idx) = rad.pss.a;
	tab.ca        = c(2);
end
end
sp_a = cvec(sp_a);
%  Sxc fxc reg cx R2Sx R2b pper cmigration

 clf;
 pssa=tab.pssa;
n = length(pssa);
 [pssa_,sdx]=sort(pssa);
 id(sdx,1)=(1:51)';
f_C = {'con','bandpass','logn','gamma','normalmirrored','phase_drift'};

Sc = [];
fc = [];
r2 = [];
p_periodic = [];

%for idx=1:51; 
%	sp_a(idx).stat.fc.xp.phase_drift=sp_a(idx).stat.Sc.fc.xp.phase_drift;
%end


for idx=1:length(f_C)
	fc=[fc,arrayfun(@(x) x.stat.fc.xp.(f_C{idx}),sp_a)];
	Sc=[Sc,arrayfun(@(x) x.stat.Sc.xp.(f_C{idx}),sp_a)];
	try
	r2(:,idx) = [arrayfun(@(x) x.stat.fit.xp.(f_C{idx}).stat.goodness.r2,sp_a)];
	catch
	end
end

%h = Sc(:,end);
%Sc(:,end) = fc(:,end);
f%c(:,end)= h;

reg = Sc.*fc;

tab.Scxp  = Sc(:,end); 
tab.fcxp  = fc(:,end); 
tab.regx  = reg(:,end);
tab.p_periodic = arrayfun(@(x) x.stat.p_periodic,sp_a);
tab.r2Sxp = r2(:,end);
% TODO migration rate


	%,arrayfun(@(x) x.stat.fc.xp.bandpass,sp_a)];
% Sc = cvec(arrayfun(@(x) x.stat.Sc.xp.con,sp_a));
% Sc(:,2) = arrayfun(@(x) x.stat.Sc.xp.bandpass,sp_a);

figure(1);
clf
subplot(2,3,1)
plot(pssa(sdx),Sc(sdx,:),'.-')
legend(f_C)
subplot(2,3,2)
plot(pssa,1./fc,'.')
ylim([0,250]);
subplot(2,3,3)
plot(pssa(sdx),Sc(sdx,:).*fc(sdx),'.-')
legend(f_C)

subplot(2,3,4)
plot(pssa(sdx),r2(sdx,:),'.-'), mean(r2), median(r2);
legend(f_C)
% plot(pssa,Sc.*fc,'.');

subplot(2,3,5)
plot(pssa(sdx),tab.p_periodic(sdx),'.-');

% TODO migration rate 
%tab.p_periodic = p_periodic;
tabs = sortrows(tab);
tabs.cva = tabs.pssa/rad.pmu.a;
tabs=[tabs(:,end),tabs(:,2:end-1)];
tabs = round(tabs,3,'significant');
writetable(tabs,[folder,'/rietkerk-striped.csv']);
writetable(tabs,'output/rietkerk-striped.csv');

figure(2);
clf;
plot(NaN)
yyaxis right
plot(cva,cc(:,1),'o','markersize',2,'color','none','markerfacecolor',[0,0,0.8])
set(gca,'ycolor',[0,0,0.8])
ylim([0,1]);
set(gca,'ytick',0:0.2:1);
xlabel('Exogeneous heterogeneity $CV(a)$','interpreter','latex');
ylabel('Fraction of exogenous heterogeneity $s_a$','interpreter','latex');
axis square;
xlim([0,0.5])
if (pflag)
	ps = 4;
	pdfprint(1,'img/anisotropic-fraction-of-exogenous-heterogeneity.pdf',ps);
end

