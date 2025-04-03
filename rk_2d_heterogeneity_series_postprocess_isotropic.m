% Wed  5 Mar 15:23:08 CET 2025

% TODO merge with analyze

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
if (pflag)
	markersize =2;
else
	markersize = 5;
end

theta = [256];
label = {'spotted','labyrinthine','gapped'};
R_ = [0.8,1,1.15]
clear sp_a; 

for idx=1:3
	fh = figure(idx);
	clf
	fh.Visible = false;
end

if (1) %~exist('tab','var'))

for tdx=1:length(theta)

for jdx=1:length(R_)

folder = sprintf('mat/server/vxh-0-eyh-100-R-%g-psl-%g-L-1024-1024-T-50000-seed-0/',R_(jdx),theta(tdx));

file_C = dir([folder,'/*-final.mat']);
tab{jdx,tdx} = table();

bb = NaN(1024,1024,51);
for idx=1:length(file_C)
disp([jdx,idx]);

f = [file_C(idx).folder,'/',file_C(idx).name];
disp(f)
analysis_str = [f(1:end-10),'-analyzed.mat'];
%sprintf('mat/server/vxh-0-eyh-100-R-%g-psl-%g-L-1024-1024-T-50000-seed-0/*final.mat',R_(jdx),theta(tdx));

%if (exist(f,'file'))
clear rad sp y
%try
load(f);
tab{jdx,tdx}.cva(idx) = rad.pss.a/rad.pmu.a;
tab{jdx,tdx}.name{idx} = f;
if (~isempty(y))
[b,w,h] = rad.extract2(y(:,2));
bb(:,:,idx) = b;
tab{jdx,tdx}.b(idx) = mean(b,'all');
%$f
%[idx,rms(rad.pss.a),rms(y),rms(b,'all'),rms(bb(:,:,idx),'all')]
else
%	idx
	tab{jdx,tdx}.b(idx) = NaN;
end % else of isemtpy(y)

if (mean(b,'all')>0)

if (isfield(rad.p,'a') && ~isscalar(rad.p.a))
	a = rad.p.a;
	a = reshape(a,rad.nx);
	Sa = rad.psS.a;
	Sa = max(0,real(Sa));
	kmax = 2;

ie = rad.infiltration_enhancement(flat(b));
infiltration = flat(a).*rad.infiltration_enhancement(flat(b)).*flat(h);
s2a = var(a,[],'all');
s2i = var(infiltration,[],'all');
mui = mean(infiltration,'all');
mua = mean(a,'all');
bara = rad.pmu.a;
infiltration_ = flat(bara.*rad.infiltration_enhancement(b)).*flat(h);
s2i_ = var(infiltration_);
%mui_ = mean(infiltration);
%infiltration = infiltration.*flat(h);
%cva_ = std(rad.p.a,[],'all')./mean(rad.p.a,'all');
%cvi  = std(infiltration,[],'all')./mean(infiltration,'all');
%tab{jdx,tdx}.relstd(idx) = cva_./cvi;
tab{jdx,tdx}.relstd(idx) = sqrt(s2a)/mua/(sqrt(s2i)/mui);
tab{jdx,tdx}.relstd_(idx) = (s2i - s2i_)/s2i;
tab{jdx,tdx}.corr_a_I(idx) = corr(a(:),infiltration(:));
tab{jdx,tdx}.corr_h_I(idx) = corr(h(:),infiltration(:));
tab{jdx,tdx}.corr_ie_I(idx) = corr(ie(:),infiltration(:));
tab{jdx,tdx}.bmax(idx) = max(b,[],'all');
a_ =a(:);
ie_ = ie(:);
h_ = h(:);
b_ = b(:);
infiltration_ = infiltration(:);
A = [ones(prod(rad.nx),1),(a_-mean(a_))/std(a_),(ie_-mean(ie_))/std(ie_),(h_-mean(h_))/std(h_)];
c = A \ ((infiltration_-mean(infiltration_))/std(infiltration_))
ccc(idx,1:4,jdx) = c;
tab{jdx,tdx}.fraction_exogenous(idx) = c(2);

A = [ones(prod(rad.nx),1),(a_-mean(a_))/std(a_),(b_-mean(b_))/std(b_),(h_-mean(h_))/std(h_)];
c = A \ ((infiltration_-mean(infiltration_))/std(infiltration_))
ccc_(idx,1:4,jdx) = c;
%(s2i - s2i_)/s2i;
%sqrt(s2a)/mua/(sqrt(s2i)/mui);

else
	ccc(idx,:,jdx) = NaN;
	ccc_(idx,:,jdx) = NaN;
	tab{jdx,tdx}.relstd(idx) = NaN;
	tab{jdx,tdx}.relstd_(idx) = 0;
	kmax = 1;
end
 

clear sp
ff = dir(f);
fa = dir(analysis_str);

if (exist(analysis_str,'file') && (fa.datenum>ff.datenum))
	load(analysis_str,'sp');
else
	disp('reanalyzing');
	for kdx=1 %:kmax
	
	sp(kdx) = Spatial_Pattern();
	sp.opt.suppress_low_frequency_components = 0;
	sp(kdx).L = rad.L;
	sp(kdx).b = b;
	if (~isscalar(a))
		sp(kdx).source = a;
		% sp(kdx).source.S = Sa;
	end % if
	sp(kdx).analyze_grid();
	sp(kdx).fit_parametric_densities();
	sp(kdx).predict_pattern();
	end % for kdx
	save(analysis_str,'sp');
end % else of if exist analyze_str


tab{jdx,tdx}.runtime(idx) = out.runtime(end);
tab{jdx,tdx}.n_step(idx) = out.n_step(end);
tab{jdx,tdx}.adapt_time_step(idx) = rad.opt.adapt_time_step;
tab{jdx,tdx}.r2_bp(idx)  = sp(1).stat.fit.radial.bandpass.stat.goodness.r2;
tab{jdx,tdx}.p_periodic(idx)  = sp(1).stat.p_periodic;
tab{jdx,tdx}.Sc_hat(idx) = sp(1).stat.Sc.radial.hat;
tab{jdx,tdx}.fc_hat(idx) = sp(1).stat.fc.radial.hat;
tab{jdx,tdx}.Sc_con(idx) = sp(1).stat.Sc.radial.con;
tab{jdx,tdx}.fc_con(idx) = sp(1).stat.fc.radial.con;
tab{jdx,tdx}.Sc_bar(idx) = sp(1).stat.Sc.radial.bar;
tab{jdx,tdx}.fc_bar(idx) = sp(1).stat.fc.radial.bar;
tab{jdx,tdx}.Sc_bp(idx)  = sp(1).stat.Sc.radial.bandpass;
tab{jdx,tdx}.fc_bp(idx)  = sp(1).stat.fc.radial.bandpass;
tab{jdx,tdx}.R(idx)      = rad.pmu.R;
tab{jdx,tdx}.coherence(idx) = sp(1).stat.coherence.radial;

c = corr(a(:),b(:));
tab{jdx,tdx}.r2ab(idx)  = c.^2;
for kdx=1 %:kmax
if (0)
S = interp1(sp(kdx).f.r,sp(kdx).S.fit.radial.bandpass,sp(kdx).f.rr,'linear',0);
T = sqrt(S);
b_ = ifft2(T.*fft2(a));

bt=graythresh_scaled(b(:));
bt = b>bt;
q=quantile(b_(:),1-mean(bt(:)));
bt_=b_>q;

ct = corr(bt_(:),bt(:));
c  = sign_to_pearson(ct);
if (1==kdx)
	tab{jdx,tdx}.r2fwb(idx) = c.^2;
else
	tab{jdx,tdx}.r2fab(idx) = c.^2;
end % else of if kdx == 1
end

tab{jdx,tdx}.r2fwb(idx) = NaN;
tab{jdx,tdx}.r2fab(idx) = sp.stat.fit.b_.lin_thresh.r2;

end % for kdx
end % if mean(b)>0

if (0)
btt = bt + 2*bt_;
figure(jdx)
subplot(6,9,idx);
imagesc(btt);
axis equal
axis square;
end % if 0

if (0)
figure(100+jdx)
subplot(6,9,idx);
imagesc(bt_);
axis equal
axis square;
end % if 0

if (0)
subplot(2,2,1);
imagesc(bt);
subplot(2,2,2);
imagesc(bt_);
end % if 0

sp_a(idx,jdx) = sp;

end % for idx

[tab{jdx,tdx},sdx] = sortrows(tab{jdx,tdx},'cva');
bb_ = bb;
bb = bb_(:,:,sdx);
ccc(:,:,jdx) = ccc(sdx,:,jdx);
ccc_(:,:,jdx) = ccc_(sdx,:,jdx);
sp_a(:,jdx) = sp_a(sdx,jdx);

% quick fix
tab{jdx,tdx}.coherence = tab{jdx,tdx}.coherence /1024;
tab{jdx,tdx}.Scr = arrayfun(@(x) x.stat.Sc.radial.bandpass,sp_a);
tab{jdx,tdx}.fcr = arrayfun(@(x) x.stat.fc.radial.bandpass,sp_a);
tab{jdx,tdx}.regr = tab{jdx,tdx}.Scr./tab{jdx,tdx}.fcr;
tab{jdx,tdx}.r2Sr = arrayfun(@(x) x.stat.fit.radial.bandpass.stat.goodness.r2,sp_a);
tab{jdx,tdx}.p_periodic = arrayfun(@(x) x.stat.p_periodic,sp_a);

%set(0,'CurrentFigure',jdx)
%id = round(100*tab{jdx,tdx}.cva/0.2);
for idx=1:length(file_C)
%fdx = find(id == idx)
%if (~isempty(fdx))
try
subplot(6,9,idx);
imagesc(bb(:,:,idx));
axis equal
axis square;
catch e
e
end % catch of try
end % for idx file_C (sa)

%end

% TODO compute and migration rate
% TODO compute and write fraction of exogenous heterogeneity
tabw{jdx,tdx} = [tab{jdx,tdx}(:,'cva'),...
        tab{jdx,tdx}(:,'Scr'), ...
	tab{jdx,tdx}(:,'fcr'), ...
	tab{jdx,tdx}(:,'regr'), ...
	tab{jdx,tdx}(:,'coherence'), ...
	tab{jdx,tdx}(:,'r2Sr'), ...
	tab{jdx,tdx}(:,'fraction_exogenous') ...
];
tabw{jdx,tdx} = round(tabw{jdx,tdx},3,'significant');
writetable(tabw{jdx,tdx},[folder,'/rietkerk-',label{jdx},'.csv']);
writetable(tabw{jdx,tdx},['output/rietkerk-',label{jdx},'.csv']);

end % for jdx % R 

end % for tdx

save('mat/pattern-formation-statistics-quick.mat','tab','ccc','ccc_'); 
else
	load('mat/pattern-formation-statistics-quick.mat')
end % if (0,1)

for idx=1:3
	fh = figure(idx);
	fh.Visible = true;
end

meta = pattern_formation_metadata;
%figure(2)
%clf;
cmap = meta.colororder;
 xlim_ =[-1e-2,0.51];
 bara = 0.2;
for tdx=1:length(theta)
for idx=1:3
	splitfigure([2,3],[2,1],fflag);
	if (1==idx) cla; end
	plot(tab{idx,tdx}.cva,tab{idx,tdx}.r2ab, ...
		'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
	hold on;
	xlim(xlim_);
	axis square
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Goodness of pattern fit $R_{(a,RK)}^2$','interpreter','latex');
	%legend('Spotted','Labyrinthine','Gapped','location','east');
 	%%set(gca,'colororder',meta.colororder);
	ylim([0,1])
	%title('Unfiltered heterogeneity','interpreter','latex');

	splitfigure([2,3],[2,2],fflag);
	if (1==idx) cla; end
	plot(tab{idx,tdx}.cva,tab{idx,tdx}.r2fwb, ...
		'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
	hold on
	xlim(xlim_);
	axis square
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Goodness of pattern fit $R_{(BP(W),RK)}^2$','interpreter','latex');
	%legend('Spotted','Labyrinthine','Gapped','location','east');
 	%set(gca,'colororder',meta.colororder);
	ylim([0,1])
	%title({'Filtered heterogeneity','with unknown spectrum'},'interpreter','latex');

	splitfigure([2,3],[2,3],fflag);
	if (1==idx) cla; end
	plot(tab{idx,tdx}.cva,tab{idx,tdx}.r2fab, ...
		'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
	hold on
	xlim(xlim_);
	ylim([0,1])
	axis square
 	%set(gca,'colororder',meta.colororder);
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Goodness of pattern fit $R_{BP(a),b}^2$','interpreter','latex');
	%legend('Spotted','Labyrinthine','Gapped','location','southeast');
 	%set(gca,'colororder',meta.colororder);
	%title({'Filtered heterogeneity','with known Spectrum'},'interpreter','latex');

	splitfigure([2,3],[2,4],fflag);
	if (1==idx) cla; end
	plot(tab{idx,tdx}.cva,squeeze(ccc(:,2,idx)), ...
		'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
	hold on
	xlim(xlim_);
	ylim([0,1]);
	axis square;
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel({'Fraction of exogeneous heterogeneity'})

	splitfigure([2,3],[2,5],fflag);
	if (1==idx) cla; end
	plot(tab{idx,tdx}.cva, tab{idx,tdx}.coherence, ...
		'o', 'markerfacecolor', cmap(idx,:), ...
		'color','none', 'linewidth', 1, ...
		'markersize', markersize);
	hold on
	ylim([0,1]);
	xlim(xlim_);
	axis square
	xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
	ylabel('Average Coherence $\bar c$','interpreter','latex');
	

 splitfigure([2,3],[10,1],fflag);
 if (1 == idx) cla; end
 plot(tab{1}.cva,tab{idx,1}.Sc_bp.*tab{idx,1}.fc_bp, ...
		'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
 %set(gca,'colororder',meta.colororder);
 hold on;
 xlim(xlim_);
 ylim([0,5.25]);
 axis square
 ylabel('Regularity $S_{rc}/\lambda_c$','interpreter','latex');
 xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
% fit a power law
 A = log(tab{idx,1}.cva).^[0,1];
 reg=tab{idx,1}.Sc_bp.*tab{idx,1}.fc_bp;
 rhs = log(reg);
 fdx = isfinite(rhs) & isfinite(A(:,2));
 c=A(fdx,:) \ rhs(fdx);
 cc(:,idx) =c;
 regp = exp(A*c);
 r2(idx) = 1 - rms(regp(fdx)-reg(fdx)).^2./var(reg(fdx));
if (0)
 %set(gca,'colororderindex',idx);
 plot(tab{1}.cva, exp(A*c)); 
 %set(gca,'colororder',meta.colororder);
end
 splitfigure([2,3],[10,2],fflag);
 if (1 == idx) cla; end
 plot(tab{1}.cva,1./tab{idx,1}.fc_bp, ...
	       'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
 xlim(xlim_);
 ylabel('Wavelength $\lambda_c$','interpreter','latex');
 hold on
 xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
 %set(gca,'colororder',meta.colororder);
 axis square
 if (0)
 A = (tab{1}.cva).^[0,1];
 b = 1./tab{idx,1}.fc_bp;
 fdx = isfinite(b) & isfinite(A(:,2));
 c= A(fdx,:)\b(fdx);
 %set(gca,'colororderindex',idx);
 wavelength_p = A*c;
 plot(tab{1}.cva,wavelength_p);
 end

 splitfigure([2,3],[10,3],fflag);
 if (1 == idx) cla; end
 plot(tab{1}.cva,tab{idx,1}.r2_bp,'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
 xlim(xlim_);
 ylabel('Goodness of density fit $R_{S_r}^2$','interpreter','latex');
 hold on
 xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
 ylim([0,1]);
 %set(gca,'colororder',meta.colororder);
 axis square

 splitfigure([2,3],[10,4],fflag);
 plot(tab{1}.cva,tab{idx,1}.p_periodic, ...
	'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1,'markersize',markersize);
 xlim(xlim_);
 hold on
 axis square

 splitfigure([2,3],[10,5],fflag);
 plot(tab{idx}.cva,tab{idx}.Sc_bp, ...
	'o','markerfacecolor',cmap(idx,:),'color','none','linewidth',1, ...
	'color','none','markersize',markersize);
 xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
 ylabel('Density Maximum $S_{rc}$','interpreter','latex');
 hold on
 xlim(xlim_);
 ylim([0,1.05*max(cellfun(@(x) max(x.Sc_bp),tab))]);
 axis square

 splitfigure([2,3],[10,6],fflag);
 if (1 == idx) cla(); end
 plot(tab{1}.cva,tab{idx}.relstd, ...
	'o','markerfacecolor',cmap(idx,:),'color','none', ...
	'linewidth',1, ...
	'color','none','markersize',1);
 hold on
 xlabel('Exogenous heterogeneity $CV(a)$','interpreter','latex');
 ylabel('Exogenous / Total heterogeneity $CV(a)/CV(I)$','interpreter','latex');
 xlim([-0.01,0.51]);


 end; % for idx
 end; % for tdx

figure(1e4);
for idx=1:3;
	subplot(2,3,idx);
	plotyy(tab{1}.cva,tab{idx}.runtime,tab{1}.cva,tab{idx}.adapt_time_step);
end

figure(1e5);
%plot(NaN(3),'o','color','none','markerfacecolor',cmap)
cmap=meta.colormap;
clf;
 for idx=1:3;
 plot(NaN(1),'o','color','none','markerfacecolor',cmap(idx,:));
 hold on;
 end;
 legend('Spotted','Gapped','Labyrinthine');
 box off;
 axis off;

if (pflag)
	ps = 4;
	pdfprint( 21,'img/rk-2d-isotropic-pattern-fit-unfiltered-heterogeneity.pdf',ps);
	pdfprint( 22,'img/rk-2d-isotropic-pattern-fit-filtered-heterogeneity-assuming-flat-noise-spectrum.pdf',ps);
	pdfprint( 23,'img/rk-2d-isotropic-pattern-fit-filtered-heterogeneity.pdf',ps);
	pdfprint( 24,'img/rk-2d-isotropic-fraction-of-exogenous-heterogeneity.pdf',ps);
	pdfprint( 25,'img/rk-2d-isotropic-spectral-coherence.pdf',ps);
	pdfprint(101,'img/rk-2d-isotropic-regularity',ps);
	pdfprint(102,'img/rk-2d-isotropic-wavelength',ps);
	pdfprint(103,'img/rk-2d-isotropic-density-goodness-of-fit',ps);
	pdfprint(105,'img/rk-2d-isotropic-density-maximum',ps);
 	pdfprint(1e5,'img/pattern-type-legend.pdf',ps);
end

