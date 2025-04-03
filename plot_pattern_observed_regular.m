% Tue  7 Mar 21:28:32 CET 2023
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
%% plot an isotropic and an anisotropic natural pattern
%
function sp = plot_pattern_observed_regular(meta) 
	if (nargin<1)
		meta = pattern_formation_metadata();
	end
	pflag = meta.pflag;
	fflag = pflag;
	ps = meta.plotscale;
	fcmap = meta.fcmap;
	field = 'hp';

	file_C = {
		  'patterns/anisotropic_28.26432_11.17174.png', 2.5,
		  'patterns/isotropic_-114.866787_38.385718_mercator_cropped.png', 2.5
		  %'patterns/isotropic_-107.150197_31.320446.png', 3.5
	}
	C = {'bandpass','phase_drift','logn','gamma'};

	
	if (~exist(meta.filename.observed_patterns,'file') || ~meta.reload)
		for idx=1:length(file_C)
		sp(idx) = Spatial_Pattern();
		sp(idx).imread(file_C{idx,1});
		if (1==idx)
		        sp(idx).L   = sp(idx).L([2,1]);
			sp(idx).b = sp(idx).b';
        		sp(idx).msk.b = sp(idx).msk.b';  
		end
 		sp(idx).analyze_grid();
		sp(idx).fit_parametric_densities();
%		sp(idx).b_square = sp(idx).b_square';
%		sp(idx).b = sp(idx).b';
		end
		save(meta.filename.observed_patterns,'-v7.3','sp');
	else
		load(meta.filename.observed_patterns,'sp');
	end

	for idx=1:length(sp)
	sp(idx).opt.scalefield = field;
	
	if (sp(idx).stat.isisotropic)
		printf('Src/lc: %g\n',sp(idx).stat.fc.radial.(field)*sp(idx).stat.Sc.radial.(field));
		printf('lc: %g\n',1./sp(idx).stat.fc.radial.(field));
		printf('R2 %g\n',sp(idx).stat.fit.radial.bandpass.stat.goodness.r2);
	else
		printf('Sxc/lc: %g\n',sp(idx).stat.fc.x.(field)*sp(idx).stat.Sc.x.(field));
		printf('Syc/lc: %g\n',sp(idx).stat.fc.x.(field)*sp(idx).stat.Sc.y.(field));
		printf('lc: %g\n',1./sp(idx).stat.fc.x.(field));
		printf('R2 %g\n',sp(idx).stat.fit.x.phase_drift.stat.goodness.r2);
	end
	printf('p-periodic %g\n',sp(idx).stat.p_periodic)

	% plot pattern
	splitfigure([2,3],[idx,1],fflag);
	sp(idx).plot('b');
	colormap(flipud(gray));
	caxis(0.5+0.0125*[-1,1]);
	axis(min(sp(idx).L)*sp(idx).stat.fc.radial.hp*[0,1,0,1])
	colormap(fcmap(256))

	% plot 2d periodogram
	splitfigure([2,3],[idx,2],fflag);
	sp(idx).S.rot.hp = sp(idx).S.rot.hp';
	sp(idx).plot('S.rot.hp');
	if (sp(idx).stat.isisotropic)
		caxis([0,1.5*sp(1).stat.Sc.radial.hp*sp(1).stat.fc.radial.hp/6])
	else
		caxis([0,12*sp(1).stat.Sc.radial.hp*sp(1).stat.fc.radial.hp/6])
	end
	axis(file_C{idx,2}*[-1,1,-1,1])
%	caxis([0,])
	colormap(fcmap(16))

	splitfigure([2,3],[idx,3],fflag);
	sp(idx).R.rot.hp = sp(idx).R.rot.hp';
	sp(idx).plot('R.rot.hp');
	xlim([-2.5,2.5]);	
	ylim([-2.5,2.5]);
	if (idx==1)	
	%caxis([-0.3355,0.5]);
	caxis([-0.4,1])
	colormap(fcmap(14));
	else
	%caxis([-0.0512,0.125]);
	%caxis([-0.13,0.97]);
	xlim(2.5*[-1,1]);
	ylim(2.5*[-1,1]);
	caxis([-0.1,1.0]);
	colormap(fcmap(11));
	end
	cbh = colorbar();
	title(cbh,'$\hat R$','interpreter','latex');
	

	% density along primary axis
	splitfigure([2,3],[idx,4],fflag);
	cla();
	if (~sp(idx).stat.isisotropic)
		sp(idx).plot('S.rot.x.hp','linewidth',1);
		hold on
		fdx = sp(idx).f.x>0;
		%S = sp(idx).w.x(fdx).*sp(idx).S.rot.x.phase_drift;
		%S = S/spectral_density_area(sp(idx).f.x(fdx),S);
		%plot(sp(idx).f.x(fdx)/sp(idx).stat.fc.x.hp,S*sp(idx).stat.fc.x.hp,'linewidth',1);
		sp(idx).plot('S.rot.x.phase_drift','linewidth',1);
		legend('empirical','PNI-fit');
	else
		sp(idx).plot('S.radial.hp','linewidth',1);
		hold on
		plot(sp(idx).f.r/sp(idx).stat.fc.radial.hp,sp(idx).S.radial.bandpass*sp(idx).stat.fc.radial.hp,'linewidth',1)
		legend('empirical','BP-fit');
	end
	xlim(file_C{idx,2}*[0,1])
	axis square

	% density along secondary axis
	splitfigure([2,3],[idx,5],fflag);
	cla();
	if (~sp(idx).stat.isisotropic)
		sp(idx).plot('S.rot.y.hp','linewidth',1);
		hold on;
		sp(idx).plot('S.rot.y.phase_drift_parallel','linewidth',1);
		%plot(sp(idx).f.y(fdx)/sp(idx).stat.fc.x.hp,sp(idx).S.rot.y.phase_drift_parallel(fdx)*sp(idx).stat.fc.x.hp,'linewidth',1)
		xlim(file_C{idx,2}*[0,1])
		%legend('empirical','PNI-fit');
	else
		sp(idx).plot('S.rot.angular.hp','linewidth',1);
		hold on;
		ylim([0,0.5]);
		hline(1/pi,'color','r','linewidth',1)
		%legend('empirical','BP-fit');
	end
	axis square
if (0)
	splitfigure([2,3],[idx+20,2],fflag);
	cla
	plot(NaN,NaN);
	hold on

	fdx = sp(idx).f.x>=0;	

	if (~pflag)

	err = [];
	leg_C{1} = 'orig';
 	for jdx=1:length(C)
	splitfigure([2,3],[idx,4],fflag);
	S = sp(idx).S.rot.x.(C{jdx})(fdx);
	S = sp(idx).w.x(fdx).*S;
	S = S/spectral_density_area(f(fdx),S);
	plot(sp(idx).f.x(fdx)/sp(idx).stat.fc.x.hp,S*sp(idx).stat.fc.x.hp)

	% plot residual
	%splitfigure([2,3],[idx+20,2],fflag);
	%plot(sp(idx).f.x(fdx)/sp(idx).stat.fc.x.hp,(S-sp(idx).S.rot.x.hp(fdx))); %*sp(idx).stat.fc.x.hp))

		err(1,jdx) = sp(idx).stat.fit.x.(C{jdx}).stat.goodness.rmse; %./sp(idx).stat.fc.x.hp;
		err(2,jdx) = sp(idx).stat.fit.x.(C{jdx}).stat.goodness.cdf_l1; %./sp(idx).stat.fc.x.hp;
		err(3,jdx) = sp(idx).stat.fit.x.(C{jdx}).stat.goodness.mise_cramer; %./sp(idx).stat.fc.x.hp;
		leg_C{jdx+1} = [C{jdx},num2str([err(1,jdx),err(2,jdx) err(3,jdx)])];
	end
	%splitfigure([2,2],[idx,3],fflag);
	%plot(sp(idx).f.x/sp(idx).stat.fc.x.hp,sp(idx).w.x);
	err
	legend(leg_C{:})
	else
		splitfigure([2,3],[idx,4],fflag);
		plot(sp(idx).f.x(fdx)/sp(idx).stat.fc.x.hp,sp(idx).S.rot.x.phase_drift(fdx)*sp(idx).stat.fc.x.hp,'linewidth',1)	
		legend('empirical','PNI');
	end
	xlim(file_C{idx,2}*[0,1]);
	axis square

	splitfigure([2,3],[idx,5],fflag);
	plot(sp(idx).f.y(fdx)/sp(idx).stat.fc.x.hp,sp(idx).S.rot.y.phase_drift_parallel(fdx)*sp(idx).stat.fc.x.hp,'linewidth',1)
	xlim([0,2.5]);
	axis square		

	% radial density
%	splitfigure([2,3],[idx,6],fflag);
%	cla();
%	sp(idx).plot('S.radial.hp','linewidth',1);
%	hold on

	% residual
%	splitfigure([2,3],[idx+20,3],fflag);
%	cla
%	plot(NaN,NaN);
%	hold on


	if (~pflag)
 	for jdx=1:length(C)
		splitfigure([2,3],[idx,6],fflag);
		plot(sp(idx).f.r/sp(idx).stat.fc.radial.hp,sp(idx).S.radial.(C{jdx})*sp(idx).stat.fc.radial.hp)
		splitfigure([2,3],[idx+20,3],fflag);
	%plot(sp(idx).f.x(fdx)/sp(idx).stat.fc.x.hp,(sp(idx).S.radial.hp(fdx)-sp(idx).S.radial.hp(fdx))); %*sp(idx).stat.fc.x.hp))
		err(1,jdx) = sp(idx).stat.fit.radial.(C{jdx}).stat.goodness.cdf_l1./sp(idx).stat.fc.radial.hp;
		err(2,jdx) = sp(idx).stat.fit.radial.(C{jdx}).stat.goodness.mise_cramer./sp(idx).stat.fc.radial.hp;
		leg_C{jdx+1} = num2str([err(1,jdx),err(2,jdx) err(3,jdx)]);
	end
	err
	legend(leg_C);
	else
		splitfigure([2,3],[idx,6],fflag);
	plot(sp(idx).f.r/sp(idx).stat.fc.radial.hp,sp(idx).S.radial.bandpass*sp(idx).stat.fc.radial.hp,'linewidth',1)
		legend('empirical','bandpas');
	end
	xlim(file_C{idx,2}*[0,1]);
	axis square
end

	if (pflag)
		f = basename(file_C{idx,1});
		f = f(1:end-4);
		pdfprint(10*idx+1,['img/',f,'-pattern.pdf'],ps);
		pdfprint(10*idx+2,['img/',f,'-S2d.pdf'],ps);
		pdfprint(10*idx+3,['img/',f,'-correlogram.pdf'],ps);
if (sp(idx).stat.isisotropic)
		pdfprint(10*idx+4,['img/',f,'-Sr.pdf'],ps);
		pdfprint(10*idx+5,['img/',f,'-Sa.pdf'],ps);
else
		pdfprint(10*idx+4,['img/',f,'-Sx.pdf'],ps);
		pdfprint(10*idx+5,['img/',f,'-Sy.pdf'],ps);
end
		%pdfprint(10*idx+6,['img/',f,'-Sr.pdf'],ps);
	end
	end % for idx
end

