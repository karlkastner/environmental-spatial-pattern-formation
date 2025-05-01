% 2024-12-09 22:32:26.141104355 +0100
if (~exist('pflag','var'))
	pflag = 0;
end
	% spatial extent
	L = 1024;
	% spatial resolution
	dx =1;
	% number of grid points
	nx = L/dx;
	% spectral axes
	[fx,fy,frr] = fourier_axis_2d(L*[1,1],nx*[1,1]);

	% mean of exogenous heterogeneity
	mu_a = 0.2;
	% coefficient of variation of exogenous heterogeneity
	cv_a = 0.1;
	% standard deviation of exogenous heterogeneity
	sd_a = cv_a*mu_a;
	% correlation length of heterogeneity
        theta = L/4;

	% spatial oversampling factor
	m_spatial = 5;
	% spectral oversampling factor
        m_spectral = 3;
	% relative window size
        p_window = 0.5;

	% reset random number generator for reproducibility
	rng(0);

	% generate heterogeneity
        [a,C2d,S2d] = geometric_ou_2d_grid_cell_averaged_generate(mu_a,sd_a,theta,L*[1,1],nx*[1,1],m_spatial,m_spectral,p_window);
	% radial periodogram
	[Sr,fr] = periodogram_radial(S2d,L*[1,1]);
	Sr = Sr.normalized;


	figure(1e3);
	clf;
	cv_a_=[0.1,0.2,0.5];
	xa=linspace(0,0.2*3,1e3)';
	p=[];
	for idx=1:length(sd_a)
		[pa,pb]  = lognpdf_moment2par(mu_a,mu_a*cv_a_(idx));
		p(:,idx) = lognpdf(xa,pa,pb);
		hold on;
	end
	plot(xa/mu_a,p*mu_a,'linewidth',1);
	if (~pflag)
		[y,x]=ksdensity(a(:),'Function','pdf');
		hold on
		plot(x,y,'--','linewidth',2);
	end

	lh=legend(num2str(cvec(sd_a)));
	title(lh,'$CV(a)$','interpreter','latex');
	xlabel('Infiltration coefficient $a / \bar a$','interpreter','latex');
	ylabel('Probability density $P{\cdot}\bar a$','interpreter','latex');
	set(gca,'colororder',[0,0,0;
	0.8,0,0;
	0,0,0.8]);
	ylim([0,4.5]);
	xlim([0,xa(end)/mu_a]);
	axis square

	figure(1e3+1);
	loglog(fr,Sr,'linewidth',1,'color',[0,0,0.8]);
	ylim([1e-3,1e2]);
	xlabel('Wavenumber $k_r/(2 \pi)$','interpreter','latex');
	ylabel('Radial density $S_{a,r}$','interpreter','latex');
	axis square
	set(gca,'ytick',10.^(-3:3));
	set(gca,'xtick',10.^(-3:3));

	a0=0.2;
	n=1024;
	x=0:n-1;
	figure(1e3+2);
	%imagesc(x,x,reshape(aa(:,11)/a0,1024*[1,1]));
	imagesc(x,x,a/mu_a);
	colormap gray;
	axis square;
	c=colorbar();
	title(c,'$a/\bar a$','interpreter','latex');
	caxis(1+5*[-1,1]*sd_a/mu_a); %0-aa(:,11))/a0))
	xlabel('Position $x$ / m','interpreter','latex');
	ylabel('Position $y$ / m','interpreter','latex');
	axis square
	%set(gca,'xtick',10.^(-3:3));
	%set(gca,'ytick',10.^(-3:3));

	if (pflag)
		ps = 4;
		pdfprint(1e3+0,'img/bare-soil-infiltration-distribution.pdf',4);
		pdfprint(1e3+1,'img/bare-soil-infiltration-spectral-density.pdf',4);
		pdfprint(1e3+2,'img/bare-soil-infiltration-spatial-map.pdf',4);
	end

