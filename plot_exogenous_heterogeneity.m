% 2024-12-09 22:32:26.141104355 +0100
if (~exist('pflag','var'))
	pflag = 0;
end
dx =1;
L = 1024;
nx = L/dx;
% TODO use OU
[fx,fy,frr] = fourier_axis_2d(L*[1,1],nx*[1,1]);
L = 1024;

mu = 1;
sd = 0.1;
rng(0);
if (0)
	Sr=[1./fx.^2,1./(fx.^2+1./L.^2)];
	Sr(1)=0;
	Sr=Sr(:,2);
	df = 1/L;
	Sr=Sr./(sum(Sr)*df);
	S2d = 1./(frr.^2 + 1./L.^2);
	T2d = sqrt(S2d);
	e = randn(nx,nx);
	a = ifft2(T2d.*fft2(e));
	a = mu+sd*a/std(a,[],'all');
	fr = fx;
else
	integration_order = 5;
        oversampling_factor_spectral = 2;
        p_window = 0.5;         
        theta = 256;                                     
        [a,C2d,S2d] = geometric_ou_2d_grid_cell_averaged_generate(mu,sd,theta,L*[1,1],nx*[1,1],integration_order,oversampling_factor_spectral,p_window);
	[Sr,fr] = periodogram_radial(S2d,L*[1,1]);
	Sr = Sr.normalized;
end


figure(1e3);
 clf;
 a0=0.2;
 sa=[0.1,0.2,0.5];
 x=linspace(0,0.2*3,1e3)';
 p=[];
 for idx=1:length(sa)
	[pa,pb]  = lognpdf_moment2par(a0,0.2*sa(idx));
	p(:,idx) = lognpdf(x,pa,pb);
	hold on;
 end
 plot(x/a0,p*a0,'linewidth',1);
 if (~pflag)
	[y,x]=ksdensity(a(:),'Function','pdf')
	hold on
	plot(x,y,'--','linewidth',2);
 end

 lh=legend(num2str(cvec(sa)));
 title(lh,'$CV(a)$','interpreter','latex');
 xlabel('Infiltration coefficient $a / \bar a$','interpreter','latex');
 ylabel('Probability density $P{\cdot}\bar a$','interpreter','latex')
set(gca,'colororder',[0,0,0;
        0.8,0,0;
        0,0,0.8]);
ylim([0,4.5])
xlim([0,x(end)/a0]);
axis square

figure(1e3+1)
loglog(fr,Sr,'linewidth',1,'color',[0,0,0.8]);
ylim([1e-3,1e2]);
xlabel('Wavenumber $k_r/(2 \pi)$','interpreter','latex');
ylabel('Radial density $S_{a,r}$','interpreter','latex')
axis square
set(gca,'ytick',10.^(-3:3));
set(gca,'xtick',10.^(-3:3));

a0=0.2;
n=1024;
x=0:n-1;
figure(1e3+2)
%imagesc(x,x,reshape(aa(:,11)/a0,1024*[1,1]));
imagesc(x,x,a/mu);
colormap gray;
axis square;
c=colorbar();
title(c,'$a/\bar a$','interpreter','latex');
caxis(1+5*[-1,1]*sd/mu); %0-aa(:,11))/a0))
xlabel('Position $x$ / m','interpreter','latex');
ylabel('Position $y$ / m','interpreter','latex')
axis square
%set(gca,'xtick',10.^(-3:3));
%set(gca,'ytick',10.^(-3:3));

if (pflag)
	ps = 4;
	pdfprint(1e3+0,'img/bare-soil-infiltration-distribution.pdf',4);
	pdfprint(1e3+1,'img/bare-soil-infiltration-spectral-density.pdf',4);
	pdfprint(1e3+2,'img/bare-soil-infiltration-spatial-map.pdf',4);
end

