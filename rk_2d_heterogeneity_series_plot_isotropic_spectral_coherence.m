% 2025-03-26 14:26:51.476953306 +0100
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

tab_filename = 'mat/pattern-formation-statistic-isotropic.mat';

load(tab{1}.name{11});
% pattern
 b    = rad.extract2(y(:,2));
% heterogenous infiltration coefficient
 a    = reshape(rad.p.a,rad.nx);
% discrete spectrum of pattern
 fb   = fft2(b-mean(b(:)));
% discrete spectrum of heterogeneity
 fa   = fft2(a-mean(a,'all'));
% peropdpgrams and cross-periodogram
 hSab = fa.*conj(fb);
 hSa  = fa.*conj(fa);
 hSb  = fb.*conj(fb);

% estimated densities and cross-densities
if (1)
 nf = 1.5;
 hSa = ifftshift(gaussfilt2(fftshift(hSa),nf));
 hSab = ifftshift(gaussfilt2(fftshift(hSab),nf));
 hSb = ifftshift(gaussfilt2(fftshift(hSb),nf));
end

% radial densities and cross-densities
 [Sra,fr]  = periodogram_radial(hSa,rad.L);
 Srab = periodogram_radial(hSab,rad.L);
 Srb  = periodogram_radial(hSb,rad.L);

 [Srbc,mdx]=max(Srb.normalized);
 fc = fr(mdx);
 coherence=abs(Srab.mu).^2./(Srb.mu.*Sra.mu);

 figure(1);
 clf();
 plot(fr/fc,Srb.normalized*fc,'linewidth',1);
 ylim([0,1.2]);
 ylabel('Spectral density $S_r/\lambda_c$','interpreter','latex');
 yyaxis right
 plot(fr/fc,coherence,'linewidth',1);
 ylim([0,1.2]);
 xlim([0,3]);
 ylabel('Spectral coherence $|S_{ber}|^2/(S_{br} S_{er})$','interpreter','latex');
 xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
 axis square

 cva = rad.pss.a/rad.pmu.a;

if (pflag)
	ps = 4;
	name = sprintf('img/rk-isotropic-spectral-coherence-cva-%g.pdf',cva);
	pdfprint(1,name,ps);
end 

