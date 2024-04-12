% Thu 19 Oct 12:31:30 CEST 2023
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
%% schematic figure of desities for regular and irregular patterns
%
function plot_schematic_pattern_formation(meta)
if (nargin()<1)
	meta = pattern_formation_metadata();
end
pflag = meta.pflag;
fflag = pflag;
fcmap = meta.fcmap;
ps    = meta.plotscale;

% spatial extent
L = 50;
% number of grid points
n = L^2;
% characteristic frequency fc = 1/\lambda_c
fc = 1;
% maximum of density
Sc0 = 1.0;
Sc = 1.5;
% density parameters
[par] = bandpass2d_continuous_pdf_mode2par(fc,Sc)
[par1] = bandpass2d_continuous_pdf_mode2par(fc,Sc0)
par = [1/par(1),par(2)];

% spatial axis
x = linspace(0,L,n);

% reset random number generator for reproducibility
rng(0);

%[fx,fy,ffr] = fourier_axis(L*
% spatial noise
e = randn(n);

type_C = {'periodic','regular','irregular'};
nt = length(type_C)
for idx=1:nt
% determine densities
switch(type_C{idx})
case{'regular'}
	S = bandpass2d_discrete_pdf(L*[1,1],n*[1,1],par(1),par(2));
%	[Sr,fr] = periodogram_radial(S,L*[1,1]);
	T = sqrt(S);
	b = ifft2(T.*fft2(e));
case{'periodic'}
	b=generate_isotropic_pattern(fc,n*[1,1],L*[1,1],0);
	S = fft2(b-mean(b,'all')).^2;
case{'irregular'}
	% L40, n 2*L^2 0.065
	p2 = par1(2);
	p1 = 0.065;
	% 1.18 0.85
	p1 = 0.073;
	S = lowpass2d_discrete_pdf(L*[1,1],n*[1,1],p1,p2);
%	[Sr,fr] = periodogram_radial(S,L*[1,1]);
	T = sqrt(S);
	b = ifft2(T.*fft2(e));
end
[Sr,fr] = periodogram_radial(S,L*[1,1]);

splitfigure([2,3],[1,idx],fflag);
%q = quantile(b,0.66);
imagesc(x,x,b);
q = quantile(b(:),[0.7,0.8]);
caxis(q)
axis(10*fc*[0,1,0,1])
%axis equal
axis square
if (idx<3)
	xlabel('$x/\lambda_c$','interpreter','latex');
	ylabel('$y/\lambda_c$','interpreter','latex');
else
	xlabel('$x/\lambda_l$','interpreter','latex');
	ylabel('$y/\lambda_l$','interpreter','latex');
end
colormap(fcmap(256));

splitfigure([2,3],[1,idx+3],fflag);
%plot(fr,Sr.normalized/Sr.normalized(1))
if (idx>1)
	plot(fr/fc,fc*Sr.normalized,'linewidth',1)
else
	plot(fr/fc,Sr.normalized/L,'linewidth',1)
end
switch (idx)
case{1}
	ylabel('$S/L$','interpreter','latex')
	xlabel('$k/k_c$','interpreter','latex');
	ylim([0,1+sqrt(eps)]);
case{2}
	ylabel('$S/\lambda_c$','interpreter','latex')
	xlabel('$k/k_c$','interpreter','latex');
	ylim([0,Sc+0.05]);
case{3}
	ylabel('$S/\lambda_l$','interpreter','latex')
	xlabel('$k/k_c$','interpreter','latex');
end
xlim([0,3.5])
axis square
%vline(1); hline(0.5)
figure(2);
if (1==idx) clf; end
loglog(fr,Sr.normalized/max(Sr.normalized))
hold on

end
if (pflag)
	for idx=1:3
		pdfprint(10+idx,['img/pattern-synthetic-',type_C{idx},'.pdf'],ps);
		pdfprint(13+idx,['img/density-synthetic-',type_C{idx},'.pdf'],ps);
	end
end
