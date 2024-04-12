% Tue  2 Apr 16:05:27 CEST 2024
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
%% plot entropy of several distributions deoending on regularity
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% distributions to illustrate
d = {'lognormal','BP','PNI'}
%d = {'normal','gamma','BP','PNI','lognormal'}

% spatial extend of pattern
L = 500;
% number of grid points, chosen that resolution in frequency and reals space are equal
n = L^2;
% spectral resolution 
df = 1./L;
% characteristic frequency 1/lambda_c
fc = 1;
% regularity Sc/lc
reg = logspace(-1,1,10);
% axis in frequency space
fx = fourier_axis(L,n)';

% entropy
e =[];
% distribution parameter
aa = [];
bb = [];
figure(1);
clf();
Sc_ = [];
% for each regularity
for idx=1:length(reg)
    % for each distribution
    for ddx=1:length(d)
	switch (d{ddx})
	case {'normal'}
		% determine parameters
		[mu,sd] = normpdf_wrapped_mode2par(fc,reg(idx));
		aa(idx,ddx) = mu;
		bb(idx,ddx) = sd;
		% density
		S = normpdf_wrapped(fx,mu,sd);
	case {'gamma'}
		% determine parameters
		[a,b] = gamma_mode2par(fc,reg(idx));
		aa(idx,ddx) = a;
		bb(idx,ddx) = b;
		% density
		S = gampdf(fx,a,b);
	case {'BP'}
		% determine parameters
		[a,b] = bandpass2d_continuous_pdf_mode2par(fc,reg(idx),[],L,n);
		aa(idx,ddx) = a;
		bb(idx,ddx) = b;
		[fc_(idx,ddx),Sc_(idx,ddx)] = bandpass2d_continuous_pdf_mode(a,b,L,n);
		% density
		S = bandpass2d_continuous_pdf_hankel(L,n,a,b);
		%[fc_,Sc_] = bandpass2d_continuous_pdf_mode(fc,b)l
		%S = bandpass2d_continuous_pdf_exact(fx,a,b);
	case {'lognormal'}
		% determine parameters
		[a,b] = logn_mode2par(fc,reg(idx));
		aa(idx,ddx) = a;
		bb(idx,ddx) = b;
		% density
		S = lognpdf(fx,a,b);
	case {'PNI'}
		% determine parameters
		[a,b] = phase_drift_pdf_mode2par(fc,reg(idx),true);
		%[fc_(idx,ddx),Sc_(idx,ddx)] = phase_drift_pdf_mode(a,b);
		aa(idx,ddx) = a;
		bb(idx,ddx) = b;
		% density
		S = phase_drift_pdf(fx,a,b);
	end
	% normalize the density
	S = S(1:end/2);
	S = S/(sum(S)*df);
	fdx = (S>0);
	% entropy
	e(idx,ddx) = -sum(S(fdx).*log2(S(fdx)))*df;	
	% plot distributions
	if (idx == idx)
	figure(1)
	nrow = round(sqrt(length(reg)*1/2));
	ncol = ceil(length(reg)/nrow);
	subplot(nrow,ncol,idx)
	if (1==ddx)
	cla
	end
	plot(fx(1:end/2),S)
	hold on
	xlim([0,3])
	end
    end % for jdx
end % for idx

splitfigure([2,2],[2,1],fflag);
semilogx(reg,aa);
ylabel('Parameter a');

splitfigure([2,2],[2,2],fflag);
semilogx(reg,bb);
ylabel('Parameter bb');

splitfigure([2,2],[2,3],fflag);
semilogx(reg,e,'linewidth',1);
legend(d{:},'location','northeast');
xlabel('Regularity $S_c/\lambda_c$','interpreter','latex');
ylabel('Entropy $-\frac{1}{2 \pi} \int S \log_2(S) \mathrm{d} k$','interpreter','latex');
set(gca,'colororder',[0,0,0;0.75,0,0;0,0,0.75]);

if (pflag)
	ps = 3.5;
	pdfprint(23,'img/entropy-vs-regularity.pdf',ps);
end

