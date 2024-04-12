% Mon 25 Mar 13:31:12 CET 2024
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
%% illustrate, that the anisotropic RK-model responds to noise by
%% phase noise integration

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

% reset random number generator for reproducibility
rng(0);
% spatial extent
L=10;
% number of grid points
n=1e3;
% spatial axis
x=(0:n-1)'*L/n;
% spatial resolution
dx=L/n;

% number of patterns to generate
m=1e4;

% id of example pattern to plot
id=11;

% characteristic frequency fc = 1/\lambda_c
fc=1;
% maximum of the spectral density 
Sc = 2.5;
% parameters of the phase_drift distribution
[f0,s] = phase_drift_pdf_mode2par(fc,Sc)

% noise along x for m-patterns
e = randn(n,m);
% deterministic frequency
omega=1;
% deterministic phase
phase = omega*x;
% phase deviation by random walk
dp = s*cumsum(e)*sqrt(dx);
% perturbe phase
phase = phase + dp; 

% quantiles of the distribution of the phase
q = (x + norminv([0.16,0.84],0,s*sqrt(x)));
%quantile(phase(end,:),[0.16,0.84])
%q(end,:)

 splitfigure([2,2],[1,1],fflag);
 cla
 errorarea2(x,[q(:,1),x,q(:,2)],[1,1,1]/2);
 hold on;
 set(gca,'colororderindex',1)
 plot(x,[x,phase(:,id)]);
 ylim([0,6]);
 xlabel('Distance $x/\lambda_c$','interpreter','latex');
 ylabel('Phase \phi/(2 \pi)');
 xlim([0,5]);
 legend('Deterministic','Stochastic','location','northwest');

 splitfigure([2,2],[1,2],fflag);
 cla
 b=0.5*(1+-cos(2*pi*f0*phase(:,8)));
 b2 = 0.5*(1 + -cos(2*pi*f0*x));
 b2 = b2/rms(b2);
 b = b/rms(b);
% area(x,b,'facecolor',[0,0.5,0],'linewidth',1);
 plot(x,b2,'k','linewidth',1)
 hold on
 plot(x,b,'r','linewidth',1);
% plot(x,b/rms(b);
 xlabel('Distance x/\lambda_c');
 ylabel('Biomass b/rms(b)');
 xlim([0,5])
 xi = (0.5:4.5)';
 qi = xi + norminv([0.16,0.84],0,s*sqrt(xi));
 %vline(qi)
 phase = 1.7*ones(size(qi));
 line(qi',phase','color','k')
 w = 0.1;
 for idx=1:size(qi,1)
	line(qi(idx,1)*[1,1],1.7+w/2*[-1,1],'color','k')
	line(qi(idx,2)*[1,1],1.7+w/2*[-1,1],'color','k')
 end

if (pflag)
	ps = 3.5;
	pdfprint(11,'img/phase-noise-integration-illustration-phase.pdf',ps);
	pdfprint(11,'img/phase-noise-integration-illustration-phase.png',ps);
	pdfprint(12,'img/phase-noise-integration-illustration-pattern.pdf',ps);	
end

