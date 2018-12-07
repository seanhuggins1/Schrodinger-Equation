function [outputArg1,outputArg2] = seq_crank(Nspace,Ntime,tau,varargin)
    %
% If called with no arguments, echo a useage line. 
%
% if nargin == 0
%    disp(' ')
%    disp('y=function_skeleton(x,[varargin])')
%    disp(' ')
%    funcOutput=[];
%    return
% end

%
% Check input typing.
%
% if ~isnumeric(funcInput)
%    disp(' ')
%    disp('ERROR(function_skeleton): funcInput must be numeric')
%    funcOutput=[];
%    return
% end

%
% Check that all varargin come in pairs.
%
if mod(length(varargin),2) ~= 0
  disp(' ')
  disp('Error: mis-match (odd number) of vargargin inputs')
  disp(' ')
  funcOutput=[];
  return
end

%
% Set defaults and parse the varargin arguments.
%
method='implicit';
L=200;
plotType='psi';
V=zeros(Nspace, 1);
potential_spikes = [];
WParam=[10 0 0.5];

for j=1:2:length(varargin)

    switch lower(varargin{j})

    case 'potential',
       potential_spikes=varargin{j+1};
    case 'length',
       L=varargin{j+1};

    case 'method',
       method=lower(varargin{j+1});
       
    case 'plottype',
       plotType=lower(varargin{j+1});
       
    case 'wparam',
       WParam =varargin{j+1};
       
    otherwise
       disp(' ')
       disp(sprintf('WARNING: unknown varargin <%s> ignored',varargin{j}))
       disp(' ')

    end

end

%Place spikes of potential in potential array
V(potential_spikes,1) = 1;


sig0 = WParam(1);
x0 = WParam(2);
k0 = WParam(3);
hbar = 1;
m = 1/2;

h = L/(Nspace-1);

%%%%%%%%%%Nspace+ 1 or Nspace?
x = linspace(-L/2,L/2, Nspace);

%Compute the Hamiltonian Matrix
I = eye(Nspace);
for j = 1:Nspace
    %Logical indexing to handle periodic boundary conditions
    jm = j-1;
    jp = j+1;
    if jm == 0
        jm = Nspace;
    end
    if jp == Nspace+1
        jp = 1;
    end
    
    for k = 1:Nspace
        H(j,k) = ((-hbar^2)/m)*(I(jp,k) + I(jm,k) - 2*I(j,k))/(h^2) + V(j,1)*I(j,k);
    end
end





%Initialize the Wave Function 
psi = zeros(length(x), Ntime);
psi(:,1) = (1/sqrt(sig0*(pi^(1/2)))).*exp((1i*k0).*x).*exp((-(x-x0).^2)./(2*(sig0^2)));

%Initialize the probabliity distribution
prob = zeros(length(x), Ntime);
prob(:,1) = psi(:,1).*conj(psi(:,1));


% plot(x, phi(:,1));
% axis([-L/2 L/2 -.5 .5]);
% drawnow;

CN = (inv((I + ((1i*tau)/(2*hbar)).*H)))*(I - ((1i*tau)/(2*hbar)).*H);
FTCS = (I + ((1i*tau)/(2*hbar)).*H);

if strcmp(method,'implicit')
    M=CN;
elseif strcmp(method,'explicit')
    M=FTCS;
    %Grab spectral radius
    r = max(eig(M))
    %Check for stability
    if r > 1
        error('Spectral radius is greater than unity. Please choose a smaller time step for a stable solution.')
        return
    end
        
end



clf
figure(1);

if strcmp(plotType, 'psi')
    for t = 1:Ntime
        psi(:,t+1) = M*psi(:,t);
        plot(x, real(psi(:,t)));
        axis([-L/2 L/2 -.5 .5]);
        xlabel('x');
        ylabel('psi(x,t)');
        drawnow;
    end
elseif strcmp(plotType, 'prob')
    for t = 1:Ntime
        psi(:,t+1) = M*psi(:,t);
        prob(:,t+1) = psi(:,t+1).*conj(psi(:,t+1));
        plot(x, prob(:,t), 'r');
        axis([-L/2 L/2 0 .1]);       
        xlabel('x');
        ylabel('P(x,t)');
        drawnow;
    end
end


