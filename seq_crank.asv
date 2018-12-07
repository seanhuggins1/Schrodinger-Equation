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
Method='Implicit';
L=200;
PlotType='Psi';
Potential=[];
WParam=[10 0 0.5];

for j=1:2:length(varargin)

    switch lower(varargin{j})

    case 'verbose',
       verbose=varargin{j+1};

    case 'plot',
       doplot=varargin{j+1};

    otherwise
       disp(' ')
       disp(sprintf('WARNING: unknown varargin <%s> ignored',varargin{j}))
       disp(' ')

    end

end


sig0 = WParam(1);
x0 = WParam(2);
k0 = WParam(3);
hbar = 1;
m = 1/2;
h = L/(Nspace-1);




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
        H(j,k) = ((-hbar^2)/m)*(I(jp,k) + I(jm,k) - 2*I(j,k))/(h^2);
    end
end


%Compute the Crank-Nicholson matrix
CN = (inv((I + ((1i*tau)/(2*hbar)).*H)))*(I - ((1i*tau)/(2*hbar)).*H);

%%%%%%%%%%Nspace+ 1 or Nspace?
x = linspace(-L/2,L/2, Nspace);
%Initialize the Wave Function 
phi(:,1) = (1/sqrt(sig0*(pi^(1/2)))).*exp((1i*k0).*x).*exp((-(x-x0).^2)./(2*(sig0^2)));
plot(x, phi(:,1));
axis([-L/2 L/2 -.5 .5]);
drawnow;


for t = 1:Ntime
    phi(:,t+1) = CN*phi(:,t);
    plot(x, phi(:,t));
    axis([-L/2 L/2 -.5 .5]);
    drawnow;
    %pause(0.00001);
end



end

