function seq_crank(Nspace,Ntime,tau,varargin)
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
       potential_spikes = round(varargin{j+1});
       
       %Check validity of potential spikes input
       if ~isnumeric(potential_spikes) || ~isreal(potential_spikes) ...
               || any(isinf(potential_spikes)) || any(isnan(potential_spikes)) ...
               || any(potential_spikes <= 0)  || any(potential_spikes > Nspace)
           error('Potential spikes must be finite, positive integers within the grid.'); 
       end
    case 'length',
       L=varargin{j+1};
       
       %Check validity of length input
       if ~isnumeric(L) || ~isreal(L) || isinf(L) || isnan(L) || L <= 0
           error('Length must be a finite, positive, real number.'); 
       end
           
       
       
       
    case 'method',
       method=lower(varargin{j+1});
           
       allowedMethods = {'implicit','explicit'};
       
       %Check validity of method type input
       if ~any(strcmp(method,allowedMethods)) 
           error('Method type not recognized. Please choose "implicit" for Crank-Nicholson or "explicit" for FTCS.');      
       end           
    case 'plottype',
       plotType=lower(varargin{j+1});
       
       %Allowed strings for plot type
       allowedPlots = {'psi','prob','party'};
       
       %Check validity of plot type input
       if ~any(strcmp(plotType,allowedPlots))
           error('Plot type not recognized. Please choose "psi" to plot the wave function or "prob" to plot probability density and total probability.');      
       end       
       
    case 'wparam',
       WParam = varargin{j+1};  
       if length(WParam) ~= 3;
           error('WParam must be a vector with 3 parameters [sigma0 x0 k0].');   
       end
       
       
    otherwise
       disp(' ')
       disp(sprintf('WARNING: unknown varargin <%s> ignored',varargin{j}))
       disp(' ')

    end

end

%Place spikes of potential in potential array
V(potential_spikes,1) = 1;

%Grab WParam constants from inputs
sig0 = WParam(1);
x0 = WParam(2);
k0 = WParam(3);

%Set important constants
hbar = 1;
m = 1/2;

%Set the spatial grid spacing 
h = L/(Nspace-1);



%Compute the Hamiltonian Matrix
I = eye(Nspace);
for j = 1:Nspace
    %Logical indexing to handle periodic boundary conditions
    jm = j-1;
    jp = j+1;
    %If our left end is out of range, set it to the right end
    if jm == 0
        jm = Nspace;
    end
    %If our right end is out of range, set it to the left end
    if jp == Nspace+1
        jp = 1;
    end
    
    %Build the tridiagonal Hamiltonian with periodic boundary conditions
    for k = 1:Nspace
        H(j,k) = ((-hbar^2)/m)*(I(jp,k) + I(jm,k) - 2*I(j,k))/(h^2) + V(j,1)*I(j,k);
    end
end


%Crank-Nicholson Scheme
if strcmp(method,'implicit')
    %Calculate the matrix needed for forward time integration
    M=(inv((I + ((1i*tau)/(2*hbar)).*H)))*(I - ((1i*tau)/(2*hbar)).*H);
    % NO STABILITY CHECK REQUIRED. CRANK NICHOLSON IS UNCONDITIONALLY
    % STABLE
    
%FTCS Scheme
elseif strcmp(method,'explicit')
    %Calculate the matrix needed for forward time integration
    M = (I - ((1i*tau)/(hbar)).*H);
    %Grab spectral radius
    r = max(abs(eig(M)));
    %Check for stability
    if r > 1
        error('Spectral radius is greater than unity. Please choose a smaller time step for a stable solution.')
        %DO NOT PERFORM INTEGRATION IF UNSTABLE. EXIT CODE.
        return
    end
        
end

%Initialize x spatial array
x = linspace(-L/2,L/2, Nspace);

%Initialize the Wave Function 
psi = zeros(length(x), Ntime);
psi(:,1) = (1/sqrt(sig0*(pi^(1/2)))).*exp((1i*k0).*x).*exp((-(x-x0).^2)./(2*(sig0^2)));

%Clear whatever figure is open 
clf
hold off
%Open up a new figure for plotting
figure(1);

%If the user chooses to plot psi
if strcmp(plotType, 'psi')  
    %Integrate from 1 to the total time
    for t = 1:Ntime
        
        %Calculate forward time psi using matrix multiplication
        psi(:,t+1) = M*psi(:,t);
        
        %Plot psi at the current time
        plot(x, real(psi(:,t)), 'b');
        
        %Set the limits of the axes
        axis([-L/2 L/2 -.5 .5]);
        
        %Label the axes
        xlabel('x');
        ylabel('psi(x,t)');
        
        %Draw the current plot
        drawnow;
    end
%If the user chooses to plot probability
elseif strcmp(plotType, 'prob')
    %Initialize the probabliity distribution
    prob = zeros(length(x), Ntime);
    prob(:,1) = psi(:,1).*conj(psi(:,1));

    %Initialize the total probability
    totalProb = zeros(length(x), Ntime);
    totalProb(:,1) = trapz(x, prob(:,1));
    
    %Integrate from 1 to the total time
    for t = 1:Ntime
        %Calculate forward time psi using matrix multiplication
        psi(:,t+1) = M*psi(:,t);
        
        %Calculate probablility density
        prob(:,t+1) = psi(:,t+1).*conj(psi(:,t+1));
        
        %Calculate the total probability using the inner product 
        totalProb(:,t+1) = trapz(x, prob(:,t+1));
        
        %Plot the probability density at the current time
        plot(x, prob(:,t), 'r');
        
        %Set the axes
        axis([-L/2 L/2 0 .1]);  
        
        %Label the axes
        xlabel('x');
        ylabel('P(x,t)');
        
        %Draw the current plot
        drawnow;
    end
    %Open up a new figure
    figure(2);
    
    %Make the time array again for plotting total probability
    t = 1:Ntime;
    
    %Plot the base 10 logarithm of the absolute value of 1 - total
    %probability as a function of time, to show the difference between the
    %expected value of probability vs. the observed value.
    plot(t, log10(abs(1-totalProb(:,t))), 'k');
    
    %Label the axes
    xlabel('t');
    ylabel('log10 |1 - |P(x,t)|^2|');
%If the user chooses to party
elseif strcmp(plotType, 'party')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'Color','k')
    
    [y, Fs] = audioread('boogiewonderland.mp3');
    sound(y, Fs, 16);
    %Integrate from 1 to the total time
    for t = 1:Ntime
        
        %Calculate forward time psi using matrix multiplication
        psi(:,t+1) = M*psi(:,t);
        
        %Plot psi at the current time
        plot(x, real(psi(:,t)), 'LineWidth', 2 + 5*abs(cos(0.5*t)), 'color', [ abs(cos(0.25*t)) (1 - abs(sin(0.5*t))) abs(sin(0.1*t))]);
        
        %Party Mode
        set(gca,'Color','k')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        %plot(x, V(x))
        txt = 'PARTY ON!!!';
        txt = text((sin(0.1*t)*L)/4,0.3*cos(0.1*t),txt, 'color', [round(abs(cos(0.5*t))) round(abs(sin(0.25*t))) round(abs(cos(0.75*t)))])
        txt.FontSize = 40 + 10*sin(0.2*t);
        txt.Rotation = t
        txt.HorizontalAlignment = 'center';
        
        %Set the limits of the axes
        axis([-L/2 L/2 -.5 .5]);
        
        %Draw the current plot
        drawnow;
    end
end

%End Party
clear sound;

