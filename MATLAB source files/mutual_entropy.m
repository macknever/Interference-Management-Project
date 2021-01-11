
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  This snippet show the calculation of the mutual information 
%%%  random variable X go through a Guanssian channel(N(0,delta^2))
%%%  We get Y. SO,
%%%  Y = mod(X+Z,2), Z ~ N(0,delta^2)
%%%  X is uniform in[+p,-p]
%%%  The question is I(X;Y) = ??
%%%  
%%%
%%%  I(X;Y) = H(Y)-H(Y|X)
%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% play with p and delta, and observe the plot
%%% test1:set p = [0.01,10], with delta = 12
%%% test2:set delta = [1:100], with p = 0.1 or 0.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% test1:
%%% note: when testing this part, comment out the other part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p=0.01:0.1:15.91
% y=[]
% delta = 6
% for i = 1:160
%     y(i) = mutual_info(@pdf_4pam,p(i),delta)
% end
% plot(p,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% test2:
%%% note: when testing this part, comment out the other part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = 1:0.2:20.8
p=delta-100
y=[]

for i = 1:100
    y(i) = mutual_info(@pdf_4pam,p(i),delta(i))
end
plot(delta,y)
plot(p,y)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculate the final answer
%  H(Y)- H(Y|X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = mutual_info(f,p,delta)
    h1 = entropy_channel(f,p,delta)
    h2 = 0.5*log(2*pi*(delta^2))
    h = h1-h2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculate the entropy of Y
%  Which is the integration of the -pdf(y)log(pdf(y))
%  H(Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = entropy_channel(f,p,delta)
 %f=@func;
 
 
 y = integral(@(x)f(x,p,delta),-30,30)
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculate the pdf of Y
%  And the part inside the integration
%  pdf(y)
%  -pdf(y)log(pdf(y))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = pdf_4pam(x,p,delta)
    %p=0
    %delta=1
    con = 1/(delta*sqrt(2*pi))
    A = con.*exp(-((x-p)./delta).^2./2)
    B = con.*exp(-((x+p)./delta).^2./2)
    C = A + B
    y = -C.*log(C)
end

