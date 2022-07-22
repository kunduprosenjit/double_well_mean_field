function [x_eff,beta_eff] = betaspace(A,x)

deg1 = sum(A);
x_n = A*x;
if sum(sum(A)) == 0 % if the network has no connectivity
    % beta xnn are not well defined when there is no network
    beta_eff = 0;
    x_eff = 0;
else % the network has its connectivity
    beta_eff = full(sum(sum(A*A))/sum(sum(A))); 
    x_eff = sum(x_n)/sum(deg1); 
end
