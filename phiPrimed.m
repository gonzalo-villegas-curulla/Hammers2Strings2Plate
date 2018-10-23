function y = phiPrimed (K, alpha, eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functionm phiPrimed() returns the derivative of phi()
% for Newton-Raphson solver purposes. The variables-in
% remain the same.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y = K*(eta/2 + abs(eta)/2)^alpha;
% NOTE: Truncation error potentially could
% explain the wiggling in dHdt history, below a threshold of 1e-11.

end
