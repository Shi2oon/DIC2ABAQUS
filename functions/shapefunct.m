function [H,dHex] = shapefunct(eta, xi, ~)

%Shape functions
H(1) = 1/4*(1-eta)*(1-xi);
H(2) = 1/4*(1+eta)*(1-xi);
H(3) = 1/4*(1+eta)*(1+xi);
H(4) = 1/4*(1-eta)*(1+xi);
%Differentiation of the shape functions
%[dNi/d_eta ... ; dNi/d_xi ...]
dHex(1,1) = xi/4 - 1/4;
dHex(2,1) = eta/4 - 1/4;
dHex(1,2) = 1/4 - xi/4;
dHex(2,2) = - eta/4 - 1/4;
dHex(1,3) = xi/4 + 1/4;
dHex(2,3) = eta/4 + 1/4;
dHex(1,4) = - xi/4 - 1/4;
dHex(2,4) = 1/4 - eta/4;

