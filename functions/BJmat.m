function [B,J] = BJmat(dHex, el, Ue)

%Jacobian matrix - CHECKED
 J = dHex*[Ue(1:2:end);Ue(2:2:end)]';
%  Jequivalent = [sum(dHex(1,:).*Ue(1:2:end)) sum(dHex(1,:).*Ue(2:2:end));sum(dHex(2,:).*Ue(1:2:end)) sum(dHex(2,:).*Ue(2:2:end))];

%shape function for x and y
% shape function derivative in the real coordinate system
%deltaNi/deltax and deltaNi?deltay
 dHxy = J\dHex; %or inv(J)*dHex
%Assemble B matrix
 B = zeros(3,size(el.n,2)*2);
 B(1,1:2:end) = dHxy(1,:);
 B(2,2:2:end) = dHxy(2,:);
 B(3,1:2:end) = dHxy(2,:);
 B(3,2:2:end) = dHxy(1,:);