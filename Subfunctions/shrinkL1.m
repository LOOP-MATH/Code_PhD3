 % This code solves the following problem:
%    min 0.5*||x-c||^2 + mu*||x||_1
% input:
% c      - the given point
% mu     - the l1 penalty parameter
% output
% d      - soft thresholding
% Download from Tingkei Pong's homepage.
function d = shrinkL1(c,mu)

d = sign(c).*max(0,abs(c) - mu);
end