function [obj, grad,FF] = obj_EPG13_exc_wrapper(params,exc,ESP,T1,T2,c,B1,target,frequencies,klim)
% wrapping obj_EPG13 to provide arbitrary excitation pattern and exclude
% from optimization.
% Here expect params to be without excitation. Everything else with
% 'excitation'.

nt = size(target,1);
[Ns, Nch] = size(B1);
np = length(params)/(2*Nch);

x = reshape(params(1:np*Nch),np,Nch).'; %real parts
y = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts

% mz add excitation back
x = [zeros(Nch, 1) x];
y = [zeros(Nch, 1) y];

params0 = x +1i*y;
params1 = zeros(Nch,sum(frequencies));
sumfreq = cumsum(frequencies);

for j = (length(frequencies)):-1:2
    params1(:,sumfreq(j-1)+1:sumfreq(j)) = repmat(params0(:,j),[1 frequencies(j)]);
end

params1(:,1:frequencies(1)) = repmat(params0(:,1),[1 frequencies(1)]);
params1 = params1.';
params2 = [real(params1(:).') , imag(params1(:).')];
params = params2;


np = length(params)/(2*Nch);

x = reshape(params(1:np*Nch),np,Nch).'; %real parts
y = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts
th = zeros(Ns,np,Nch);
for ch = 1:Nch
    th(:,:,ch) = B1(:,ch)*(x(ch,:)+1i*y(ch,:)); % effective theta
end
th(:,1,:) = exc;

[obj, grad,FF] = obj_EPG13_fromeff(th,ESP,T1,T2,c,B1,target,frequencies,klim);
grad = reshape(grad, np, Nch, 2);
grad = grad(2:end, :, :);
grad = grad(:);
