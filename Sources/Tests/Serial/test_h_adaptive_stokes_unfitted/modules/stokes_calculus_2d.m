clear
close all
clc

syms q

N = 2;
x = sym('x%d',[N 1]);

u = [x(2) ; x(1) ] ;% [ x(1)^q*x(2)^q+x(1)^q ; x(1)^q*x(2)^q+x(2)^q];
p = 1  ;% x(1)^(q-1)*x(2)^(q-1);


grad_u = sym(zeros(N,N));
for i=1:N
    for j=1:N
        grad_u(i,j) = diff(u(j),x(i));
    end
end

div_u = sym(zeros(1,1));
for i=1:N
    div_u(1) = div_u(1) + diff(u(i),x(i));
end

div_grad_u = sym(zeros(N,1));
for i=1:N
    for j=1:N
        div_grad_u(i) = div_grad_u(i) + diff(grad_u(j,i),x(i));
    end
end

grad_p = sym(zeros(N,1));
for i=1:N
    grad_p(i) = diff(p,x(i));
end


for i = 1:N
    fprintf('u(%i)= ',i);
    disp(u(i))
end

for i = 1:N
    for j = 1:N
        fprintf('grad_u(%i,%i)= ',i,j);
        disp(grad_u(i,j))
    end
end

fprintf('div_u = ')
disp(div_u)

for i = 1:N
    fprintf('div_grad_u(%i)= ',i);
    disp(div_grad_u(i))
end

fprintf('p = ')
disp(p)

for i = 1:N
    fprintf('grad_p(%i)= ',i);
    disp(grad_p(i))
end


