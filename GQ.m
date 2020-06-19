function R = GQ() % Gaussian Quadrature

str = input('Enter the equation f(x): ','s');
f = inline(str);
n = input ('Enter the number of points: ');
A = zeros(1,n);

for i = 1:n
    A(i) = 1/2*((1-(2*i)^(-2))^(-0.5));
end

B = zeros(n);
for i=1:n-1
    B(i,i+1) = A(i); 
    B(i+1,i) = A(i);
end

x = eig(B); 
[C, D] = eig(B);

w = zeros(n,1);
for i = 1:n
    w(i,1) = 2*((C(1,i)^2));
end

sum = 0;
for i = 1:n
    sum = sum + (w(i)*f(x(i)));
end
fprintf('The value of the integral of f(x) within limits -1 to 1 :');
disp(sum)
