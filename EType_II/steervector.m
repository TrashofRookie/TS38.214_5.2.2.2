function [a, x, y, r, r_0] = steervector(phi, mean_phi, theta, lambda, Nt, N1, N2, d)
steer_h = zeros(N1,1);
steer_v = zeros(N2,1);
r = zeros(Nt,1);
r_0 = zeros(Nt,1);
x = zeros(Nt,1);
y = zeros(Nt,1);
for i = 1 : N2
    for j = 1 : N1
        x((i-1)*N2 + j) = (j-1)*d;
        y((i-1)*N2 + j) = (i-1)*d;
    end
end
for i = 1 : N1
    steer_h(i) = (1/sqrt(N1))*exp(-1i*2*pi*d/lambda*sin(theta)*(i-1)*sin(mean_phi-phi));
end
for i = 1 : N2
    steer_v(i) = (1/sqrt(N2))*exp(-1i*2*pi*d/lambda*(i-1)*cos(mean_phi-phi));
end
a = kron(steer_h,steer_v);
r(1) = 0;
r_0(1) = 0;
for i = 2 : Nt
    r(i) = sqrt(x(i)^2+y(i)^2);
    r_0(i) = atan(y(i)/x(i));
end
end
