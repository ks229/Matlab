clear
%properties
den = 3000;
wd = 0.1; ht = 0.2;  len = 0.5;
area = wd * ht;
E = 10^10;
I = (len*wd^3)/12;
n=input('Enter no of nodes :');
l = len/(n-1);
M1 = (den*l*area/420)*[156 22*l 54 -13*l;22*l 4*l*l 13*l -3*l*l;54 13*l 156 -22*l;-13*l -3*l*l -22*l 4*l*l];
K1 = (E*I/l^3)*[12 6*l -12 6*l;6*l 4*l*l -6*l 2*l*l;-12 -6*l 12 -6*l;6*l 2*l*l -6*l 4*l*l];
K = zeros(2*n,2*n);
M = zeros(2*n,2*n);
%Matrix stiffness and mass
for i=1:2:2*(n-1)
    for j=1:4
        for l=1:4
            K(i+j-1,i+l-1) = K(i+j-1,i+l-1) + K1(j,l);
            M(i+j-1,i+l-1) = M(i+j-1,i+l-1) + M1(j,l);
        end
    end
end
%boundary conditions
    K(1,:)=[];
    K(1,:)=[]; % Second row
    K(:,1)=[];
    K(:,1)=[]; % Second column
    M(1,:)=[];
    M(1,:)=[]; % Second row
    M(:,1)=[];
    M(:,1)=[]; % Second column 
%initialise values
a = 10;
w = 100;
alpha = 0.25;
gamma = 0.5;
dt = 1;     %timestep
tf = 20;     %time end
fprintf('Freq is :');
Freq = w / (2*pi);
disp(Freq);
%integration constants
a0 = 1/(alpha*dt^2); 
a1 = gamma/(alpha*dt);
a2 = 1/(alpha*dt);  
a3 = 1/(2*alpha)-1;
a4 = (gamma/alpha)-1; 
a5 = dt/2*(gamma/alpha-2);
a6 = dt*(1-alpha);  
a7 = (gamma/alpha-2);
t = 0;
i = 1;
U = zeros(2*(n-1),tf/dt);
V = zeros(2*(n-1),tf/dt);
A = zeros(2*(n-1),tf/dt);
V(:,i+1) = V(:,i) + [(1-gamma)*A(:,i) + gamma*A(:,i+1)]*dt;
U(:,i+1) = U(:,i) + dt*V(:,i) + [(0.5-alpha)*A(:,i) + alpha*A(:,i+1)]*dt^2;
F = zeros(2*(n-1),tf/dt);
Keff= K + a0*M;
while t<tf
    F(2*(n-1),i) = a*sin(w*t);
    Feff = F(:,i) + M*(a0*U(:,i) + a2*V(:,i) + 2*A(:,i));
    U(:,i+1) =inv(Keff)*Feff ;
    A(:,i+1) = a0*(U(:,i+1) - U(:,i)) - a2*V(:,i) - a3*A(:,i);
    V(:,i+1) = V(:,i) + a6*A(:,i) + a7*A(:,i+1);
    t = t + dt;
    i = i + 1;
end
plot(U(1,:),'red','lineWidth',2); hold on;
plot(U(2,:),'blue','lineWidth',2); hold on;
plot(U(3,:),'magenta','lineWidth',2);

 


