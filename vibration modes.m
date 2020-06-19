clear;
clc;

N=300; % Number of elements

rho=2700; % Density
E=7E10; % Young Modulus
A=0.5*0.5; % Cross section area
I=(0.5*0.5^3)/12; % Second moment of inertia
Le=5; % Length

L=Le/N;

% Stiffness Matrix

k=(E*I/(L^3))*[12, 6*L, -12, 6*L; 6*L, 4*L^2, -6*L, 2*L^2; -12, -6*L, 12, -6*L; 6*L, 2*L^2, -6*L, 4*L^2];
k1=(E*I/(L^3))*[24, 0, -12, 6*L; 0, 8*L^2, -6*L, 2*L^2; -12, -6*L, 24, 0; 6*L, 2*L^2, 0, 8*L^2];

% Mass Matrix

m=(rho*A*L/420)*[156, 22*L, 54, -13*L; 22*L, 4*L^2, 13*L, -3*L^2; 54, 13*L, 156, -22*L; -13*L, -3*L^2, -22*L, 4*L^2];
m1=(rho*A*L/420)*[312, 0, 54, -13*L; 0, 8*L^2, 13*L, -3*L^2; 54, 13*L, 312, 0; -13*L, -3*L^2, 0, 8*L^2];

for i=1:4
    for j=1:4
            K(i,j)=k(i,j);
            K(2*N-2+i,2*N-2+j)=k(i,j);
    end
    
end
for n=1:(N-2)
    for i=1:4
          for j=1:4
          K(i+2*n,j+2*n)=k1(i,j);
          end
    end
end
K
for i=1:4
    for j=1:4
            M(i,j)=m(i,j);
            M(2*N-2+i,2*N-2+j)=m(i,j);
    end
end

for n=1:(N-2)
    for i=1:4
          for j=1:4
          M(i+2*n,j+2*n)=m1(i,j);
          end
    end
end
M
% Boundary conditions (cantilever beam)

    K(1,:)=[]; 
    K(1,:)=[]; % Second row
    K(:,1)=[];
    K(:,1)=[]; % Second column

    M(1,:)=[];
    M(1,:)=[]; % Second row
    M(:,1)=[];
    M(:,1)=[]; % Second column
    
% Eigenvectores and eigenvalues

[V,D]=eig(K,M);

omega=sqrt(D);
f=omega/(2*pi);

% Comparison between finite element and exact solution

f1 = ((1.8751^2)/(2*pi*Le^2))*sqrt((E*I)/(rho*A))
f(1,1)

% Vibration modes (first, second and third modes)

for i=1:(size(V,1)/2)
    X(i,1)=V(2*i-1,1);
    X(i,2)=V(2*i-1,2);
    X(i,3)=V(2*i-1,3);
end

plot(X(:,1),'red','lineWidth',2); hold on;
plot(X(:,2),'blue','lineWidth',2); hold on;
plot(X(:,3),'magenta','lineWidth',2);