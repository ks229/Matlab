clear
%Material properties
den = 3000;
wd = 0.1;
ht = 0.2;
len = 0.5;
area = wd * ht;
E = 10^10;
I = (len*wd^3)/12;
n=input('Enter no of nodes :');
l = len/(n-1);
%Stiffness matrix and Mass matrix
M1 = (den*l*area/420)*[156 22*l 54 -13*l;22*l 4*l*l 13*l -3*l*l;54 13*l 156 -22*l;-13*l -3*l*l -22*l 4*l*l];
K1 = (E*I/l^3)*[12 6*l -12 6*l;6*l 4*l*l -6*l 2*l*l;-12 -6*l 12 -6*l;6*l 2*l*l -6*l 4*l*l];
K = zeros(2*n,2*n);
M = zeros(2*n,2*n);
for i=1:2:2*(n-1)
    for j=1:4
        for l=1:4
            K(i+j-1,i+l-1) = K(i+j-1,i+l-1) + K1(j,l); %Stiffness matrix
            M(i+j-1,i+l-1) = M(i+j-1,i+l-1) + M1(j,l); %Mass matrix
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
    %Eigen values
   [A,D] = eig(K,M);
   w = sqrt(D);
   f = w/(2*pi);
   fprintf('Freq =');
   disp(f(1,1));    %frequency
   if (size(A,1)/2 ~= 1)
        for i=1:(size(A,1)/2)
        X(i,1)=A(2*i-1,1);
        X(i,2)=A(2*i-1,2);
        X(i,3)=A(2*i-1,3);
        end
   plot(X(:,1),'red','lineWidth',2); hold on;
   plot(X(:,2),'blue','lineWidth',2); hold on;
   plot(X(:,3),'magenta','lineWidth',2);
   else
       for i=1:(size(A,1)/2);
        X(i,1)=A(2*i-1,1);
        X(i,2)=A(2*i-1,2);
       end
   plot(X(:,1),'red','lineWidth',2); hold on;
   plot(X(:,2),'blue','lineWidth',2);
   end
  