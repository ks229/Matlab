function S = GE() % Gaussian Elimination

a = 'No. of Rows and Columns of Matrix A : ';
n = input(a);


for i = 1:n
    for j = 1:n   
        fprintf('A%d%d',i,j);
        A(i,j) = input('\n ');
    end
end
disp ('A Matrix: ')
disp(A)


for k = 1:n
    fprintf('B%d1',k);
    B(k,1) = input('\n ');
end
disp ('B Matrix: ')
disp(B)


for j = 1:n
    for i = 1:n
        if(i>j)
            
                A(i,j) = A(i,j)-((A(i,i)/A(i,i))*A(i,j));
            
        end
    end
end
disp('A Matrix after Forward Elimination')
disp(A)


X = zeros(n,1);
for i = n:-1:1
    s=0;
    for j = i+1:n
            s = s+(A(i,j)*Y(j,1));
    end
    Y(i,1) = (B(i,1)-s)/A(i,i);
    X = Y;
end
disp('The Solution X Matrix: ')
disp(X)