function [C] = MatrixMultiplication (A, B)
n = input('Enter no of rows of Matrix A');
m = input('Enter no of columns of Matric A');
A = zeros(n,m)
column = 0;
row = 0;
for column = 1:n
    for row = 1:m
        c = ('Enter value of A');
        A(column,row) = input(c)
    end
end
k = input('Enter no of rows of Matrix B');
l = input('Enter no of columns of Matric B');
B=zeros(k,l);
column = 0;
row = 0;
for row = 1:k
    for column = 1:l
        s = ('Enter value of B');
        B(row,column) = input(s)
    end
end
column=0;
row=0;
for column =1:n
       for row= 1:l
           C(column,row)=A(column,:)*B(:,row);
       end
end
return 