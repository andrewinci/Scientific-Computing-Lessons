clear all
close all
clc

fprintf('Gauss method \n')
A=[1,2,3;
    4,5,6;
    7,8,0];
b=[1;1;1];
xe=A\b % exact solution computed by MATLAB
x=Gauss(A,b)

fprintf('Thomas method \n')
A=[3,-1,0;
    -1,3,-1;
    0,-1,3];
d=[1;1;1];
xe=Gauss(A,d)
a=[-1,-1,-1];
b=[3,3,3];
c=[-1,-1,-1];
x=Thomas(a,b,c,d)

fprintf('Conjugate Gradient method \n')
xe=Gauss(A,d)
x=CG(A,d)

fprintf('Matrix-free conjugate Gradient method \n')
xe=Gauss(A,d)
x=CGop(A,d)