% BOUNDARY VALUE PROBLEM - SOLUTION WITH LINEAR FINITE ELEMENTS

clc; clear variables; close all;

% INITIALIZATIONS

L = 1; % Domain
N = 3; % No. of elements - Adjust assembly of global SM for more elements
n = 2*N+1; % No. of nodes
he = L/(N); % Element length
K = zeros(n,n); % Global Stiffness matrix initialization
syms w1(x) w2(x) w3(x) f(x)  h;
w1(x) = (1 - x/h)*(1-2*x/h);   % Shape functions at local coordinate system - Hermite polynomials
w2(x) = (4*x/h)*(1 - x/h);
w3(x)=(-x/h)*(1-2*x/h);
w = [w1(x) w2(x) w3(x)];
f(x) = -(1-2*(2*x+cos(2*x)));    % fmass

% ELEMENT STIFFNESS MATRIX
 
for i=1:3
    for j=1:3
        ke(i,j) = int(diff(w(i),x)*diff(w(j),x)+w(i)*diff(w(j),x) -4* w(i)*w(j),x,0,h);
    end
end
ke = vpa(subs(ke,h,he));

% GLOBAL STIFFNESS MATRIX ASSEMBLE

K=QuadraticBarAssemble(K,ke,1,2,3);
K=QuadraticBarAssemble(K,ke,3,4,5);
K=QuadraticBarAssemble(K,ke,5,6,7);
%K=QuadraticBarAssemble(K,ke,7,8,9);
%K=QuadraticBarAssemble(K,ke,9,10,11);
%K=QuadraticBarAssemble(K,ke,11,12,13);

% ELEMENT FORCE VECTOR

syms p1(x) p2(x) xb xa
p1(x) = ((xb-x)/h)*((xb+xa-2*x)/h);
p2(x) = (4*(x-xa)/h)*((xb-x)/h);   
p3(x)=(-(x-xa)/h)*((xb+xa-2*x)/h);
p = [p1(x) p2(x) p3(x)];
dp1(x)= diff(p1,x);
dp2(x)=diff(p2,x);
dp3(x)=diff(p3,x);
dp=[dp1 dp2 dp3];

f1 = int(f(x)*p1(x),x,xa,xb); % Integration
f2 = int(f(x)*p2(x),x,xa,xb);
f3 = int(f(x)*p3(x),x,xa,xb);
f=[f1 f2 f3];
xaa=0;
xbb=he;
f_el=zeros(N,3);
for i=1:N
    for j=1:3
        f_el(i,j)= vpa(subs(f(j),{xa,xb,h},{xaa,xbb,he}));
    end
    xaa=xaa+he;
    xbb=xbb+he;
end

% FORCES AT EACH NODE

fe=zeros(n,1);
fe(1,1)=f_el(1,1);
fe(n,1)=f_el(N,3);    
j=[1 2 3];
m= 1:0.5:2*N;
for i= 2:2:n-1
        fe(i,1)=f_el(i/2,j(2));
        for l = 3:2:n-2
               fe(l,1) = f_el(l-1-m(l-2),j(3)) + f_el(l-m(l-2),j(1));
        end
end

% IMPOSE BOUNDARY CONDITIONS - REDUCE STIFFNESS MATRIX 

K_final=zeros(n-1,n-1);
for i=1:n-1
    for j=1:n-1
        K_final(i,j)=K(i+1,j+1);    
    end
end

% REDUCE FORCE VECTOR

f_final=zeros(n-1,1);
for i=1:n-1
     f_final(i)=fe(i+1);
end

% SOLUTION OF THE LINEAR SYSTEM 

dirichlet=sin(2)-1; % Nonhomogeneous boundary condition
f1_final=zeros(n-1,1);
for i = 1:n-1
    f1_final(i) = f_final(i)-K_final(i,n-1)*dirichlet;  % Changes in right hand side vector
end
f_solve=zeros(n-2,1);
for i = 1:n-2
    f_solve(i)=f1_final(i);
end
K_solve=zeros(n-2,n-2);
for i=1:n-2
    for j=1:n-2
            K_solve(i,j)= K_final(i,j); % Changes in stiffness matrix
    end
end

% SOLVE FOR DISPLACEMENTS AT EACH NODE

U=K_solve\f_solve;
U_final = zeros(n,1);
U_final(1)=0;
U_final(n)= dirichlet; %boundary conditions
for i = 2:n-1
    U_final(i) = U(i-1);
end

% POSTPROCESSING - PLOT SOLUTION AT WHOLE DOMAIN - SOLUTION IS LIKE PIECEWISE FUNCTION

syms U_solution
xaa=0;
xbb=he;
for i = 1:N
     j=1:2:n+1;
     U_solution(i) = U_final(j(i))*subs(p1(x),{xa, xb, h}, {xaa,xbb,he}) +  U_final(j(i)+1)*subs(p2(x), {xa, xb, h}, {xaa,xbb,he})+U_final(j(i)+2)*subs(p3(x),{xa, xb, h}, {xaa,xbb,he});  %solution
     xaa=xaa+he;
     xbb=xbb+he;
end
figure
for i=1:N
     z=linspace((i-1)*he,i*he);
     sol=subs(U_solution(i),x,z);  %solution plot
     plot (z,sol,'r', 'LineWidth',2)
     hold on 
end 
xlabel('x')
ylabel('u')
legend('1D Quadratic Elements')

% PLOT EXACT SOLUTION AT WHOLE DOMAIN

l=linspace(0,1,100);
ExactSolution=@(l) sin(2*l)-l;
plot(l,ExactSolution(l),'b--','DisplayName','UExact', 'LineWidth',2)
title(sprintf('U(x) for N= %.0f Quadratic Finite Elements', N))
hold off
 
% PLOT DERIVATE OF U AT WHOLE DOMAIN
 
syms U_derivative
xaa=0;
xbb=he;
for i = 1:N
     j=1:2:n+1;
     U_derivative(i) = U_final(j(i))*subs(dp1(x),{xa, xb, h}, {xaa,xbb,he}) +  U_final(j(i)+1)*subs(dp2(x), {xa, xb, h}, {xaa,xbb,he})+  U_final(j(i)+2)*subs(dp3(x), {xa, xb, h}, {xaa,xbb,he});  %solution for the derivative
     xaa=xaa+he;
     xbb=xbb+he;
end
figure
for i=1:N
     z=linspace((i-1)*he,i*he);
     der=subs(U_derivative(i),x,z);  
     plot (z,der,'LineWidth',2)
     hold on 
end
xlabel('x')
ylabel('du')
legend('du(x)')
title(sprintf('dU/dx for N= %.0f Quadratic Finite Elements', N))
hold off

% FUNCTIONS NEEDED

function y = QuadraticBarAssemble(K,k,i,j,m)
K(i,i) = K(i,i) + k(1,1);
K(i,j) = K(i,j) + k(1,2);
K(i,m) = K(i,m) + k(1,3);
K(j,i) = K(j,i) + k(2,1);
K(j,j) = K(j,j) + k(2,2);
K(j,m) = K(j,m) + k(2,3);
K(m,i) = K(m,i) + k(3,1);
K(m,j) = K(m,j) + k(3,2);
K(m,m) = K(m,m) + k(3,3);
y = K;
end
