% BOUNDARY VALUE PROBLEM - SOLUTION WITH LINEAR FINITE ELEMENTS

clc; clear variables; close all;

% INITIALIZATIONS

L = 1; % Domain
N =8; % No. of elements - - Adjust assembly of global SM for more elements
n = N+1; % No. of nodes
he = L/(N); % Element length
K = zeros(n,n); % Global Stiffness matrix initialization
syms w1(x) w2(x) f(x)  h;
w1(x) = 1 - x/h;   % Shape functions at local coordinate system - Lagrange Polynomials
w2(x) = x/h;
w = [w1(x) w2(x)];
f(x) = -(1-2*(2*x+cos(2*x)));    %%% fmass

% ELEMENT STIFFNESS MATRIX

for i=1:2
    for j=1:2
        ke(i,j) = int(diff(w(i),x)*diff(w(j),x)+w(i)*diff(w(j),x)-4*w(i)*w(j),x,0,h);
    end
end
ke = vpa(subs(ke,h,he));

% GLOBAL STIFFNESS MATRIX ASSEMBLY

K=LinearBarAssemble(K,ke,1,2);
K=LinearBarAssemble(K,ke,2,3);
K=LinearBarAssemble(K,ke,3,4);
K=LinearBarAssemble(K,ke,4,5);
K=LinearBarAssemble(K,ke,5,6);
K=LinearBarAssemble(K,ke,6,7);
K=LinearBarAssemble(K,ke,7,8);
K=LinearBarAssemble(K,ke,8,9);

% CALCULATE FORCE VECTOR OVER ELEMENT

syms p1(x) p2(x) xb xa
p1(x) = (xb-x)/h; % shape functions at global coordinate system
p2(x) = (x-xa)/h;
p = [p1(x) p2(x)];
dp1(x)= diff(p1,x);
dp2(x)=diff(p2,x);
dp=[dp1 dp2];

f1 = int(f(x)*p1(x),x,xa,xb); % Integration
f2 = int(f(x)*p2(x),x,xa,xb);
f=[f1 f2];

xaa=0;
xbb=he;
f_el=zeros(N,2);
for i=1:N
    for j=1:2
        f_el(i,j)= vpa(subs(f(j),{xa,xb,h},{xaa,xbb,he}));
    end
    xaa=xaa+he;
    xbb=xbb+he;
end


% FORCES AT EACH NODE

fe=zeros(n,1);
for i=2:n-1
    fe(i,1)=f_el(i-1,2)+f_el(i,1);
end
fe(1,1)=f_el(1,1);
fe(n,1)=f_el(N,2);

% IMPOSE BOUNDARY CONDITIONS - REDUCE THE STIFFNESS MATRIX

K_final=zeros(n-1,n-1);
for i=1:n-1
    for j=1:n-1
        K_final(i,j)=K(i+1,j+1);    
    end
end

% REDUCED FORCE VECTOR

f_final=zeros(n-1,1);
for i=1:n-1
     f_final(i)=fe(i+1);
end

% SOLUTION OF THE LINEAR SYSTEM

dirichlet=sin(2)-1; % Non-Homogeneous Dirchlett boundary condition
f1_final=zeros(n-1,1);
for i = 1:n-1
    f1_final(i) = f_final(i)-K_final(i,n-1)*dirichlet;  % changes in right hand side vector
end
f_solve=zeros(n-2,1);
for i = 1:n-2
    f_solve(i)=f1_final(i);
end
K_solve=zeros(n-2,n-2);
for i=1:n-2
    for j=1:n-2
            K_solve(i,j)= K_final(i,j); % changes in stiffness matrix
    end
end
U=K_solve\f_solve; % solve the system for displacements at each node
U_final = zeros(n,1);
U_final(1)=0;
U_final(n)= dirichlet; %boundary conditions
for i = 2:n-1
    U_final(i) = U(i-1);
end

% PLOT SOLUTION AT WHOLE DOMAIN - SOLUTION IS LIKE PIECEWISE FUNCTION

syms U_solution
xaa=0;
xbb=he;
for i = 1:N
     U_solution(i) = U_final(i)*subs(p1(x),{xa, xb, h}, {xaa,xbb,he}) +  U_final(i+1)*subs(p2(x), {xa, xb, h}, {xaa,xbb,he});  %solution
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
legend('1D Linear Elements')

% PLOT EXACT SOLUTION AT WHOLE DOMAIN

l=linspace(0,1,100);
ExactSolution=@(l) sin(2*l)-l;
plot(l,ExactSolution(l),'b--','DisplayName','UExact', 'LineWidth',2)
title(sprintf('U(x) for N= %.0f Linear Finite Elements', N))
hold off
 
% PLOT DERIVATE OF U AT WHOLE DOMAIN
 
syms U_derivative
xaa=0;
xbb=he;
for i = 1:N
     U_derivative(i) = U_final(i)*subs(dp1(x),{xa, xb, h}, {xaa,xbb,he}) +  U_final(i+1)*subs(dp2(x), {xa, xb, h}, {xaa,xbb,he});  %solution for the derivative
     xaa=xaa+he;
     xbb=xbb+he;
end
figure
for i=1:N
     z=linspace((i-1)*he,i*he);
     der=subs(U_derivative(i),x,z);  
     plot (z,der, 'LineWidth',2)
     hold on 
 end
 xlabel('x')
 ylabel('du')
 legend('du(x)')
 title(sprintf('dU/dx for N= %.0f Linear Finite Elements', N))
 hold off

 % FUNCTIONS NEEDED

function y = LinearBarAssemble(K,k,i,j)
K(i,i) = K(i,i) + k(1,1) ;
K(i,j) = K(i,j) + k(1,2) ;
K(j,i) = K(j,i) + k(2,1) ;
K(j,j) = K(j,j) + k(2,2) ;
y = K;
end
 
