f=@(x)sqrt(x);
a=0; 
b=1; 
f_exact= integral(f,a,b);
Error_Trap=zeros(11,1);
Error_Simp=zeros(11,1);
for m = [10,20,40,80,160,320,640,1280,2560,5120,10240]  
    if rem(m,2)==1 % to make sure the eligibility of Simpson's Rule
        fprintf('\n Enter valid m!!!'); 
        m=input('\n Enter m as even number ');
    end
    h=(b-a)/m;
    sum=0; % here we partially calculate the Trapezoidal's rule
    for i=1:1:m-1
    x(i)=a+i*h;
    y(i)=f(x(i));
    sum=sum+y(i); 
    end
    so=0;se=0; % here we partially calculate the Simpson's rule
    for s=1:1:m-1
        if rem(s,2)==1
            so=so+y(s);%sum of odd terms
        else
            se=se+y(s); %sum of even terms
        end
    end
    % Formula for Trapezoidal's Rule:  (h/2)*[(y0+yn)+2*(y2+y3+..+yn-1)]
    Trap=h/2*(f(a)+f(b)+2*sum);
    Error_Trap(log2(m/10)+1)= abs(f_exact-Trap);
    % Formula for Simpson's Rule:  (h/3)*[(y0+yn)+2*(y3+y5+..odd term)+4*(y2+y4+y6+...even terms)]
    Simp=h/3*(f(a)+f(b)+4*so+2*se);
    Error_Simp(log2(m/10)+1)=abs(f_exact-Simp);
end
m = [10,20,40,80,160,320,640,1280,2560,5120,10240];
x = log(m);
y1 = log(Error_Trap);
y2 = log(Error_Simp);
plot(x,y1,x,y2)
title('log-log diagram of quadrature errors vs m')
xlabel('log_m') 
ylabel('log_error')
legend({'Composite Trapezoidal Rule','Composite Simpson Rule'},'Location','southwest')
ax = gca;
ax.FontSize = 13;

% Here we use least-square fitting to find k and D in two different diagram
% Trapezoidal's Rule
x_1 = x;
y_1 = y1;
coefs1 = LeastSquare(x_1,y_1)

% Simpson's Rule
x_2 = x;
y_2 = y2;
coefs2 = LeastSquare(x_2,y_2)

function coefs = LeastSquare(x,y)
 A = [x' 0*x'+1];
 b = y;
 [Q,R] = qr(A);
 coefs = R\(Q'*b);

 end