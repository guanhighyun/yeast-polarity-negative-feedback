%Written by Tim Elston and modified by Kaiyun Guan
%Model of a bistable system with fast positive feedback and slow negative feedback. 
%Signal affects positive feedback

clear all; close all;

%Set parameter values
k0 = .01;
k1 = 1.5; 
k2 = 1; 
k3 = 2; 

ep = .4;
k4 = ep*0.9;
a2 = ep*.8;

sh = 5;  %high signal
si = 2.5; %intermediate signal
sl = 0.5;  %low signal

% spatial discretization
UT = 4;
xvals = [0.00001:.01:UT];

%compute nulclines with no signal
x_nc = k4/a2*xvals;
u_nc = 1/k3*(((k0+k1*xvals.^2).*(UT-xvals))./xvals - k2);

% Figure 3A
figure(1)
clf 
plot(x_nc,xvals,'-b','linewidth',1.5)
hold on 
plot(u_nc,xvals,'-k','linewidth',1.5)
axis([0 5 0 UT])
xlabel('X')
ylabel('C_{GTP}')
set(gca,'fontsize',25)

%compute nulclines with high signal
s = sh;
k1_s = k1+s;
xh_nc = k4/a2*xvals;
uh_nc = 1/k3*(((k0+k1_s*xvals.^2).*(UT-xvals))./xvals - k2);

figure(1)
plot(uh_nc,xvals,'-r','linewidth',1.5); xlim([-0.05,5])

%compute nulclines with low signal
s = sl;
k1_s = k1+s;
xl_nc = k4/a2*xvals;
ul_nc = 1/k3*(((k0+k1_s*xvals.^2).*(UT-xvals))./xvals - k2);


dt = .01;  %time step
nits = 30000; %number of time steps
tout = [0:dt:nits*dt];  %time points 

%Run the equations to steady state with no signal
%initial conditions. Vary x1 to investigate bistability. 

x1 = 0; 
x2 = (k4/a2)*x1;
for i = 1:nits
   x1_n = x1 + dt*((k0 + k1*x1.^2 ).* (UT-x1) - (k2+k3*x2)*x1);
   x2_n = x2 + dt*(k4*x1 - a2*x2);
   x1 = x1_n;
   x2 = x2_n;
end

%Turn signal to high
s = sh;
k1_s = k1+s;
%output vectors
x1out = zeros(1, nits+1);
x2out = zeros(1, nits+1);

x1out(1) = x1;
x2out(1) = x2;

%iterate equations: Euler method. 
for i = 1:nits
   x1_n = x1 + dt*((k0 + k1_s*x1.^2 ).* (UT-x1) - (k2+k3*x2)*x1);
   x2_n = x2 + dt*(k4*x1 - a2*x2);
   x1 = x1_n;
   x2 = x2_n;
   
   x1out(i+1) = x1;
   x2out(i+1) = x2;
end 

x1h = x1;
x2h = x2;

% Figure 3B
figure(2)
clf
plot(tout, x1out,'-k',tout,x2out,'-b','linewidth',1.5)
axis([0 10 0 5])
%hold on 
%title('Turn on high signal')
ylabel('C_{GTP}, X')
xlabel('Time')
legend('C_{GTP}','X')
set(gca,'fontsize',25)

figure(1)
plot(x2out,x1out,'--r','LineWidth',3)
set(gca,'fontsize',25)

%remove signal
s = 0;
k1_s = k1+s;

x1out(1) = x1h;
x2out(1) = x2h;

for i = 1:nits
   x1_n = x1 + dt*((k0 + k1_s*x1.^2 ).* (UT-x1) - (k2+k3*x2)*x1);
   x2_n = x2 + dt*(k4*x1 - a2*x2);
   x1 = x1_n;
   x2 = x2_n;
   
   x1out(i+1) = x1;
   x2out(i+1) = x2;
end 
% Figure 3C
figure(3)
plot(tout, x1out,'-k',tout,x2out,'-b','linewidth',1.5)
axis([0 3 0 5])
%title('Turn off high signal')
ylabel('C_{GTP}, X')
xlabel('Time')
set(gca,'fontsize',25)

%title('Turn off signal')
ylabel('C_{GTP}, X')
xlabel('Time')

figure(1)
plot(x2out,x1out,'--r','LineWidth',3)

%Drop to low signal.
s = sl;
k1_s = k1+s;

x1 = x1h;
x2 = x2h;

x1out(1) = x1h;
x2out(1) = x2h;


for i = 1:nits
   x1_n = x1 + dt*((k0 + k1_s*x1.^2 ).* (UT-x1) - (k2+k3*x2)*x1);
   x2_n = x2 + dt*(k4*x1 - a2*x2);
   x1 = x1_n;
   x2 = x2_n;
   x1out(i+1) = x1;
   x2out(i+1) = x2;
end 
% Figure 3D
figure(4)
clf
plot(tout, x1out,'-r','linewidth',1.5)
axis([0 3 0 5])
hold on
ylabel('C_{GTP}')
xlabel('Time')
set(gca,'fontsize',25)


%Start with intermediate signal.
s = si;
k1_s = k1+s;

x1 = x1h;
x2 = x2h;
x1out(1) = x1;
x2out(1) = x2;

for i = 1:nits
   x1_n = x1 + dt*((k0 + k1_s*x1.^2 ).* (UT-x1) - (k2+k3*x2)*x1);
   x2_n = x2 + dt*(k4*x1 - a2*x2);
   x1 = x1_n;
   x2 = x2_n;
   
   x1out(i+1) = x1;
   x2out(i+1) = x2;
end 

x1_i = x1;
x2_i = x2;

%Drop from intermediate to low signal.
s = sl;
k1_s = k1+s;

x1 = x1_i;
x2 = x2_i;
x1out(1) = x1;
x2out(1) = x2;

for i = 1:nits
   x1_n = x1 + dt*((k0 + k1_s*x1.^2 ).* (UT-x1) - (k2+k3*x2)*x1);
   x2_n = x2 + dt*(k4*x1 - a2*x2);
   x1 = x1_n;
   x2 = x2_n;
   
   x1out(i+1) = x1;
   x2out(i+1) = x2;
end 
figure(4)
plot(tout, x1out,'color',[0.2 0.8 0.6],'LineWidth',2)
legend('s = 5','s = 2.5')
