%% Test of translation from triangular to cox 
close all

n=5;
maxInt=5;

startProb=ph_initial(n); % Starting probability
Q=ph_generatorTri(n+1,maxInt); % Upper triangular matrix

Q1=Q(1:end-1,1:end-1);
Q0=Q(1:end-1,end);

[Cox1,Cox0,Cox]=ttc_main(startProb,Q); % Canonical matrix
% Cox1=Cox(1:end-1,1:end-1);
% Cox0=Cox(1:end-1,end);

startProbCox=zeros(1,n); % Start vector for canonical
startProbCox(1)=1;

t=linspace(0,10,100);

% Plots pdf of triangular (blue) and Canonical (red)
figure(1)
clf
hold on
plot(t,ph_pdf(t,Q1,Q0,startProb),'b');
plot(t,ph_pdf(t,Cox1,Cox0,startProbCox),'rx');

sum=sum(ph_pdf(t,Q1,Q0,startProb)-ph_pdf(t,Cox1,Cox0,startProbCox));

% Plots cdf of triangular (blue) and Canonical (red)
figure(2)
clf
hold on
plot(t,ph_cdf(t,Q1,startProb),'b');
plot(t,ph_cdf(t,Cox1,startProbCox),'r');
%% Test of scaling Tri and Cox. Also comparing with optimal scaling
n=5; % all non-absorbing
maxint=2;
scale=2;

t=linspace(0,10);

%triangular generator
startProb=ph_initial(n);
[Lambda, Theta]=ph_generator_tri(n,maxint);
Q=[Lambda Theta;zeros(1,length(Lambda)+1)];

%Cox generator
startCox=[1 zeros(1,n-1)];
[LambdaCox,ThetaCox,Cox]=ttc_main(startProb,Lambda,Theta);

%triangular generator with scaled parameter
Q_scaled=ph_scale_absintensity(Q,scale);
Lambda_scaled=Q_scaled(1:end-1,1:end-1);
Theta_scaled=Q_scaled(1:end-1,end);

%Cox generator with scaled parameter
Cox_scaled=ph_scale_absintensity(Cox,scale);
[LambdaCox_scaled,ThetaCox_scaled,Cox_scaled]=ttc_main(startCox,Cox_scaled);
% LambdaCox_scaled=Cox_scaled(1:end-1,1:end-1);
% ThetaCox_scaled=Cox_scaled(1:end-1,end);

% Computes the residuals between the two distributions
residuals1=ph_pdf(t,Lambda_scaled,Theta_scaled,startProb)-...
    ph_pdf(t,LambdaCox_scaled,ThetaCox_scaled,startCox);

% Finds C2 which is the best approxiamtion for the scale parameter of the
% coxian distribution
[C2,val]=fmincon(@(x) normRes(t,Lambda_scaled,Theta_scaled,...
    startProb,Cox,startCox,x),1,-1,0,[],[]);

Cox_scaled_Opt=ph_scale_absintensity(Cox,C2);
LambdaCox_scaled_Opt=Cox_scaled_Opt(1:end-1,1:end-1);
ThetaCox_scaled_Opt=Cox_scaled_Opt(1:end-1,end);

% Residuals with the optimal scale parameter
residuals2=ph_pdf(t,Lambda_scaled,Theta_scaled,startProb)-...
    ph_pdf(t,LambdaCox_scaled_Opt,ThetaCox_scaled_Opt,startCox);

% Optimization using Kullback Liebler
t1=linspace(0,100,1000);
[C3,val1]=fmincon(@(x) kullbackLiebler(t1,Q_scaled,startProb,Cox,startCox,x),...
    1,-1,0,[],[]);

Cox_scaled_Kullback=ph_scale_absintensity(Cox,C3);
LambdaCox_scaled_Kullback=Cox_scaled_Kullback(1:end-1,1:end-1);
ThetaCox_scaled_Kullback=Cox_scaled_Kullback(1:end-1,end);

residuals3=ph_pdf(t,Lambda_scaled,Theta_scaled,startProb)-...
    ph_pdf(t,LambdaCox_scaled_Kullback,ThetaCox_scaled_Kullback,startCox);
%%
figure(1)
clf
title('pdf Tri and Cox not scaled')
hold on
plot(t,ph_pdf(t,Lambda,Theta,startProb),'b');
plot(t,ph_pdf(t,LambdaCox,ThetaCox,startCox),'r');

figure(2) % Tri and Cox scaled
clf
title('pdf Tri and Cox scaled')
hold on
plot(t,ph_pdf(t,Lambda_scaled,Theta_scaled,startProb),'b')
plot(t,ph_pdf(t,LambdaCox_scaled,ThetaCox_scaled,startCox),'r')



figure(3) % Tri and Cox scaled optimal
clf
title('pdf Tri and Cox scaled optimal')
hold on
plot(t,ph_pdf(t,Lambda_scaled,Theta_scaled,startProb),'b')
plot(t,ph_pdf(t,LambdaCox_scaled_Opt,ThetaCox_scaled_Opt,startCox),'r')

figure(4)
clf
title('pdf Tri and Cox scaled with Kullback')
hold on
plot(t,ph_pdf(t,Lambda_scaled,Theta_scaled,startProb),'b')
plot(t,ph_pdf(t,LambdaCox_scaled_Kullback,ThetaCox_scaled_Kullback,startCox),'r')


figure(5)
clf
title('residuals')
hold on
plot(t,residuals1);
plot(t,residuals2,'k');
plot(t,residuals3,'g');
legend('scaled','scaled optimal','Kullback Liebler')


%% Test of iteration to find relationship between scaled tri and scaled cox
close all

n=10;
maxint=10;
scale=linspace(2,10,100);

t=linspace(0,2);

poly=1;
P=zeros(poly+1,length(scale));

for i=1:length(scale)
%triangular generator
startProb=ph_initial(n-1);
Q=ph_generator_tri(n,maxint);
Lambda=Q(1:end-1,1:end-1);
Theta=Q(1:end-1,end);

%Cox generator
startCox=[1 zeros(1,n-2)];
Cox=ttc_main(Q,startProb);
LambdaCox=Cox(1:end-1,1:end-1);
ThetaCox=Cox(1:end-1,end);

%triangular generator with scaled parameter
Q_scaled=ph_scale_absintensity(Q,scale(i));
Lambda_scaled=Q_scaled(1:end-1,1:end-1);
Theta_scaled=Q_scaled(1:end-1,end);

%Cox generator with scaled parameter
Cox_scaled=ph_scale_absintensity(Cox,scale(i));
LambdaCox_scaled=Cox_scaled(1:end-1,1:end-1);
ThetaCox_scaled=Cox_scaled(1:end-1,end);

residuals=ph_pdf(t,Lambda_scaled,Theta_scaled,startProb)-...
    ph_pdf(t,LambdaCox_scaled,ThetaCox_scaled,startCox);

figure(1)
hold on
plot(ph_pdf(t,Lambda_scaled,Theta_scaled,startProb),...
    ph_pdf(t,LambdaCox_scaled,ThetaCox_scaled,startCox))

P(:,i)=polyfit(ph_pdf(t,Lambda_scaled,Theta_scaled,startProb),...
    ph_pdf(t,LambdaCox_scaled,ThetaCox_scaled,startCox),poly);

%pause(0.5)

% figure(3) % Tri and Cox scaled
% hold on
% plot(t,ph_pdf(t,Lambda_scaled,Theta_scaled,startProb),'b')
% 
% figure(4)
% hold on
% plot(t,ph_pdf(t,LambdaCox_scaled,ThetaCox_scaled,startCox),'r')

end

figure(2)
hold on
plot(scale,P(1,:))








