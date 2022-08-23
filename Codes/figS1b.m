clear;
close all;
clc;
global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu Nm;
PoolNumSpecies=100;
PoolNumPlasmid=10;
NumPlasmid=PoolNumPlasmid;
timespan=0:1000;
eta_star=0.01:0.002:0.04;
etas=2*eta_star;
C=linspecer(5);
NumComplexities=500;
rho=0*etas;
for jj=1:length(etas)
    etaPool=0*ones(PoolNumPlasmid,PoolNumSpecies,PoolNumSpecies);
    for i=1:PoolNumPlasmid
        for j=1:PoolNumSpecies
            etaPool(i,j,:)=etas(jj)*rand;          
        end
    end

    kappamax=0.02;
    kappaPool=kappamax*rand(PoolNumSpecies,PoolNumPlasmid);
    D=0.05;
    lambdaPool=0*ones(PoolNumSpecies,PoolNumPlasmid);
    for i=1:PoolNumSpecies
        lambdaPool(i,:)=0.1*rand;
    end
    sigmaPool=0.05*rand(PoolNumSpecies,1);
    gammaPool=0.02-0.04*rand(PoolNumSpecies,PoolNumSpecies);%(0.05-0.1*rand(PoolNumSpecies,PoolNumSpecies));
    temp=rand(PoolNumSpecies,PoolNumSpecies);
    gammaPool=gammaPool.*(temp>0.9);
    for jk=1:PoolNumSpecies
        gammaPool(jk,jk)=-1;
    end
    muPool=0.1+0.4*rand(PoolNumSpecies,1);
    NmPool=0.2+0.8*rand(PoolNumSpecies,1);
    
    CommunityComplexities=0*ones(1,NumComplexities);
    plasmid50=0*ones(NumComplexities,NumPlasmid);
    plasmid100=0*ones(NumComplexities,NumPlasmid);
    plasmid200=0*ones(NumComplexities,NumPlasmid);
    plasmid500=0*ones(NumComplexities,NumPlasmid);
    plasmid600=0*ones(NumComplexities,NumPlasmid);
    plasmid700=0*ones(NumComplexities,NumPlasmid);
    plasmid800=0*ones(NumComplexities,NumPlasmid);
    plasmid900=0*ones(NumComplexities,NumPlasmid);
    plasmid1000=0*ones(NumComplexities,NumPlasmid);
    for i=1:length(CommunityComplexities)
    NumSpecies=max(fix(rand*PoolNumSpecies),2);
    i
            rr=randperm(PoolNumSpecies);
            MapToPool=rr(1:NumSpecies);

            kappa=kappaPool(MapToPool,1:NumPlasmid);

            eta=etaPool(1:NumPlasmid,MapToPool, MapToPool);
            lambda=lambdaPool(MapToPool,1:NumPlasmid); 
            sigma=sigmaPool(MapToPool);
            gamma=gammaPool(MapToPool,MapToPool);
            mu=muPool(MapToPool);
            Nm=NmPool(MapToPool);
            Nm=Nm/sum(Nm);
            initial=0*ones(NumSpecies+NumSpecies*NumPlasmid,1);
            temp=rand(1,NumSpecies);
            initial(1:NumSpecies)=temp/sum(temp);
            for k=1:NumSpecies
                    for kk=1:NumPlasmid
                            initial(NumSpecies+(k-1)*NumPlasmid+kk)=0.01*initial(k);
                    end
            end
            [t,y]=ode45(@multi_plasmid,timespan,initial);
            s=y(end,1:NumSpecies);
            temp=s/sum(s);
            CommunityComplexities(i)=exp(sum(-temp.*log(temp)));
            for k=1:NumPlasmid
                    plasmid50(i,k)=sum(y(51,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(51,1:NumSpecies));
                    plasmid100(i,k)=sum(y(101,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(101,1:NumSpecies));
                    plasmid200(i,k)=sum(y(201,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(201,1:NumSpecies));
                    plasmid500(i,k)=sum(y(501,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(501,1:NumSpecies));
                    plasmid600(i,k)=sum(y(601,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(601,1:NumSpecies));
                    plasmid700(i,k)=sum(y(701,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(701,1:NumSpecies));
                    plasmid800(i,k)=sum(y(801,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(801,1:NumSpecies));
                    plasmid900(i,k)=sum(y(901,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(y(901,1:NumSpecies));
                    plasmid1000(i,k)=sum(y(end,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(s);
            end
    end
    
    X_input=CommunityComplexities;
    Y_input50=sum(plasmid50,2);
    Y_input100=sum(plasmid100,2);
    Y_input200=sum(plasmid200,2);
    Y_input500=sum(plasmid500,2);
    Y_input600=sum(plasmid600,2);
    Y_input700=sum(plasmid700,2);
    Y_input800=sum(plasmid800,2);
    Y_input900=sum(plasmid900,2);
    Y_input1000=sum(plasmid1000,2);
    rho50(jj)=corr(X_input',Y_input50,'type','pearson');
    rho100(jj)=corr(X_input',Y_input100,'type','pearson');
    rho200(jj)=corr(X_input',Y_input200,'type','pearson');
    rho500(jj)=corr(X_input',Y_input500,'type','pearson');
    rho600(jj)=corr(X_input',Y_input600,'type','pearson');
    rho700(jj)=corr(X_input',Y_input700,'type','pearson');
    rho800(jj)=corr(X_input',Y_input800,'type','pearson');
    rho900(jj)=corr(X_input',Y_input900,'type','pearson');
    rho1000(jj)=corr(X_input',Y_input1000,'type','pearson');
end
subplot(1,2,1);
plot(etas,rho1000,'.-','markersize',30,'color',C(5,:));
set(gca,'fontsize',16);
xlabel('mean transfer rate','fontsize',24);
ylabel('correlation coefficient','fontsize',24);
subplot(1,2,2);
plot(etas,rho50,'.-','markersize',30,'color',C(1,:));hold on;
plot(etas,rho100,'.-','markersize',30,'color',C(2,:));hold on;
plot(etas,rho200,'.-','markersize',30,'color',C(3,:));hold on;
plot(etas,rho500,'.-','markersize',30,'color',C(4,:));hold on;
plot(etas,rho1000,'.-','markersize',30,'color',C(5,:));hold on;
set(gca,'fontsize',16);
legend('50','100','200','500','1000');
legend boxoff;
set(gcf,'position',[100 100 700 300]);
save('MultiPlasmids_Rho.mat');
saveas(gcf,'MultiPlasmids_Rho.fig');
saveas(gcf,'MultiPlasmids_Rho.png');

function dydt=multi_plasmid(t,y)
        global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu Nm;
        dydt(NumSpecies*(1+NumPlasmid),1)=0;
        for i=1:NumSpecies
            sum=0;
            for j=1:NumPlasmid
            sum=sum+lambda(i,j)*y(NumSpecies+NumPlasmid*(i-1)+j);
            end
            Neg=0;
            Pos=0;
            for j=1:NumSpecies
                if gamma(j,i)<0
                Neg=Neg-gamma(j,i)*y(j);
                end
                if gamma(j,i)>0
                Pos=Pos+gamma(j,i)*y(j);
                end
            end
            
            mui=mu(i)*(1-Neg/Nm(i)+sigma(i)*Pos/(Pos+1)/Nm(i));
            if y(i)==0
                dydt(i,1)=0;
            else
                dydt(i,1)=y(i)/(sum+y(i))*mui*y(i)-D*y(i);
            end
            
            for j=1:NumPlasmid
                if y(i)==0
                    betaij=0;
                elseif sum==0
                    betaij=1;
                else
                    betaij=(1+lambda(i,j))/(1+lambda(i,j)+(sum-lambda(i,j)*y(NumSpecies+NumPlasmid*(i-1)+j))/y(i));
                end
                muij=mui/(1+lambda(i,j));
                summ=0;
                for k=1:NumSpecies
                    summ=summ+eta(j,k,i)*y(NumSpecies+(k-1)*NumPlasmid+j);
                end
                dydt(NumSpecies+(i-1)*NumPlasmid+j,1)=betaij*muij*y(NumSpecies+(i-1)*NumPlasmid+j)+(y(i)-y(NumSpecies+(i-1)*NumPlasmid+j))*summ-(kappa(i,j)+D)*y(NumSpecies+(i-1)*NumPlasmid+j);
            end
        end      
end
