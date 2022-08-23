clear;
close all;
clc;
global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu Nm;
PoolNumSpecies=100;
PoolNumPlasmid=10;
NumPlasmid=PoolNumPlasmid;
timespan=0:2000;
etas=[0.04];
C=linspecer(3);
NumComplexities=100;
rho=0*etas;
for jj=1:length(etas)
    etaPool=0*ones(PoolNumPlasmid,PoolNumSpecies,PoolNumSpecies);
    for i=1:PoolNumPlasmid
        etamax=etas(jj)*rand^2;
        for j=1:PoolNumSpecies
            etaPool(i,j,:)=etamax*(0.8+0.4*rand);          
        end
    end

    kappamax=0.02;
    kappaPool=kappamax*rand(PoolNumSpecies,PoolNumPlasmid);
    D=0.02;
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
    plasmid=0*ones(NumComplexities,NumPlasmid);
    for i=1:length(CommunityComplexities)
    NumSpecies=max(fix(rand^2*PoolNumSpecies),2);
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
                    plasmid(i,k)=sum(y(end,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))/sum(s);
                    mu_ave=0;
                    for ijk=1:NumSpecies
                        mu_ave=mu_ave+mu(ijk)*y(end,ijk)/sum(s);
                    end

                    kappa_ave=0;
                    for ijk=1:NumSpecies
                        kappa_ave=kappa_ave+kappa(ijk,k)*y(end,ijk)/sum(s);
                    end

                    lambda_ave=0;
                    for ijk=1:NumSpecies
                        lambda_ave=lambda_ave+lambda(ijk,k)*y(end,ijk)/sum(s);
                    end

                    etasum=0;
                    for ijk=1:NumSpecies
                        for jik=1:NumSpecies
                        etasum=etasum+eta(k,ijk,jik)*y(end,ijk)*y(end,jik)/(sum(s))^2;
                        end
                    end
                    omega(i,k)=etasum/(mu_ave/(mu_ave-D)*(D+kappa_ave-D/(1+lambda_ave)));
            end
    end
    
end
plot(omega,plasmid,'.','markersize',20,'color',C(3,:)); 
hold on;
x=1:0.01:6;
plot(x,1-1./x,'k-','linewidth',2);hold on;
set(gca,'fontsize',16);
xlabel('\omega','fontsize',24);
ylabel('plasmid abundance','fontsize',24);
set(gcf,'position',[100 100 300 300]);
axis([0 6 0 1]);
save('EvaluateOmega.mat');
saveas(gcf,'EvaluateOmega.fig');
saveas(gcf,'EvaluateOmega.png');

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
