clc;clear all;
% #Leaders = 10
% #Followers = 14
tic

%% Initialization
rho = 5;                 %Penalty parameter 
tol = 0.001;             %Residual tolerance
BWL = 10000;                 %LWB iteration count 
n=24;                    %Total agents
m=10;                    %Leaders
for i=1:n-m
    M(:,:,i) = zeros(n-m,n-m);
    W(:,:,i) = zeros(n-m,n-m);
    Z(:,:,i) = zeros(n-m,n-m);
end
k=1;
R(1)=1;
Resid_ag = zeros(n-m,1);
norm_weight = zeros(n-m,1);

%% Agent State Initialization
SW = ones(3,n);
SW(1:2,:) = [-4 + (17).*rand(1,14) 5 3 2 2 3 5 7 8 8 7;-4 + (17).*rand(1,14) 1 2 4 5 7 8 7 5 4 2];
initx = SW(1,:);
inity = SW(2,:);
Bx = initx;
By = inity;
x_bar = mean(initx(n-m+1:n))*ones(1,n-m);
y_bar = mean(inity(n-m+1:n))*ones(1,n-m);

%% Topology
FT = ones(n-m,n-m);
for i=3:12
    FT(4,i)=0;
end
for i=2:2:13
    FT(5,i)=0;
end
for i=3:2:13
    FT(9,i)=0;
end
for i=4:13
    FT(10,i)=0;
end
for i=2:9
    FT(11,i)=0;
end
for i=1:13
    FT(12,i)=0;
end
for i=2:14
    FT(13,i)=0;
end
for i=1:n-m
    FT(i,i)=1;
end
d = sum(FT,1)-1;                                        %Out-degree
adj_sparse = sparse(FT-eye(n-m));
distances = graphallshortestpaths(adj_sparse);
D = max(distances(:));                                  %Diameter
%w = ((1/max(d))^(2*D+1))*ones(n-m,1);
w = 0.01 + 0.1*rand(n-m,1);

LT = zeros(n-m,m);
LT(1,18-14)=1; LT(2,17-14)=1; LT(3,19-14)=1; LT(4,15-14)=1; LT(5,16-14)=1; 
LT(6,21-14)=1; LT(7,20-14)=1; LT(10,24-14)=1; LT(12,23-14)=1; LT(14,22-14)=1;
        
%% Distributed-Optimized-Containment
while R(k)>tol 
    r=zeros(n-m,1);
    r1=zeros(n-m,n-m);
    r2=zeros(n-m,n-m);
        
    %% ADMM
    %Primal Update (W)
    for i=1:n-m
        cvx_begin quiet 
        variables C(n-m,n-m)    
        expressions cost
        cost = norm(C-(1/(n-m))*ones(n-m,n-m),2)/(n-m) + trace(C'*M(:,:,i)) +... 
        (rho/2)*sum(sum((C-Z(:,:,i)).*(C-Z(:,:,i))));
        minimize cost
        subject to 
        for g=1:n-m
            if FT(i,g)==0
                C(i,g)==0;
            end
        end
        C*ones(n-m,1)==ones(n-m,1);
        C'*ones(n-m,1)==ones(n-m,1);
        cvx_end
        W(:,:,i) = C;
    end
    %Primal Update (Z) [Linear Weight Balance]
    WS = W(:,:,1);
    l=1;
    for i=2:n-m
        WS = [WS;W(:,:,i)];
    end
    while l<=BWL
        Wt = eye(n-m) - (diag(d)-FT+eye(n-m))*diag(w(:,(k-1)*BWL+l));
        WS = kron(Wt,eye(n-m))*WS;
        w(:,(k-1)*BWL+l+1) = 0.5*(eye(n-m) + inv(diag(d))*(FT-eye(n-m)))*w(:,(k-1)*BWL+l);
        l=l+1;
    end
    for o=1:n-m
        Z(:,:,o) = WS((o-1)*(n-m)+1:o*(n-m),:);
    end
    %Dual Update (M)
    for i=1:n-m
        M(:,:,i) = M(:,:,i) + rho*(W(:,:,i)-Z(:,:,i));
    end
      
    %% Stopping Criterion
    for i=1:n-m
        for h=1:n-m
            r1(i,h) = (1/(n-m))*norm(W(:,:,i)-W(:,:,h),'fro')*FT(i,h);
            r2(i,h) = norm(W(i,h,i),1)*(1-FT(i,h));
        end
    end
    for i=1:n-m
        r(i) = max(max(r1(i,:)),max(r2(i,:)));
    end
    Resid_ag(:,k) = r;
    for i=1:n-m
        norm_weight(i,k) = norm(W(:,:,i)-(1/(n-m))*ones(n-m,n-m),2);
    end
    k=k+1;
    R(k) = max(r);
    [k;R(k);norm_weight(1,k-1)]
    toc
end

%% Protocol Push-Sum
%Row extraction
A = zeros(n-m,n-m);
for i=1:n-m
    A(i,:) = W(i,:,i);
end
itr=1;
max_iter=500;
A = [A LT];
while norm(Bx(itr,1:n-m)-x_bar,2)>0.05 && itr<=max_iter
    for g=1:n-m
        S = [0;0;0];
        for h=1:n
            S = S + A(g,h)*SW(:,h);
        end
        SW(:,g) = S;
    end
    itr=itr+1;
    Bx(itr,:) = SW(1,:)./SW(3,:);
    By(itr,:) = SW(2,:)./SW(3,:);
end

%% Laplacian and Perron 
A = [A;zeros(m,n-m) eye(m)];
L = zeros(n,n);
D = zeros(n,n);
Ai = zeros(n,n);
A1 = zeros(n-m,n-m);
A2 = zeros(n-m,m);
for i=1:size(A,1)
    D(i,i) = sum(A(i,:));
end
for i=1:n-m
    for j=1:n
        if j<=n-m
            A1(i,j) = A(i,j);
        else
            A2(i,j-n+m) = A(i,j);
        end
    end
end
L = D-A;
Ai = D-L;
af = max(abs(real(eig(Ai(1:n-m,1:n-m) - (1/(n-m))*ones(n-m,n-m)))));

%% Topology Plot
xy = [ 4.5 1 0 0 1 4.5 8 9 9 8 3 6 3 6 4.5,2,1,1,2,4.5,7,8,8,7 ; ...
      0 2 4 5 7 9 7 5 4 2 6 6 3 3 1 2 4 5 7 8 7 5 4 2]';
figure(1)
gplot(A,xy,'--dk')
hold on
plot([4.5,2,1,1,2,4.5,7,8,8,7],[1,2,4,5,7,8,7,5,4,2],'dr',4.5,0,'db',1,2,'db',0,4,'db',0,5,'db',1,7,'db',4.5,9,'db',8,7,'db',9,5,'db',9,4,'db',8,2,'db',3,6,'db',6,6,'db',3,3,'db',6,3,'db')
xlabel('x-coordinates');
ylabel('y-coordinates');
title('Graph depicting Adjacency Matrix(A)');
legend('Links','Leader','Follower');
axis square
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
print -depsc grp10PS.eps

%% Performance Plot
mx=max(max(Ai));
mn=min(min(Ai));
ne=0;                                    %Counting optimized #edges
for i=1:n-m
    for j=1:n
        if A(i,j) > 1e-4
            ne=ne+1;
        end
    end
end

figure(2)
for i=1:n
    plot([1:itr],Bx(:,i)')
    hold on
end
xlabel(['k'],'interpreter','latex','FontWeight','bold');
ylabel(['$x_{k}$'],'interpreter','latex','FontWeight','bold');
title(['Single Integrator Containment (10L/14F) (Weights \in [',num2str(mn),',',num2str(mx),']) (Fiedlar Eigenvalue of A_{ff} = ',num2str(af),') (#edges = ',num2str(ne),')']);
key1 = {'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','L1','L2','L3','L4','L5','L6','L7','L8','L9','L10'};
legend(key1,'Location','South','Fontsize',6);
print -depsc pplot10PS.eps

figure(3)
for i=1:itr
    diff_Sx(i) = norm(Bx(i,1:n-m)-x_bar,2);
    diff_Sy(i) = norm(By(i,1:n-m)-y_bar,2);
end
plot([1:itr],diff_Sx,[1:itr],diff_Sy)
xlabel(['k'],'interpreter','latex','FontWeight','bold');
ylabel(['$||x_{k}-\overline{x}||_{2}$'],'interpreter','latex','FontWeight','bold');
title(['Optimization performance curve (Distributed)'],'interpreter','latex');
legend('X','Y');
print -depsc opplot10PS.eps

figure(4)
plot([1:k],log(R))
xlabel(['k'],'interpreter','latex','FontWeight','bold','FontSize',20);
ylabel(['$log_eR(k)$'],'interpreter','latex','FontWeight','bold','FontSize',20);
%title(['Maximum Residual value'],'interpreter','latex');
print -depsc rplot10PS.eps

W_wb = eye(n-m) - (diag(d)- FT + eye(n-m))*diag(w(:,end));
slem = max(abs(real(eig(W_wb - (1/(n-m))*ones(n-m,n-m)))));
figure(5)
for i=1:n-m
    p1=plot([1:k-1],norm_weight(i,:),'-b','LineWidth',2);                          %ADMM
    hold on
end
p2=plot([1:k-1],0.7071*ones(1,k-1),'--r','LineWidth',2.7);                      %WB
hold on
p3=plot([1:k-1],slem*ones(1,k-1),'-.k','LineWidth',2.7);                      %Centralized
xlabel(['k'],'interpreter','latex','FontWeight','bold','FontSize',20);
ylabel(['$|| (A_1)_{i}(k)-\frac{1}{n-m} 11^{T}||_{2}$'],'interpreter','latex','FontWeight','bold','FontSize',20);
%title(['Objective variation plot (Directed Followers)'],'interpreter','latex');
%legend({'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','WBA','Centralized A_1^{*}'},'Location','Southeast','FontSize',8);
lgd=legend([p1,p2,p3],'Agents estimate','Centralized A_1^{*}','WBA','Location','Northeast');
lgd.FontSize=13;
print -depsc ovplot10PS.eps

toc
                  