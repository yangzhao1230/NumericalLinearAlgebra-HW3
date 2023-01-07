close all
clear all
clc
v0=-ones(32,1);
v1=ones(31,1);
v2=ones(30,1);
A=diag(v0,0)+diag(v1,1)+diag(v2,2);
[U,S,V]=svd(A);
npts=20;
s=.8*norm(A,1);
xmin=-s;xmax=s;ymin=-s;ymax=s;
x=xmin:(xmax-xmin)/(npts-1):xmax;
y=ymin:(ymax-ymin)/(npts-1):ymax;
[xx,yy]=meshgrid(x,y);
zz=xx+sqrt (-1)*yy;
[U,T]=schur(A);
if isreal(A),
    [U,T]=rsf2csf(U,T);
end
T=triu(T);
eigA=diag(T);

hold off
plot(real(eigA),imag(eigA),'.','markersize',15)
hold on
axis([xmin xmax ymin ymax]), axis square, grid on, drawnow
select=find(real(eigA)>-250);
n=length(select);
for i=1:n,
    for k=select(i)-1:-1:i,
        G([2 1],[2 1])=planerot([T(k,k+1) T(k,k)-T(k+1,k+1)]')';
        J=k:k+1;
        T(:,J)=T(:,J)*G;
        T(J,:)=G'*T(J,:);
    end
end
T=triu(T(1:n,1:n));
I=eye(n);

clear t
t=[0:0.1:50];
for k=1:length(t),
    y(k)=norm(exp(t(k).*A),2);
end
figure
grid on
semilogy(t,y)
xlabel('t')
ylabel('||e^{tA}||_2');
sigmin=Inf*ones(length(y),length(x));
for i=1:length(y),
    if isreal(A) & (ymax==-ymin) & (i>length(y)/2),
        sigmin(i,:)=sigmin(length(y)+1-i,:);
    else
        for j=1:length(x),
            z=zz(i,j);
            T1=z*I-T;
            T2=T1';
            if real(z)<100
                sigold=0;
                qold=zeros(n,1);
                beta=0;
                H=[];
                q=randn(n,1)+sqrt(-1)*randn(n,1);
                q=q/norm(q);
                for k=1:99,
                    v=T1\(T2\q)-beta*qold;
                    alpha=real(q'*v);
                    v=v-alpha*q;
                    beta=norm(v);
                    qold=q;
                    q=v/beta;
                    H(k+1,k)=beta;
                    H(k,k+1)=beta;
                    H(k,k)=alpha;
                    if alpha>1e100,
                        sig=alpha;
                    else
                        sig = max(eig(H(1:k,1:k)));
                    end
                    if (abs(sigold/sig-1)<.001) |  (sig<3 & k>2),
                        break
                    end
                    
                    sigold = sig;
                end
                sigmin(i,j) = 1/sqrt(sig);
            end
        end
    end
end
contour(x,y,log10(sigmin+1e-20),-8:-1);

