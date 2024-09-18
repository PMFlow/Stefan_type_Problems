function [c]=GRW_1D(c0,I,dx,dt,q,D,d,stepU)

%% (unbiased) GRW solution
N=10^24; n0=round(c0*N);
u=floor(q*stepU+0.5);
r=2*D*dt/d^2/dx^2*ones(1,I);
n=n0; nn=zeros(1,I);
restr=0; restjump=0;

for x=1:I
    if n(x) > 0
        xa=x+u(x);
        restr=n(x)*(1-r(x))+restr; nsta=floor(restr);
        restr=restr-nsta; njump=n(x)-nsta;
        if xa<1
            nn(1)=nn(1)+n(x);
        elseif xa>I
            nn(I)=nn(I)+n(x);
        else
            nn(xa)=nn(xa)+nsta;
            if(njump)>0
                restjump=njump/2+restjump;
                nj(1)=floor(restjump); restjump=restjump-nj(1);
                nj(2)=njump-nj(1);
                if xa==1
                    nn(2)=nn(2)+nj(2);
                elseif xa==I
                    nn(I-1)=nn(I-1)+nj(1);
                else
                    for i=1:2
                        if nj(i)>0
                            xj=xa+(2*i-3)*d;
                            nn(xj)=nn(xj)+nj(i);
                        end
                    end
                end
            end
        end
    end
end
c=nn/N;
