function [rho,info] = TwoSiteTanTRGSweep(rho,H,Para,env)
% Two-site projector splitting with explicit boundary dimensions.
L=Para.L; q=Para.d^2; delta=Para.tau*exp(H.lgnorm);
[A,W]=ThermalSweepLayout.pack(rho,H);
info=struct('TruncationError',0,'Entropy',0,'AddedStates',0, ...
    'ExpansionSeconds',0,'MaxExpandedBond',0);
for direction=[1,-1]
    if direction>0, bonds=1:L-1; else, bonds=L-1:-1:1; end
    rolling=1;
    for b=bonds
        if direction>0
            left=rolling;
            if b+1==L, right=1; else, right=permute(env.get(b+2,'R'),[1,3,2]); end
        else
            right=rolling;
            if b==1, left=1; else, left=permute(env.get(b-1,'L'),[1,3,2]); end
        end
        dl=size(A{b},1); dr=size(A{b+1},2);
        pair=tensorprod(A{b},A{b+1},2,1,NumDimensionsA=3);
        pair=reshape(permute(pair,[1,3,2,4]),dl,dr,q*q);
        apply=TensorSweep.twoSiteOperator(left,W{b},W{b+1},right,dl,dr);
        pair=evolve(apply,pair,-delta);
        matrix=reshape(permute(reshape(pair,dl,dr,q,q),[1,3,2,4]),dl*q,dr*q);
        [U,S,V]=svd(matrix,'econ'); s=diag(S); total=sum(s.^2);
        k=max(1,min(Para.MCrit,sum(s/sqrt(total)>Para.TSRGStol)));
        ns=norm(s(1:k)); rho.lgnorm=rho.lgnorm+log(ns);
        info.TruncationError=max(info.TruncationError,sum(s(k+1:end).^2)/total);
        p=s(s>0).^2/total; info.Entropy=max(info.Entropy,-sum(p.*log(p)));
        U=U(:,1:k); V=V(:,1:k)';
        if direction>0, V=(s(1:k)/ns).*V; else, U=U.*(s(1:k).'/ns); end
        A{b}=permute(reshape(U,dl,q,k),[1,3,2]); A{b+1}=reshape(V,k,dr,q);
        if direction>0
            rolling=TensorSweep.extendLeft(left,A{b},W{b}); env.put(b,'L',permute(rolling,[1,3,2]));
            center=b+1; left=rolling;
        else
            rolling=TensorSweep.extendRight(right,A{b+1},W{b+1}); env.put(b+1,'R',permute(rolling,[1,3,2]));
            center=b; right=rolling;
        end
        if (direction>0 && b<L-1) || (direction<0 && b>1)
            apply=TensorSweep.operator(left,W{center},right,size(A{center}));
            A{center}=evolve(apply,A{center},delta);
            ns=norm(A{center}(:)); A{center}=A{center}/ns; rho.lgnorm=rho.lgnorm+log(ns);
        end
    end
end
rho=ThermalSweepLayout.unpack(rho,A,Para.d);
end

function y = evolve(apply,x,delta)
[y,~]=LanczosExp(apply,x,delta,'innerprod',@(a,b)a(:)'*b(:),'norm',@(a)norm(a(:)));
end
