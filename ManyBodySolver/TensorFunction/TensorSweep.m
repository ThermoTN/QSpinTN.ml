classdef TensorSweep
    % Dense sweep kernels: A(left,right,physical), W(left,right,out,in),
    % environments(bra,ket,MPO). Physical may also combine MPO ket/ancilla.
    methods (Static)
        function B = leftAction(left,A,W)
            B = tensorprod(left,A,2,1,NumDimensionsA=3);
            B = tensorprod(B,W,[2,4],[1,4],NumDimensionsA=4);
        end
        function E = extendLeft(left,A,W)
            E = tensorprod(conj(A),TensorSweep.leftAction(left,A,W),[1,3],[1,4],NumDimensionsA=3);
        end
        function E = extendRight(right,A,W)
            E = TensorSweep.extendLeft(right,permute(A,[2,1,3]),permute(W,[2,1,3,4]));
        end
        function apply = operator(left,W,right,shape)
            dl=shape(1); dr=shape(2); d=shape(3); wr=size(W,2);
            fused=tensorprod(left,W,3,1,NumDimensionsA=3);
            fused=reshape(permute(fused,[1,4,3,2,5]),dl*d*wr,dl*d);
            rm=reshape(permute(right,[2,3,1]),dr*wr,dr);
            apply=@(x)TensorSweep.fusedAction(x,fused,rm,dl,dr,d,wr);
        end
        function Y = fusedAction(A,fused,rm,dl,dr,d,wr)
            B=fused*reshape(permute(A,[1,3,2]),dl*d,dr);
            B=reshape(permute(reshape(B,dl*d,wr,dr),[1,3,2]),dl*d,dr*wr);
            Y=permute(reshape(B*rm,dl,d,dr),[1,3,2]);
        end
        function apply = twoSiteOperator(left,W1,W2,right,dl,dr)
            d=size(W1,3); w=size(W1,2);
            f=tensorprod(left,W1,3,1,NumDimensionsA=3);
            f=reshape(permute(f,[1,4,3,2,5]),dl*d*w,dl*d);
            g=tensorprod(right,W2,3,2,NumDimensionsA=3);
            g=reshape(permute(g,[2,5,3,1,4]),dr*d*w,dr*d);
            apply=@(x)TensorSweep.twoSiteAction(x,f,g,dl,dr,d,w);
        end
        function Y = twoSiteAction(A,f,g,dl,dr,d,w)
            matrix=reshape(permute(reshape(A,dl,dr,d,d),[1,3,2,4]),dl*d,dr*d);
            B=f*matrix;
            B=reshape(permute(reshape(B,dl*d,w,dr*d),[1,3,2]),dl*d,dr*d*w);
            Y=reshape(permute(reshape(B*g,dl,d,dr,d),[1,3,2,4]),dl,dr,d*d);
        end
        function [site,bond,logNorm,entropy,discarded] = split(center,direction,cap,tol)
            if direction<0, center=permute(center,[2,1,3]); end
            [dl,dr,d]=size(center);
            [U,S,V]=svd(reshape(permute(center,[1,3,2]),dl*d,dr),'econ');
            s=diag(S); total=sum(s.^2);
            assert(isfinite(total) && total>0,'TensorSweep:Norm','Zero or nonfinite center norm.');
            keep=max(1,min(cap,sum(s/sqrt(total)>tol)));
            discarded=sum(s(keep+1:end).^2)/total;
            p=s(s>0).^2/total; entropy=-sum(p.*log(p));
            ns=norm(s(1:keep)); logNorm=log(ns);
            site=permute(reshape(U(:,1:keep),dl,d,keep),[1,3,2]);
            bond=(s(1:keep)/ns).*V(:,1:keep)';
            if direction<0
                site=permute(site,[2,1,3]); bond=bond.';
            end
        end
        function next = absorb(next,bond,direction)
            if direction<0, next=permute(next,[2,1,3]); bond=bond.'; end
            next=reshape(bond*reshape(next,size(next,1),[]),size(bond,1),size(next,2),size(next,3));
            if direction<0, next=permute(next,[2,1,3]); end
        end
    end
end
