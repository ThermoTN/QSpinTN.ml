function [center,next,added] = CBEExpandBond(center,next,left,W,Wnext,right,opts,direction)
% RSVD pre-expansion of the doubly discarded two-site H|psi> component.
% Factor X*Y avoids forming the full two-site tensor. The original state is
% unchanged: append orthogonal rows to the neighbor and zeros to the center.
% For a leftward move, left/right denote the outer environments in sweep order.
if direction<0
    center=permute(center,[2,1,3]); next=permute(next,[2,1,3]);
    W=permute(W,[2,1,3,4]); Wnext=permute(Wnext,[2,1,3,4]);
end
[dl,m,d]=size(center); dr=size(next,2); q=size(next,3);
added=0;
k=min([opts.Cap-m,dl*d-m,dr*q-m]);
if opts.Ratio>0 && k>0
    cm=reshape(permute(center,[1,3,2]),dl*d,m);
    [U,~]=qr(cm,0); B=reshape(next,m,dr*q);
    X=TensorSweep.leftAction(left,center,W);
    X=reshape(permute(X,[1,4,2,3]),dl*d,[]);
    Y=TensorSweep.leftAction(right,permute(next,[2,1,3]),permute(Wnext,[2,1,3,4]));
    Y=reshape(permute(Y,[2,3,1,4]),[],dr*q);
    scale=norm(X,'fro')*norm(Y,'fro');
    X=X-U*(U'*X); Y=Y-(Y*B')*B;
    sample=min([k+opts.Oversampling,dl*d-m,dr*q-m]);
    omega=randn(dl*d,sample,'like',X);
    if ~isreal(X) || ~isreal(Y), omega=(omega+1i*randn(size(omega),'like',omega))/sqrt(2); end
    Z=Y'*(X'*omega);
    for iteration=1:opts.PowerIterations
        [Q,~]=qr(Z,0); [P,~]=qr(X*(Y*Q),0); Z=Y'*(X'*P);
    end
    [Q,~]=qr(Z,0);
    [~,S,V]=svd(X*(Y*Q),'econ'); s=diag(S);
    % Ignore numerical residue when H has no component outside the kept spaces.
    count=min(k,sum(s>max(opts.SvdTol*max(s),64*eps(class(s))*scale)));
    if count>0
        extra=Q*V(:,1:count);
        extra=extra-B'*(B*extra); extra=extra-B'*(B*extra);
        [extra,~]=qr(extra,0);
        next=reshape([B;extra'],m+count,dr,q);
        center=cat(2,center,zeros(dl,count,d,'like',center));
        added=count;
    end
end
if direction<0
    center=permute(center,[2,1,3]); next=permute(next,[2,1,3]);
end
end
