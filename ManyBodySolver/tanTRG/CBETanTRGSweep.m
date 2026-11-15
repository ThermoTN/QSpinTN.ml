function [rho,info] = CBETanTRGSweep(rho,H,Para,env)
% One-site projector-splitting evolution with RSVD pre-expansion and SVD cap.
% Fuse MPO physical output/input as an MPS physical index; H acts on output.
L=Para.L; d=Para.d; delta=Para.tau*exp(H.lgnorm);
cbe=string(Para.TanTRGUpdateSites)=="cbe";
if cbe, opts=GetCBEOptions(Para); end
[A,W]=ThermalSweepLayout.pack(rho,H);
info=struct('TruncationError',0,'Entropy',0,'AddedStates',0, ...
    'ExpansionSeconds',0,'MaxExpandedBond',0);
for direction=[1,-1]
    if direction>0, sites=1:L; rolling=1; else, sites=L:-1:1; rolling=1; end
    for i=sites
        j=i+direction;
        if j>=1 && j<=L
            outerIndex=i+2*direction;
            if outerIndex<1 || outerIndex>L
                outer=1;
            elseif direction>0
                outer=permute(env.get(outerIndex,'R'),[1,3,2]);
            else
                outer=permute(env.get(outerIndex,'L'),[1,3,2]);
            end
            timer=tic;
            n=0;
            if cbe
                [A{i},A{j},n]=CBEExpandBond(A{i},A{j},rolling,W{i},W{j},outer,opts,direction);
            end
            if direction>0
                if cbe, opposite=TensorSweep.extendRight(outer,A{j},W{j});
                else, opposite=permute(env.get(j,'R'),[1,3,2]); end
                left=rolling; right=opposite; bondSize=size(A{i},2);
            else
                if cbe, opposite=TensorSweep.extendLeft(outer,A{j},W{j});
                else, opposite=permute(env.get(j,'L'),[1,3,2]); end
                left=opposite; right=rolling; bondSize=size(A{i},1);
            end
            if cbe
                info.ExpansionSeconds=info.ExpansionSeconds+toc(timer);
                info.AddedStates=info.AddedStates+n;
                info.MaxExpandedBond=max(info.MaxExpandedBond,bondSize);
            end
        elseif direction>0
            left=rolling; right=1;
        else
            left=1; right=rolling;
        end
        apply=TensorSweep.operator(left,W{i},right,size(A{i}));
        center=evolve(apply,A{i},-delta);
        if j>=1 && j<=L
            if cbe, cap=Para.MCrit; tol=Para.TSRGStol; else, cap=inf; tol=0; end
            [A{i},bond,ln,ee,te]=TensorSweep.split(center,direction,cap,tol);
            rho.lgnorm=rho.lgnorm+ln;
            info.TruncationError=max(info.TruncationError,te); info.Entropy=max(info.Entropy,ee);
            if direction>0
                rolling=TensorSweep.extendLeft(rolling,A{i},W{i});
                env.put(i,'L',permute(rolling,[1,3,2]));
                left=rolling; right=opposite;
            else
                rolling=TensorSweep.extendRight(rolling,A{i},W{i});
                env.put(i,'R',permute(rolling,[1,3,2]));
                left=opposite; right=rolling;
            end
            bond=evolve(@(b)bondAction(b,left,right),bond,delta);
            % Keep log amplitude explicit before absorbing the backward bond step.
            ns=norm(bond,'fro'); bond=bond/ns; rho.lgnorm=rho.lgnorm+log(ns);
            A{j}=TensorSweep.absorb(A{j},bond,direction);
        else
            ns=norm(center(:)); A{i}=center/ns; rho.lgnorm=rho.lgnorm+log(ns);
        end
    end
end
rho=ThermalSweepLayout.unpack(rho,A,d);
end

function y = bondAction(b,left,right)
y=tensorprod(left,b,2,1,NumDimensionsA=3);
y=tensorprod(y,right,[2,3],[3,2],NumDimensionsA=3);
end

function y = evolve(apply,x,delta)
[y,~]=LanczosExp(apply,x,delta,'innerprod',@(a,b)a(:)'*b(:),'norm',@(a)norm(a(:)));
end
