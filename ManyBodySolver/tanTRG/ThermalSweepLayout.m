classdef ThermalSweepLayout
    % Explicit singleton boundaries avoid losing rank-one MPO legs.
    methods (Static)
        function [A,W] = pack(rho,H)
            L=numel(rho.A); d=size(rho.A{1},2); A=cell(L,1); W=cell(L,1);
            for i=1:L
                a=rho.A{i}; h=H.A{i};
                if i==1
                    a=reshape(a,1,[],d,d); h=reshape(h,1,[],d,d);
                elseif i==L
                    a=reshape(a,[],1,d,d); h=reshape(h,[],1,d,d);
                end
                A{i}=reshape(a,size(a,1),size(a,2),d*d);
                W{i}=ThermalSweepLayout.lift(h,d);
            end
        end
        function w = lift(h,d)
            w=zeros(size(h,1),size(h,2),d*d,d*d,'like',h);
            for ancilla=1:d
                idx=(ancilla-1)*d+(1:d); w(:,:,idx,idx)=h;
            end
        end
        function rho = unpack(rho,A,d)
            for i=1:numel(A)
                if i==1 || i==numel(A), rho.A{i}=reshape(A{i},[],d,d);
                else, rho.A{i}=reshape(A{i},size(A{i},1),size(A{i},2),d,d); end
            end
        end
    end
end
