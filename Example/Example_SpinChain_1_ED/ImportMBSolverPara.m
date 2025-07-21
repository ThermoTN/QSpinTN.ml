function [ Para ] = ImportMBSolverPara( Para )

% // beta list for ED
beta_list = 0.0025.*2.^(0:1:15).*2^0.1;
for int = 0.3:0.2:0.9
    beta_list = [beta_list, 0.0025.*2.^(0:1:15).*2^int];
end
beta_list = sort(beta_list);
Para.beta_list = beta_list;

end
