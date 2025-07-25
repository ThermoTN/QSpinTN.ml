function [ Lattice, ModelConf ] = SpinModel_TLXXZ_Spin1( )
% -------------------------------------------------------------
% Triangular lattice
% XXZ model
% Parameter:
%           J1xy        Nearest neighbor SxSx+SySy term
%           J1z         Nearest neighbor SzSz term
%           J2xy        Next-nearest neighbor SxSx+SySy term
%           J2z         Next-nearest neighbor SzSz term
%           gx          Lande factor of Sx direction
%           gz          Lande factor of Sz direction
%
% Hamiltonian:
%   H = \sum_<i,j> J1xy (Sx_i Sx_j + Sy_i Sy_j) + J1z Sz_i Sz_j
%       + \sum_<<i,j>> J1xy (Sx_i Sx_j + Sy_i Sy_j) + J2z Sz_i Sz_j
%       - h\sum_i Sh_i
% -------------------------------------------------------------

% =============================================================
% DEFAULT SETTINGS
% =============================================================
ModelConf.ModelName = 'TLXXZ_Spin1';
ModelConf.ModelName_Full = 'Triangular lattice-XXZ model';
ModelConf.IntrcMap = 'IntrcMap_TLXXZ_Spin1';
ModelConf.AvlbSolver = {'ED', 'XTRG', 'tTRG'};    % can be solved by ED (high-T) and XTRG (low-T)
ModelConf.LocalSpin = '1';
% Parameters' name in Hamiltonian
ModelConf.Para_Name = {'J1xy'; 'J1z'; 'J2xy'; 'J2z'; 'Dz'};
% Parameters' unit
%  'K':   Set the unit as Kelvin.
%  'ES':  Set the unit as energy scale, i.e. the correspoding term
%         interaction strength is [*] times of energy scale.
ModelConf.Para_Unit = {'K'; 'K'; 'K'; 'K'; 'K'};
% Energy scale of this model, choose one from ModelConf.Para_Name
ModelConf.Para_EnScale = 'J1z';
ModelConf.Para_Range = cell(length(ModelConf.Para_Name), 1);

ModelConf.gFactor_Num = 3;
% Choose the type of Lande g factor
%  'xyz': The direction of g factor is given along Sx, Sy, or Sz as [1,0,0], [0,1,0], [0,0,1] respectively.
%  'dir': The direction of g factor is given along a chosen direction as [x, y, z] where x^2+y^+z^2 = 1;
%         under this setting, the CmDatagInfo and ChiDatagInfo is required
%         when importing the data to decide use which g factor for conversion.
ModelConf.gFactor_Type = 'xyz';
ModelConf.gFactor_Name = cell(ModelConf.gFactor_Num, 1);
ModelConf.gFactor_Vec = cell(ModelConf.gFactor_Num, 1);
ModelConf.gFactor_Range = cell(ModelConf.gFactor_Num, 1);
ModelConf.gFactor_Name{1} = 'gx';
ModelConf.gFactor_Vec{1} = [1,0,0];
ModelConf.gFactor_Name{2} = 'gy';
ModelConf.gFactor_Vec{2} = [0,1,0];
ModelConf.gFactor_Name{3} = 'gz';
ModelConf.gFactor_Vec{3} = [0,0,1];

% =============================================================
% LATTICE GEOMETRY SETTINGS
% =============================================================
Lattice.Lx = 6;
Lattice.Ly = 4;
Lattice.BCX = 'OBC';
Lattice.BCY = 'PBC';
Lattice.L = Lattice.Lx * Lattice.Ly;

% =============================================================
% PARAMETERS OPTIMIZATION RANGE SETTINGS
% =============================================================
% J1xy range
ModelConf.Para_Range{1} = [0, 5];
% J1z range
ModelConf.Para_Range{2} = [0.5, 5];
% J2xy range
ModelConf.Para_Range{3} = 0;
% J2z range
ModelConf.Para_Range{4} = 0;
% D range
ModelConf.Para_Range{5} = [-5,5];
% gx range
ModelConf.gFactor_Range{1} = 2;
% gy range
ModelConf.gFactor_Range{2} = 2;
% gz range
ModelConf.gFactor_Range{3} = 2;
end

