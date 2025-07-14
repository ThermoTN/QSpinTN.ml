function [ t_ab, t_ba ] = Trotter_DTNC( Para )

t_ab = cell(3,3);

t_ab{1,1} = 'Sx';
t_ab{1,2} = 'Sx';
t_ab{1,3} = Para.Model.Jxy;

t_ab{2,1} = 'Sy';
t_ab{2,2} = 'Sy';
t_ab{2,3} = Para.Model.Jxy;

t_ab{3,1} = 'Sz';
t_ab{3,2} = 'Sz';
t_ab{3,3} = Para.Model.Jz;

t_ab{4,1} = 'Sz2';
t_ab{4,2} = 'Id';
t_ab{4,3} = -Para.Model.Delta/2;

t_ab{5,1} = 'Id';
t_ab{5,2} = 'Sz2';
t_ab{5,3} = -Para.Model.Delta/2;

t_ba{1,1} = 'Sx';
t_ba{1,2} = 'Sx';
t_ba{1,3} = Para.Model.Jxy;

t_ba{2,1} = 'Sy';
t_ba{2,2} = 'Sy';
t_ba{2,3} = Para.Model.Jxy;

t_ba{3,1} = 'Sz';
t_ba{3,2} = 'Sz';
t_ba{3,3} = Para.Model.Jz;

t_ba{4,1} = 'Sz2';
t_ba{4,2} = 'Id';
t_ba{4,3} = -Para.Model.Delta/2;

t_ba{5,1} = 'Id';
t_ba{5,2} = 'Sz2';
t_ba{5,3} = -Para.Model.Delta/2;

end

