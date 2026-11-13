function tests = test_thermal_energy_scale
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.old_path=path;
addpath(fullfile(root,'UtilityFunc'));
end

function teardownOnce(testCase)
path(testCase.TestData.old_path);
end

function testAbsolutePhysicalCoefficients(testCase)
result=ResolveThermalEnergyScale([-2,1,0],[0,0,.1],[2,2,2]);
verifyEqual(testCase,result.energy_scale_K,2);
verifyFalse(testCase,result.zero_fallback_used);
end

function testAllFieldsAndPrincipalAxes(testCase)
u=QSpinThermalUnits(); axes=[0,1,0;1,0,0;0,0,1];
result=ResolveThermalEnergyScale(.1,[1,0,0;0,2,0;0,0,.1],[2,3,4],axes);
verifyEqual(testCase,result.field_energies_K,u.muB_over_kB_K_per_T*[3,4,.4],'AbsTol',1e-12);
verifyEqual(testCase,result.energy_scale_K,4*u.muB_over_kB_K_per_T,'AbsTol',1e-12);
end

function testVectorNormInsteadOfLargestComponent(testCase)
u=QSpinThermalUnits();
result=ResolveThermalEnergyScale(0,[1,1,0],[2,3,4]);
verifyEqual(testCase,result.energy_scale_K,sqrt(13)*u.muB_over_kB_K_per_T,'AbsTol',1e-12);
end

function testOnlyExactlyZeroUsesFallback(testCase)
zero=ResolveThermalEnergyScale([0,0],[0,0,0],[2,2,2]);
verifyEqual(testCase,zero.energy_scale_K,1); verifyTrue(testCase,zero.zero_fallback_used);
weak=ResolveThermalEnergyScale(-1e-12,[0,0,0],[2,2,2]);
verifyEqual(testCase,weak.energy_scale_K,1e-12); verifyFalse(testCase,weak.zero_fallback_used);
end

function testRejectNonorthogonalAxes(testCase)
verifyError(testCase,@() ResolveThermalEnergyScale(1,[0,0,0],[2,2,2],ones(3)), ...
    'qspintn:EnergyScaleAxes');
end
