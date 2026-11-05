function tests = test_svdT
% Cross-release regression tests for epsilon/Nkeep singular-value truncation.
% Run from an external working directory with this tests directory on the path:
%   results = runtests('test_svdT'); assertSuccess(results);
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.original_path = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'ManyBodySolver', 'TensorFunction'), '-begin');
testCase.assertEqual(which('svdT'), fullfile(root, 'ManyBodySolver', 'TensorFunction', 'svdT.m'));
end

function teardownOnce(testCase)
path(testCase.TestData.original_path);
end

function testMultipleDiscardedValues(testCase)
% The old find(...)-1 yields [2;3], invalid as a colon endpoint in R2025a+.
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'Nkeep', 4, 'epsilon', 1e-8);
verifyTruncation(testCase, U, s, V, err, entropy, values, 2);
end

function testEpsilonBeforeNkeep(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'epsilon', 1e-8, 'Nkeep', 4);
verifyTruncation(testCase, U, s, V, err, entropy, values, 2);
end

function testNkeepIsMoreRestrictive(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'Nkeep', 1, 'epsilon', 1e-8);
verifyTruncation(testCase, U, s, V, err, entropy, values, 1);
end

function testEpsilonOnly(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'epsilon', 1e-8);
verifyTruncation(testCase, U, s, V, err, entropy, values, 2);
end

function testNoValuesBelowCutoff(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'epsilon', 1e-15);
verifyTruncation(testCase, U, s, V, err, entropy, values, 4);
end

function testSingleDiscardedValue(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'epsilon', 1e-11);
verifyTruncation(testCase, U, s, V, err, entropy, values, 3);
end

function testAllValuesBelowCutoff(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'epsilon', 2);
verifyTruncation(testCase, U, s, V, err, entropy, values, 0);
end

function testStrictThresholdEquality(testCase)
% A scalar singular value normalized to one is retained at epsilon=1.
[U,s,V,err,entropy] = svdT(2, 'epsilon', 1);
verifyTruncation(testCase, U, s, V, err, entropy, 2, 1);
end

function testNkeepOnly(testCase)
values = [1;.1;1e-10;1e-12];
[U,s,V,err,entropy] = svdT(diag(values), 'Nkeep', 3);
verifyTruncation(testCase, U, s, V, err, entropy, values, 3);
end

function testRectangularComplexAndSingleReconstruction(testCase)
inputs = {[1,2;3,5;7,11], [1,2,3;5,7,11], [1+1i,2;3i,4-2i;5,6], ...
          [1,2,3], [1;2;3], single([1,2;3,4;5,6])};
for i=1:numel(inputs)
    A = inputs{i};
    [U,s,V] = svdT(A);
    testCase.verifySize(s, [min(size(A)),1]);
    if isa(A,'single'), tolerance=2e-6; else, tolerance=2e-14; end
    % svdT intentionally returns V already conjugate-transposed.
    testCase.verifyLessThan(norm(A-U*diag(s)*V,'fro')/norm(A,'fro'), tolerance);
end
end

function verifyTruncation(testCase,U,s,V,err,entropy,values,keep)
values = values(:);
testCase.verifyEqual(numel(s),keep);
testCase.verifySize(U,[numel(values),keep]);
testCase.verifySize(V,[keep,numel(values)]);
testCase.verifyEqual(s(:),values(1:keep),'AbsTol',1e-15);
expected = diag([values(1:keep); zeros(numel(values)-keep,1)]);
testCase.verifyEqual(U*diag(s)*V,expected,'AbsTol',1e-15);
weights = values.^2 / sum(values.^2);
testCase.verifyEqual(err,sum(weights(keep+1:end)),'AbsTol',1e-30,'RelTol',1e-13);
testCase.verifyEqual(entropy,-sum(weights.*log(weights)),'AbsTol',1e-14);
end
