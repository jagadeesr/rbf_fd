% UT_find_stencil_weights unit test for find_stencil_weights
classdef ut_find_stencil_weights < matlab.unittest.TestCase
    
    methods(Test)
        % ut_find_stencil_weights tests find_stencil_weights

        function test1(testCase)
            tol = 1e-4;
			epsilon = 0.0146;
			RBFQR_flag = 1;
			stencil_support = stencil_support_selection(dm1, p1, 25);
			[o1, o2, o3] = find_stencil_weights(p1(:,stencil_support)', epsilon, RBFQR_flag);
			weights = [-39.3643   10.0898    7.0783    8.1554    6.0287    7.0845    0.9276];
			v = 2.3017e-016;
			stable_flag = 1;
            error = norm(weights-o1);
            testCase.verifyLessThanOrEqual(error,tol);
			testCase.verifyLessThanOrEqual(abs(v-o2),tol);
			testCase.verifyEqual(o3,stable_flag);
        end

        function test2(testCase)
            tol = 1e-4;
			epsilon = 0.02386;
			RBFQR_flag = 1;
			stencil_support = stencil_support_selection(dm2, p2, 45);
			[o1, o2, o3] = find_stencil_weights(p2(:,stencil_support)', epsilon, RBFQR_flag);
			weights = [-184.3131   29.6522   27.1442   31.1328   35.8469   34.8143   25.7227];
			v = 2.3057e-015;
			stable_flag = 1;
            error = norm(weights-o1);
            testCase.verifyLessThanOrEqual(error,tol);
			testCase.verifyLessThanOrEqual(abs(v-o2),tol);
			testCase.verifyEqual(o3,stable_flag);
        end
        
    end
end
