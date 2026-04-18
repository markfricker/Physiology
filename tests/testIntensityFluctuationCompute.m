classdef testIntensityFluctuationCompute < matlab.unittest.TestCase
%TESTINTENSITYFLUCTUATIONCOMPUTE  Unit tests for intensityFluctuationCompute.
%
%   Run with:
%       results = runtests('testIntensityFluctuationCompute');
%       disp(results)

    % ------------------------------------------------------------------
    % Shared helpers
    % ------------------------------------------------------------------
    methods (Static)

        function p = defaultParams(method)
            %DEFAULTPARAMS  Minimal param struct with sensible defaults.
            if nargin < 1, method = 'dff'; end
            p.channelIdx = 1;
            p.method     = method;
            p.baseline   = struct('frames', [1 3]);
        end

        function imIn = constantSignal(nY, nX, nT, level)
            %CONSTANTSIGNAL  Flat signal at a fixed intensity level.
            if nargin < 4, level = 0.5; end
            imIn = repmat(single(level), [nY nX 1 1 nT]);
        end

        function imIn = stepSignal(nY, nX, nT, baseline, stepFrame, stepLevel)
            %STEPSIGNAL  Baseline then step increase from stepFrame onward.
            imIn = repmat(single(baseline), [nY nX 1 1 nT]);
            imIn(:,:,1,1,stepFrame:end) = stepLevel;
        end

    end

    % ------------------------------------------------------------------
    % Output size and type
    % ------------------------------------------------------------------
    methods (Test)

        function testOutputSizeSeries(tc)
            nY = 16; nX = 16; nC = 2; nZ = 2; nT = 8;
            imIn = rand(nY, nX, nC, nZ, nT, 'single');
            p = tc.defaultParams();
            [series, mn, bl] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(size(series), [nY nX nC nZ nT]);
            tc.verifyEqual(size(mn,1), nY);  tc.verifyEqual(size(mn,2), nX);
            tc.verifyEqual(size(mn,3), nC);  tc.verifyEqual(size(mn,4), nZ);
            tc.verifyEqual(size(bl,1), nY);  tc.verifyEqual(size(bl,2), nX);
        end

        function testOutputIsSingle(tc)
            imIn = rand(16, 16, 1, 1, 6, 'single');
            p = tc.defaultParams();
            [series, mn, bl] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(class(series), 'single');
            tc.verifyEqual(class(mn),     'single');
            tc.verifyEqual(class(bl),     'single');
        end

        function testUnselectedChannelsAreZero(tc)
            % channelIdx=1 → channel 2 output should remain zero
            imIn = rand(16, 16, 2, 1, 6, 'single');
            p = tc.defaultParams(); p.channelIdx = 1;
            [series, mn, bl] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(series(:,:,2,:,:), zeros(16,16,1,1,6,'single'), ...
                'Unselected channel should be zero in series.');
            tc.verifyEqual(mn(:,:,2,:), zeros(16,16,1,1,'single'));
            tc.verifyEqual(bl(:,:,2,:), zeros(16,16,1,1,'single'));
        end

    end

    % ------------------------------------------------------------------
    % dff method: known ground truth
    % ------------------------------------------------------------------
    methods (Test)

        function testDffConstantSignalIsZero(tc)
            % F(t) = F0 at all t → ΔF/F₀ = 0 everywhere
            imIn = tc.constantSignal(16, 16, 8, 0.5);
            p = tc.defaultParams('dff');
            [series, mn, ~] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(double(mean(abs(series(:)))), 0, 'AbsTol', 1e-5, ...
                'Constant signal: dF/F should be zero.');
            tc.verifyEqual(double(mean(abs(mn(:)))), 0, 'AbsTol', 1e-5);
        end

        function testDffStepSignalCorrectMagnitude(tc)
            % Baseline = 0.4 (frames 1-3), step to 0.6 (frames 4-8)
            % ΔF/F₀ at step frames = (0.6 - 0.4) / 0.4 = 0.5
            imIn = tc.stepSignal(16, 16, 8, 0.4, 4, 0.6);
            p = tc.defaultParams('dff');
            p.baseline.frames = [1 3];
            [series, ~, ~] = intensityFluctuationCompute(imIn, p);
            dffStep = double(series(1,1,1,1,4));
            tc.verifyEqual(dffStep, 0.5, 'AbsTol', 1e-4, ...
                'Step of 50%% above baseline should give dF/F = 0.5.');
        end

        function testDffBaselineFramesAreNearZero(tc)
            % Baseline frames should have ΔF/F₀ ≈ 0 (by definition F=F0)
            imIn = tc.stepSignal(16, 16, 8, 0.4, 6, 0.8);
            p = tc.defaultParams('dff');
            p.baseline.frames = [1 4];
            [series, ~, ~] = intensityFluctuationCompute(imIn, p);
            baselineDff = series(:,:,:,:,1:4);
            tc.verifyEqual(double(mean(abs(baselineDff(:)))), 0, 'AbsTol', 1e-5, ...
                'Baseline frames should have zero dF/F.');
        end

        function testDffMeanEqualsTemporalMeanOfSeries(tc)
            imIn = rand(16, 16, 1, 1, 10, 'single');
            p = tc.defaultParams('dff');
            [series, mn, ~] = intensityFluctuationCompute(imIn, p);
            expected = mean(series, 5);
            tc.verifyEqual(mn, expected, 'AbsTol', 1e-5, ...
                'fluctuationMean should equal temporal mean of series.');
        end

    end

    % ------------------------------------------------------------------
    % variance method
    % ------------------------------------------------------------------
    methods (Test)

        function testVarianceConstantSignalIsZero(tc)
            imIn = tc.constantSignal(16, 16, 8, 0.5);
            p = tc.defaultParams('variance');
            [~, mn, ~] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(double(max(abs(mn(:)))), 0, 'AbsTol', 1e-5, ...
                'Constant signal: variance should be zero.');
        end

        function testVarianceKnownValue(tc)
            % Signal alternates 0.2 / 0.8 every frame → variance = 0.09
            nT = 10;
            imIn = zeros(4, 4, 1, 1, nT, 'single');
            for iT = 1:nT
                if mod(iT,2)==0, imIn(:,:,1,1,iT) = 0.8;
                else,            imIn(:,:,1,1,iT) = 0.2; end
            end
            p = tc.defaultParams('variance');
            [~, mn, ~] = intensityFluctuationCompute(imIn, p);
            % var([0.2 0.8 0.2 0.8 ...]) = 0.09
            tc.verifyEqual(double(mn(1,1,1,1)), 0.09, 'AbsTol', 0.002, ...
                'Alternating signal: variance should be 0.09.');
        end

    end

    % ------------------------------------------------------------------
    % cv method
    % ------------------------------------------------------------------
    methods (Test)

        function testCvConstantSignalIsZero(tc)
            imIn = tc.constantSignal(16, 16, 8, 0.6);
            p = tc.defaultParams('cv');
            [~, mn, ~] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(double(max(abs(mn(:)))), 0, 'AbsTol', 1e-5, ...
                'Constant signal: CV should be zero.');
        end

        function testCvIsNonnegative(tc)
            imIn = rand(16, 16, 1, 1, 8, 'single') * 0.5 + 0.1;
            p = tc.defaultParams('cv');
            [~, mn, ~] = intensityFluctuationCompute(imIn, p);
            tc.verifyGreaterThanOrEqual(double(min(mn(:))), 0, ...
                'CV should be non-negative.');
        end

    end

    % ------------------------------------------------------------------
    % Background / zero-intensity pixels
    % ------------------------------------------------------------------
    methods (Test)

        function testZeroPixelsProduceNoInfOrNaN(tc)
            imIn = rand(16, 16, 1, 1, 8, 'single') * 0.5;
            imIn(1:4, 1:4, :, :, :) = 0;   % background patch
            p = tc.defaultParams('dff');
            [series, mn, ~] = intensityFluctuationCompute(imIn, p);
            tc.verifyFalse(any(isinf(series(:))), 'Series should have no Inf.');
            tc.verifyFalse(any(isinf(mn(:))),     'Mean should have no Inf.');
        end

        function testZeroPixelSeriesIsZeroNotNaN(tc)
            % Background pixels (F0≈0) should be zeroed out, not left as NaN
            imIn = zeros(8, 8, 1, 1, 6, 'single');
            imIn(4:6, 4:6, :, :, :) = 0.5;   % only inner patch has signal
            p = tc.defaultParams('dff');
            [series, ~, ~] = intensityFluctuationCompute(imIn, p);
            bgSeries = series(1:3, 1:3, :, :, :);
            tc.verifyFalse(any(isnan(bgSeries(:))), ...
                'Background pixels should be 0, not NaN.');
            tc.verifyEqual(double(max(abs(bgSeries(:)))), 0, 'AbsTol', 1e-6);
        end

    end

    % ------------------------------------------------------------------
    % Error handling
    % ------------------------------------------------------------------
    methods (Test)

        function testUnknownMethodErrors(tc)
            imIn = rand(8, 8, 1, 1, 4, 'single');
            p = tc.defaultParams('bogus');
            tc.verifyError(@() intensityFluctuationCompute(imIn, p), ...
                'intensityFluctuationCompute:unknownMethod');
        end

    end

    % ------------------------------------------------------------------
    % Multi-Z
    % ------------------------------------------------------------------
    methods (Test)

        function testMultiZRunsAndOutputsFilledForAllZ(tc)
            imIn = rand(12, 12, 1, 3, 8, 'single') * 0.5 + 0.1;
            p = tc.defaultParams('dff');
            [series, mn, bl] = intensityFluctuationCompute(imIn, p);
            tc.verifyEqual(size(series,4), 3);
            % Each Z plane should be non-zero in the series
            for iZ = 1:3
                tc.verifyFalse(all(series(:,:,1,iZ,:) == 0, 'all'), ...
                    sprintf('Z=%d series should be non-zero.', iZ));
            end
        end

    end

end
