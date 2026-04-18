classdef testFrapCompute < matlab.unittest.TestCase
%TESTFRAPCOMPUTE  Unit tests for frapCompute.
%
%   Run with:
%       results = runtests('testFrapCompute');
%       disp(results)

    % ------------------------------------------------------------------
    % Shared helpers
    % ------------------------------------------------------------------
    methods (Static)

        function p = defaultParams()
            p.channelIdx  = 1;
            p.method      = 'frap';
            p.bleachFrame = 4;   % explicit bleach at frame 4
            p.preFrames   = 0;   % use all pre-bleach frames
        end

        function [imIn, skIn, cistIn] = makeInputs(nY, nX, nT, ...
                                                    bleachFrame, recoveryTau)
            %MAKEINPUTS  Synthetic FRAP image series.
            % Pre-bleach: uniform intensity 0.8.
            % Bleach frame: intensity inside circular ROI drops to 0.1.
            % Recovery: exponential toward 0.8 with time constant recoveryTau.
            % skIn is a thin horizontal band through the bleach ROI.
            if nargin < 4, bleachFrame  = 4;  end
            if nargin < 5, recoveryTau  = 3;  end

            imIn  = zeros(nY, nX, 1, 1, nT, 'single');
            skIn  = false(nY, nX, 1, 1, nT);
            cistIn = false(nY, nX, 1, 1, nT);

            % ROI: small square in centre
            cy = round(nY/2); cx = round(nX/2); r = max(2, round(min(nY,nX)/6));
            roi = false(nY, nX);
            for iy = max(1,cy-r):min(nY,cy+r)
                for ix = max(1,cx-r):min(nX,cx+r)
                    if (iy-cy)^2 + (ix-cx)^2 <= r^2
                        roi(iy,ix) = true;
                    end
                end
            end

            Fpre    = 0.8;
            Fbleach = 0.1;

            for iT = 1:nT
                frame = repmat(Fpre, nY, nX);
                if iT >= bleachFrame
                    t  = iT - bleachFrame;
                    Ft = Fbleach + (Fpre - Fbleach) * (1 - exp(-t / recoveryTau));
                    frame(roi) = Ft;
                end
                imIn(:,:,1,1,iT) = single(frame);
            end

            % Skeleton: horizontal midline through ROI
            skIn(cy, :, 1, 1, :) = true;
        end

    end

    % ------------------------------------------------------------------
    % Output size and type
    % ------------------------------------------------------------------
    methods (Test)

        function testOutputSizeRecoverySeries(tc)
            nY = 20; nX = 20; nC = 1; nZ = 1; nT = 10;
            [imIn, skIn, cistIn] = tc.makeInputs(nY, nX, nT);
            p = tc.defaultParams();
            [rec, rm, pre, ~] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyEqual(size(rec), [nY nX nC nZ nT]);
            tc.verifyEqual(size(rm,1), nY); tc.verifyEqual(size(rm,2), nX);
            tc.verifyEqual(size(rm,3), nC); tc.verifyEqual(size(rm,4), nZ);
            tc.verifyEqual(size(pre,1), nY); tc.verifyEqual(size(pre,2), nX);
        end

        function testOutputIsSingle(tc)
            [imIn, skIn, cistIn] = tc.makeInputs(16, 16, 8);
            p = tc.defaultParams();
            [rec, ~, ~, ~] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyEqual(class(rec), 'single');
        end

    end

    % ------------------------------------------------------------------
    % fitResults struct
    % ------------------------------------------------------------------
    methods (Test)

        function testFitResultsHasRequiredFields(tc)
            [imIn, skIn, cistIn] = tc.makeInputs(16, 16, 8);
            p = tc.defaultParams();
            [~, ~, ~, fr] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyTrue(isfield(fr, 'bleachFrame'),   'fitResults missing bleachFrame.');
            tc.verifyTrue(isfield(fr, 'method'),        'fitResults missing method.');
            tc.verifyTrue(isfield(fr, 'skeletonCurve'), 'fitResults missing skeletonCurve.');
            tc.verifyTrue(isfield(fr, 'cisternaCurve'), 'fitResults missing cisternaCurve.');
        end

        function testFitResultsMethodMatchesParam(tc)
            [imIn, skIn, cistIn] = tc.makeInputs(16, 16, 8);
            p = tc.defaultParams(); p.method = 'flip';
            [~, ~, ~, fr] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyEqual(fr.method, 'flip');
        end

        function testFitResultsBleachFrameMatchesParam(tc)
            [imIn, skIn, cistIn] = tc.makeInputs(16, 16, 10, 5);
            p = tc.defaultParams(); p.bleachFrame = 5;
            [~, ~, ~, fr] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyEqual(fr.bleachFrame, 5);
        end

        function testSkeletonCurveSize(tc)
            nT = 10;
            [imIn, skIn, cistIn] = tc.makeInputs(16, 16, nT);
            p = tc.defaultParams();
            [~, ~, ~, fr] = frapCompute(imIn, skIn, cistIn, p);
            % skeletonCurve: [nCh × nZ × nT]
            tc.verifyEqual(size(fr.skeletonCurve, 3), nT, ...
                'skeletonCurve third dimension should be nT.');
        end

    end

    % ------------------------------------------------------------------
    % Normalisation: pre-bleach frames = 1
    % ------------------------------------------------------------------
    methods (Test)

        function testPrebleachFramesAreOne(tc)
            nT = 10; bleachFrame = 4;
            [imIn, skIn, cistIn] = tc.makeInputs(20, 20, nT, bleachFrame, 5);
            p = tc.defaultParams(); p.bleachFrame = bleachFrame;
            [rec, ~, ~, ~] = frapCompute(imIn, skIn, cistIn, p);
            % Interior pixels (outside bleach ROI) in pre-bleach frames = 1
            preRec = rec(1, 1, 1, 1, 1:bleachFrame-1);
            tc.verifyEqual(double(preRec(:)), ones(bleachFrame-1, 1), ...
                'AbsTol', 1e-5, ...
                'Pre-bleach frames should normalise to 1.');
        end

        function testBleachFrameDropsBelowOne(tc)
            % The bleach frame itself should show a drop in the bleach ROI
            nT = 10; bleachFrame = 4;
            [imIn, skIn, cistIn] = tc.makeInputs(20, 20, nT, bleachFrame, 5);
            p = tc.defaultParams(); p.bleachFrame = bleachFrame;
            [rec, ~, ~, ~] = frapCompute(imIn, skIn, cistIn, p);
            % Centre pixel is in bleach ROI — should have rec < 1 at bleachFrame
            cy = 10; cx = 10;
            recCentre = double(rec(cy, cx, 1, 1, bleachFrame));
            tc.verifyLessThan(recCentre, 0.5, ...
                'Centre pixel at bleach frame should show strong drop below 1.');
        end

    end

    % ------------------------------------------------------------------
    % Bleach auto-detection
    % ------------------------------------------------------------------
    methods (Test)

        function testAutoDetectBleachFrame(tc)
            % Drop from 0.8 to 0.1 between frame 4 and 5 → detected bleach = 5
            nT = 10; trueBleach = 5;
            [imIn, skIn, cistIn] = tc.makeInputs(20, 20, nT, trueBleach, 5);
            p = tc.defaultParams(); p.bleachFrame = 0;   % auto-detect
            [~, ~, ~, fr] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyEqual(fr.bleachFrame, trueBleach, ...
                'Auto-detected bleach frame should match true bleach frame.');
        end

    end

    % ------------------------------------------------------------------
    % prebleachMap
    % ------------------------------------------------------------------
    methods (Test)

        function testPrebleachMapMatchesPreBleachMean(tc)
            nT = 10; bleachFrame = 4;
            [imIn, skIn, cistIn] = tc.makeInputs(20, 20, nT, bleachFrame, 4);
            p = tc.defaultParams(); p.bleachFrame = bleachFrame;
            [~, ~, pre, ~] = frapCompute(imIn, skIn, cistIn, p);
            % Outside bleach ROI, pre-bleach mean ≈ 0.8
            tc.verifyEqual(double(pre(1,1,1,1)), 0.8, 'AbsTol', 0.01, ...
                'prebleachMap should equal mean pre-bleach intensity.');
        end

    end

    % ------------------------------------------------------------------
    % No Inf or NaN in normal operating conditions
    % ------------------------------------------------------------------
    methods (Test)

        function testNoInfOrNaNInRecoverySeries(tc)
            [imIn, skIn, cistIn] = tc.makeInputs(20, 20, 10, 4, 3);
            p = tc.defaultParams();
            [rec, ~, ~, ~] = frapCompute(imIn, skIn, cistIn, p);
            tc.verifyFalse(any(isinf(rec(:))), 'recoverySeries should have no Inf.');
        end

        function testEmptyMasksRunWithoutError(tc)
            % All-false skeleton and cisternae masks should not cause errors
            [imIn, ~, ~] = tc.makeInputs(16, 16, 8, 4, 3);
            skIn   = false(16, 16, 1, 1, 8);
            cistIn = false(16, 16, 1, 1, 8);
            p = tc.defaultParams();
            tc.verifyWarningFree(@() frapCompute(imIn, skIn, cistIn, p));
        end

    end

end
