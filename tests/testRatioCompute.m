classdef testRatioCompute < matlab.unittest.TestCase
%TESTRATIOCOMPUTE  Unit tests for ratioCompute.
%
%   Run with:
%       results = runtests('testRatioCompute');
%       disp(results)

    % ------------------------------------------------------------------
    % Shared helpers
    % ------------------------------------------------------------------
    methods (Static)

        function p = defaultParams()
            %DEFAULTPARAMS  Two-channel ratio, no AF, no calibration.
            p.numeratorChannel   = 1;
            p.denominatorChannel = 2;
        end

        function imIn = twoChannelImage(nY, nX, nT, numLevel, denLevel)
            %TWOCHANNELIMAGE  Spatially uniform two-channel image.
            imIn = zeros(nY, nX, 2, 1, nT, 'single');
            imIn(:,:,1,:,:) = numLevel;
            imIn(:,:,2,:,:) = denLevel;
        end

        function imIn = threeChannelImage(nY, nX, nT, numLevel, denLevel, afLevel)
            %THREECHANNELIMAGE  Three-channel image with AF channel.
            imIn = zeros(nY, nX, 3, 1, nT, 'single');
            imIn(:,:,1,:,:) = numLevel;
            imIn(:,:,2,:,:) = denLevel;
            imIn(:,:,3,:,:) = afLevel;
        end

    end

    % ------------------------------------------------------------------
    % Output size and type
    % ------------------------------------------------------------------
    methods (Test)

        function testOutputSizeIsNYNX1NZnT(tc)
            nY = 16; nX = 16; nZ = 2; nT = 5;
            imIn = rand(nY, nX, 2, nZ, nT, 'single');
            p = tc.defaultParams();
            [rs, rm, cs, cm] = ratioCompute(imIn, p);
            tc.verifyEqual(size(rs), [nY nX 1 nZ nT]);
            tc.verifyEqual(size(rm,1), nY);  tc.verifyEqual(size(rm,2), nX);
            tc.verifyEqual(size(rm,3), 1);   tc.verifyEqual(size(rm,4), nZ);
            tc.verifyEqual(size(cs), [nY nX 1 nZ nT]);
        end

        function testOutputIsSingle(tc)
            imIn = rand(8, 8, 2, 1, 4, 'single');
            p = tc.defaultParams();
            [rs,~,~,~] = ratioCompute(imIn, p);
            tc.verifyEqual(class(rs), 'single');
        end

    end

    % ------------------------------------------------------------------
    % Basic ratio arithmetic
    % ------------------------------------------------------------------
    methods (Test)

        function testSimpleRatioCorrectValue(tc)
            % num=0.6, den=0.3 → ratio = 2.0
            imIn = tc.twoChannelImage(8, 8, 4, 0.6, 0.3);
            p = tc.defaultParams();
            [rs, rm, ~, ~] = ratioCompute(imIn, p);
            tc.verifyEqual(double(rs(1,1,1,1,1)), 2.0, 'AbsTol', 1e-4, ...
                'ratio should be num/den = 2.0.');
            tc.verifyEqual(double(rm(1,1,1,1)), 2.0, 'AbsTol', 1e-4, ...
                'ratioMean should equal ratio for constant input.');
        end

        function testRatioMeanEqualsTemporalMean(tc)
            % ratioMean should equal the mean of ratioSeries over T
            imIn = rand(12, 12, 2, 1, 8, 'single') * 0.4 + 0.1;
            imIn(:,:,2,:,:) = imIn(:,:,2,:,:) + 0.3;   % keep den > num
            p = tc.defaultParams();
            [rs, rm, ~, ~] = ratioCompute(imIn, p);
            expected = mean(rs, 5, 'omitnan');
            tc.verifyEqual(rm, expected, 'AbsTol', 1e-5, ...
                'ratioMean should equal temporal mean of ratioSeries.');
        end

        function testRatioUnityWhenChannelsEqual(tc)
            % num == den → ratio = 1 everywhere
            imIn = tc.twoChannelImage(8, 8, 4, 0.5, 0.5);
            p = tc.defaultParams();
            [rs, ~, ~, ~] = ratioCompute(imIn, p);
            tc.verifyEqual(double(mean(rs(:))), 1.0, 'AbsTol', 1e-4);
        end

    end

    % ------------------------------------------------------------------
    % Autofluorescence correction
    % ------------------------------------------------------------------
    methods (Test)

        function testAfCorrectionDenominatorOnly(tc)
            % num=0.6, den=0.6, AF=0.3, scaleD=1, scaleN=0
            % → corrected: num=0.6, den=0.6-0.3=0.3 → ratio=2.0
            imIn = tc.threeChannelImage(8, 8, 4, 0.6, 0.6, 0.3);
            p = tc.defaultParams();
            p.autofluorescence.use               = true;
            p.autofluorescence.channel           = 3;
            p.autofluorescence.scaleDenominator  = 1.0;
            p.autofluorescence.scaleNumerator    = 0.0;
            [rs, ~, ~, ~] = ratioCompute(imIn, p);
            tc.verifyEqual(double(rs(1,1,1,1,1)), 2.0, 'AbsTol', 1e-3, ...
                'Denominator-only AF correction: ratio should be 2.0.');
        end

        function testAfCorrectionBothChannels(tc)
            % num=0.8, den=0.6, AF=0.2, scaleN=1, scaleD=1
            % → corrected: num=0.6, den=0.4 → ratio=1.5
            imIn = tc.threeChannelImage(8, 8, 4, 0.8, 0.6, 0.2);
            p = tc.defaultParams();
            p.autofluorescence.use               = true;
            p.autofluorescence.channel           = 3;
            p.autofluorescence.scaleDenominator  = 1.0;
            p.autofluorescence.scaleNumerator    = 1.0;
            [rs, ~, ~, ~] = ratioCompute(imIn, p);
            tc.verifyEqual(double(rs(1,1,1,1,1)), 1.5, 'AbsTol', 1e-3, ...
                'Both-channel AF correction: ratio should be 1.5.');
        end

        function testAfOvercorrectionGivesNaN(tc)
            % AF > denominator → corrected den ≤ 0 → ratio should be NaN
            imIn = tc.threeChannelImage(8, 8, 3, 0.5, 0.3, 0.5);
            p = tc.defaultParams();
            p.autofluorescence.use              = true;
            p.autofluorescence.channel          = 3;
            p.autofluorescence.scaleDenominator = 1.0;
            p.autofluorescence.scaleNumerator   = 0.0;
            [rs, ~, ~, ~] = ratioCompute(imIn, p);
            tc.verifyTrue(all(isnan(rs(:))), ...
                'Over-corrected denominator should produce NaN ratio.');
        end

    end

    % ------------------------------------------------------------------
    % Calibration: none
    % ------------------------------------------------------------------
    methods (Test)

        function testNoCalibrationGivesNaNCalibrated(tc)
            imIn = tc.twoChannelImage(8, 8, 4, 0.5, 0.25);
            p = tc.defaultParams();
            [~, ~, cs, cm] = ratioCompute(imIn, p);
            tc.verifyTrue(all(isnan(cs(:))), ...
                'calibratedSeries should be all-NaN when method=''none''.');
            tc.verifyTrue(all(isnan(cm(:))));
        end

    end

    % ------------------------------------------------------------------
    % Calibration: linear
    % ------------------------------------------------------------------
    methods (Test)

        function testLinearCalibration(tc)
            % ratio=2.0, slope=3, intercept=-1 → calibrated=5.0
            imIn = tc.twoChannelImage(8, 8, 4, 0.6, 0.3);
            p = tc.defaultParams();
            p.calibration.method    = 'linear';
            p.calibration.slope     = 3;
            p.calibration.intercept = -1;
            [~, ~, cs, cm] = ratioCompute(imIn, p);
            tc.verifyEqual(double(cs(1,1,1,1,1)), 5.0, 'AbsTol', 1e-4, ...
                'Linear calibration: 3*2.0 - 1 should give 5.0.');
            tc.verifyEqual(double(cm(1,1,1,1)), 5.0, 'AbsTol', 1e-4);
        end

    end

    % ------------------------------------------------------------------
    % Calibration: Henderson-Hasselbalch
    % ------------------------------------------------------------------
    methods (Test)

        function testHHCalibrationAtMidpointGivesPKa(tc)
            % R = (Racid + Rbase)/2 → pH = pKa + log10(1) = pKa
            Racid = 0.2; Rbase = 1.0; pKa = 7.0;
            Rmid  = (Racid + Rbase) / 2;    % 0.6
            imIn  = tc.twoChannelImage(8, 8, 4, Rmid, 1.0);
            % ratio = Rmid / 1.0 = Rmid
            p = tc.defaultParams();
            p.calibration.method = 'henderson-hasselbalch';
            p.calibration.pKa    = pKa;
            p.calibration.Racid  = Racid;
            p.calibration.Rbase  = Rbase;
            [~, ~, cs, ~] = ratioCompute(imIn, p);
            tc.verifyEqual(double(cs(1,1,1,1,1)), pKa, 'AbsTol', 0.01, ...
                'HH: at midpoint ratio, calibrated value should equal pKa.');
        end

        function testHHOutsideRangeIsNaN(tc)
            % R < Racid → pH is undefined
            Racid = 0.3; Rbase = 0.9;
            imIn  = tc.twoChannelImage(8, 8, 3, 0.1, 1.0);  % ratio=0.1 < Racid
            p = tc.defaultParams();
            p.calibration.method = 'henderson-hasselbalch';
            p.calibration.pKa    = 7.0;
            p.calibration.Racid  = Racid;
            p.calibration.Rbase  = Rbase;
            [~, ~, cs, ~] = ratioCompute(imIn, p);
            tc.verifyTrue(all(isnan(cs(:))), ...
                'Ratio below Racid should give NaN calibrated value.');
        end

        function testHHMissingParamsErrors(tc)
            imIn = tc.twoChannelImage(8, 8, 3, 0.5, 1.0);
            p = tc.defaultParams();
            p.calibration.method = 'henderson-hasselbalch';
            p.calibration.pKa    = 7.0;
            % Missing Racid / Rbase → error
            tc.verifyError(@() ratioCompute(imIn, p), ...
                'ratioCompute:missingCalibration');
        end

    end

    % ------------------------------------------------------------------
    % Multi-Z
    % ------------------------------------------------------------------
    methods (Test)

        function testMultiZOutputSize(tc)
            imIn = rand(10, 10, 2, 3, 5, 'single') * 0.4 + 0.1;
            p = tc.defaultParams();
            [rs, rm, ~, ~] = ratioCompute(imIn, p);
            tc.verifyEqual(size(rs,4), 3, 'nZ=3 should be preserved in output.');
            tc.verifyEqual(size(rm,4), 3);
        end

    end

end
