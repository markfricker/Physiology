classdef testRedoxCalibrate < matlab.unittest.TestCase
%TESTREDOXCALIBRATE  Unit tests for redoxCalibrate.
%
%   Run with:
%       results = runtests('testRedoxCalibrate');
%       disp(results)

    % ------------------------------------------------------------------
    % Shared helpers
    % ------------------------------------------------------------------
    methods (Static)

        function p = defaultParams()
            %DEFAULTPARAMS  Typical grx1-roGFP2 calibration parameters.
            p.rMin     = 0.1;    % fully-reduced endpoint ratio
            p.rMax     = 1.0;    % fully-oxidised endpoint ratio
            p.iFactor  = 5.5;    % instrument factor (typical for roGFP2)
            p.midpoint = -280;   % E°' for roGFP2 vs SHE at reference pH (mV)
            p.temp     = 293;    % 20 °C in K
            p.pH       = 7.0;
            p.pH0      = 7.0;    % same → no pH correction
            p.n        = 2;
        end

        function ratioIn = uniformRatio(nY, nX, nT, value)
            %UNIFORMRATIO  Spatially uniform ratio image at a single value.
            ratioIn = repmat(single(value), [nY nX 1 1 nT]);
        end

        function ratioIn = gradientRatio(nY, nX, nT, rMin, rMax)
            %GRADIENTRATIO  Horizontal gradient from rMin (left) to rMax (right).
            x = linspace(rMin, rMax, nX);
            ratioIn = repmat(reshape(single(x),[1 nX 1 1 1]), [nY 1 1 1 nT]);
        end

    end

    % ------------------------------------------------------------------
    % Output size and type
    % ------------------------------------------------------------------
    methods (Test)

        function testOutputSizeMatchesInput(tc)
            nY = 16; nX = 16; nC = 1; nZ = 2; nT = 5;
            ratioIn = rand(nY, nX, nC, nZ, nT, 'single') * 0.9 + 0.1;
            p = tc.defaultParams();
            [oxd, eh] = redoxCalibrate(ratioIn, p);
            tc.verifyEqual(size(oxd), [nY nX nC nZ nT]);
            tc.verifyEqual(size(eh),  [nY nX nC nZ nT]);
        end

        function testOutputIsSingle(tc)
            ratioIn = tc.uniformRatio(8, 8, 4, 0.5);
            p = tc.defaultParams();
            [oxd, eh] = redoxCalibrate(ratioIn, p);
            tc.verifyEqual(class(oxd), 'single');
            tc.verifyEqual(class(eh),  'single');
        end

    end

    % ------------------------------------------------------------------
    % OxD: range, monotonicity, endpoints
    % ------------------------------------------------------------------
    methods (Test)

        function testOxdBoundedBetweenZeroAndOne(tc)
            ratioIn = tc.gradientRatio(4, 32, 3, 0.05, 1.1);  % deliberately outside [rMin rMax]
            p = tc.defaultParams();
            [oxd, ~] = redoxCalibrate(ratioIn, p);
            % Compare against the single-precision representations of the
            % clamp bounds -- oxd is single, and single(0.999) round-trips
            % to a value fractionally above the double literal 0.999.
            tc.verifyGreaterThanOrEqual(double(min(oxd(:))), double(single(0.001)), ...
                'OxD should be ≥ 0.001 (clamped).');
            tc.verifyLessThanOrEqual(double(max(oxd(:))), double(single(0.999)), ...
                'OxD should be ≤ 0.999 (clamped).');
        end

        function testOxdAtRminIsNearZero(tc)
            % R = rMin → fully reduced → OxD ≈ 0.001
            ratioIn = tc.uniformRatio(4, 4, 3, 0.1);   % rMin = 0.1
            p = tc.defaultParams();
            [oxd, ~] = redoxCalibrate(ratioIn, p);
            tc.verifyEqual(double(oxd(1,1,1,1,1)), 0.001, 'AbsTol', 1e-4, ...
                'OxD at rMin should be clamped to 0.001.');
        end

        function testOxdAtRmaxIsNearOne(tc)
            % R = rMax → fully oxidised → OxD ≈ 0.999
            ratioIn = tc.uniformRatio(4, 4, 3, 1.0);   % rMax = 1.0
            p = tc.defaultParams();
            [oxd, ~] = redoxCalibrate(ratioIn, p);
            tc.verifyEqual(double(oxd(1,1,1,1,1)), 0.999, 'AbsTol', 1e-4, ...
                'OxD at rMax should be clamped to 0.999.');
        end

        function testOxdMonotonicallyIncreasing(tc)
            % As R increases from rMin to rMax, OxD should increase
            p = tc.defaultParams();
            Rs  = single(linspace(p.rMin, p.rMax, 20));
            oxds = zeros(1, 20, 'single');
            for k = 1:20
                r = repmat(Rs(k), [2 2 1 1 2]);
                [oxd, ~] = redoxCalibrate(r, p);
                oxds(k) = oxd(1);
            end
            diffs = diff(double(oxds));
            tc.verifyTrue(all(diffs >= 0), ...
                'OxD should increase monotonically with ratio.');
        end

        function testOxdWithIFactorOneMatchesSimpleFormula(tc)
            % When iFactor=1, OxD = (R-rMin)/(rMax-rMin)  (linear normalisation)
            p = tc.defaultParams(); p.iFactor = 1;
            rMin = p.rMin; rMax = p.rMax;
            Rval = single(0.55);
            expected = (Rval - rMin) / (rMax - rMin);
            expected = max(0.001, min(0.999, double(expected)));
            ratioIn = tc.uniformRatio(4, 4, 2, Rval);
            [oxd, ~] = redoxCalibrate(ratioIn, p);
            tc.verifyEqual(double(oxd(1,1,1,1,1)), expected, 'AbsTol', 1e-4, ...
                'iFactor=1: OxD should equal linear normalisation.');
        end

    end

    % ------------------------------------------------------------------
    % Eh: Nernst equation analytical checks
    % ------------------------------------------------------------------
    methods (Test)

        function testEhAtMidpointEqualsE0prime(tc)
            % OxD = 0.5 → ln(1/1) = 0 → Eh = E0'
            % With pH=pH0, E0' = midpoint = -280 mV
            % Find R that gives OxD=0.5: (R-rMin)/[delta*(rMax-R)+(R-rMin)] = 0.5
            % → R-rMin = delta*(rMax-R)/2 + (R-rMin)/2 ... solve numerically
            p = tc.defaultParams();
            % For iFactor=5.5: R_half = (rMin + delta*rMax)/(1 + delta)
            delta = p.iFactor;
            Rhalf = (p.rMin + delta * p.rMax) / (1 + delta);

            ratioIn = tc.uniformRatio(4, 4, 2, single(Rhalf));
            [oxd, eh] = redoxCalibrate(ratioIn, p);

            tc.verifyEqual(double(oxd(1,1,1,1,1)), 0.5, 'AbsTol', 0.01, ...
                'OxD at Rhalf should be 0.5.');
            tc.verifyEqual(double(eh(1,1,1,1,1)), double(p.midpoint), 'AbsTol', 0.5, ...
                'Eh at OxD=0.5 should equal midpoint E°''.');
        end

        function testEhMoreOxidisedIsLessNegative(tc)
            % Higher OxD → more oxidised → Eh is less negative (less reducing)
            p = tc.defaultParams();
            rLow  = tc.uniformRatio(2, 2, 2, single(0.3));   % more reduced
            rHigh = tc.uniformRatio(2, 2, 2, single(0.7));   % more oxidised
            [~, ehLow]  = redoxCalibrate(rLow,  p);
            [~, ehHigh] = redoxCalibrate(rHigh, p);
            tc.verifyGreaterThan(double(ehHigh(1)), double(ehLow(1)), ...
                'More oxidised (higher R) should give less negative Eh.');
        end

        function testPHCorrectionShiftsMidpoint(tc)
            % pH = pH0 + 1 → E0' shifts by -(2.303*R*T/F) per pH unit ≈ 59 mV at 293K
            p0 = tc.defaultParams();  p0.pH = 7.0;  p0.pH0 = 7.0;
            p1 = tc.defaultParams();  p1.pH = 8.0;  p1.pH0 = 7.0;

            % Use Rhalf for p0 (OxD≈0.5) so that Eh ≈ E0'
            delta = p0.iFactor;
            Rhalf = single((p0.rMin + delta * p0.rMax) / (1 + delta));

            [~, eh0] = redoxCalibrate(tc.uniformRatio(2,2,2, Rhalf), p0);
            [~, eh1] = redoxCalibrate(tc.uniformRatio(2,2,2, Rhalf), p1);

            % Expected shift: −2.303 * 8.314 * 293 / 96485 * 1000 ≈ −58.2 mV per pH unit
            expectedShift = -(2.303 * 8.314 * 293 / 96485) * 1000;
            actualShift   = double(eh1(1)) - double(eh0(1));
            tc.verifyEqual(actualShift, expectedShift, 'AbsTol', 1.0, ...
                'pH shift should match Nernst one-electron pH correction.');
        end

        function testNernstSymmetry(tc)
            % Eh(OxD) + Eh(1-OxD) should equal 2*E0' (symmetric around E0')
            p = tc.defaultParams();
            p.iFactor = 1;   % linear normalisation simplifies R↔OxD mapping
            rMin = p.rMin; rMax = p.rMax;
            % R1 gives OxD=0.3, R2 gives OxD=0.7
            oxd1 = 0.3; oxd2 = 0.7;
            R1 = single(rMin + oxd1*(rMax-rMin));
            R2 = single(rMin + oxd2*(rMax-rMin));
            [~, eh1] = redoxCalibrate(tc.uniformRatio(2,2,2,R1), p);
            [~, eh2] = redoxCalibrate(tc.uniformRatio(2,2,2,R2), p);
            tc.verifyEqual(double(eh1(1)+eh2(1)), 2*double(p.midpoint), 'AbsTol', 0.5, ...
                'Nernst: Eh(OxD) + Eh(1-OxD) should equal 2*E0''.');
        end

    end

    % ------------------------------------------------------------------
    % No Inf/NaN in typical operating range
    % ------------------------------------------------------------------
    methods (Test)

        function testNoInfOrNaNInNormalRange(tc)
            p = tc.defaultParams();
            ratioIn = tc.gradientRatio(4, 20, 3, p.rMin + 0.01, p.rMax - 0.01);
            [oxd, eh] = redoxCalibrate(ratioIn, p);
            tc.verifyFalse(any(isinf(oxd(:))) || any(isnan(oxd(:))), ...
                'OxD should have no Inf or NaN in normal operating range.');
            tc.verifyFalse(any(isinf(eh(:))) || any(isnan(eh(:))), ...
                'Eh should have no Inf or NaN in normal operating range.');
        end

    end

end
