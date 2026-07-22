classdef TestDiagnosticGroupCompletion < matlab.unittest.TestCase
    properties
        tempFolder
        simulationPath
        diagnosticsPath
        wvd
    end

    methods (TestMethodSetup)
        function createStaticDiagnosticFixture(testCase)
            testCase.tempFolder = string(tempname);
            mkdir(testCase.tempFolder);
            testCase.simulationPath = fullfile(testCase.tempFolder,"simulation.nc");
            testCase.diagnosticsPath = fullfile(testCase.tempFolder,"simulation-diagnostics.nc");

            N0 = 5.2e-3;
            wvt = WVTransformHydrostatic([1e4,1e4,1000],[8,8,5],N2=@(z) N0*N0*ones(size(z)),latitude=30);
            wvt.addForcing(WVAdaptiveDamping(wvt));
            model = WVModel(wvt);
            outputFile = model.createNetCDFFileForModelOutput(testCase.simulationPath,outputInterval=1);
            outputFile.outputTimesForIntegrationPeriod(0,2);
            for t = 0:2
                wvt.t = t;
                outputFile.writeTimeStepToOutputFile(t);
            end
            outputFile.closeNetCDFFile();

            testCase.wvd = WVDiagnostics(testCase.simulationPath);
            testCase.wvd.createDiagnosticsFile();
        end
    end

    methods (TestMethodTeardown)
        function removeStaticDiagnosticFixture(testCase)
            if ~isempty(testCase.wvd)
                if ~isempty(testCase.wvd.diagfile)
                    testCase.wvd.diagfile.close();
                end
                if ~isempty(testCase.wvd.wvfile)
                    testCase.wvd.wvfile.close();
                end
            end
            if isfolder(testCase.tempFolder)
                rmdir(testCase.tempFolder,"s");
            end
        end
    end

    methods (Test)
        function testRepeatedGeostrophicCallPreservesCompletedData(testCase)
            testCase.wvd.createGeostrophicFluxGroup();
            group = testCase.wvd.diagfile.groupWithName("geostrophic-flux");
            ggg = group.variableWithName("ggg");
            sentinel = ggg.valueAlongDimensionAtIndex("t",1) + 12345;
            ggg.setValueAlongDimensionAtIndex(sentinel,"t",1);

            testCase.wvd.createGeostrophicFluxGroup();

            testCase.verifyEqual(ggg.valueAlongDimensionAtIndex("t",1),sentinel);
            testCase.verifyWarningFree(@() testCase.wvd.createGeostrophicFluxGroup(timeIndices=[]));
        end

        function testGeostrophicMissingAndOverwritePolicies(testCase)
            testCase.wvd.createGeostrophicFluxGroup();
            group = testCase.wvd.diagfile.groupWithName("geostrophic-flux");
            ggg = group.variableWithName("ggg");
            expectedFirst = ggg.valueAlongDimensionAtIndex("t",1);
            expectedSecond = ggg.valueAlongDimensionAtIndex("t",2);
            sentinelFirst = expectedFirst + 12345;
            sentinelSecond = expectedSecond + 12345;
            ggg.setValueAlongDimensionAtIndex(sentinelFirst,"t",1);
            ggg.setValueAlongDimensionAtIndex(sentinelSecond,"t",2);
            TestDiagnosticGroupCompletion.markTimeIncomplete(group,2);

            testCase.wvd.createGeostrophicFluxGroup(timeIndices=[1 2]);

            testCase.verifyEqual(ggg.valueAlongDimensionAtIndex("t",1),sentinelFirst);
            testCase.verifyEqual(ggg.valueAlongDimensionAtIndex("t",2),expectedSecond,"AbsTol",1e-12);

            testCase.wvd.createGeostrophicFluxGroup(timeIndices=1,shouldOverwriteExisting=true);
            testCase.verifyEqual(ggg.valueAlongDimensionAtIndex("t",1),expectedFirst,"AbsTol",1e-12);
        end

        function testLegacyGeostrophicGroupMigration(testCase)
            diagfile = testCase.wvd.diagfile;
            group = diagfile.addGroup("geostrophic-flux");
            variableNames = ["ggg", "ggw", "ggw_tx", "wwg_tx"];
            dimensions = ["j", "kRadial", "t"];
            value = ones(length(testCase.wvd.j),length(testCase.wvd.kRadial));
            for iVariable = 1:length(variableNames)
                variable = group.addVariable(variableNames(iVariable),dimensions,type="double",isComplex=false);
                for outputIndex = 1:length(testCase.wvd.t_diag)
                    variable.setValueAlongDimensionAtIndex(value,"t",outputIndex);
                end
            end

            testCase.wvd.createGeostrophicFluxGroup();

            testCase.verifyEqual(group.readVariables("t"),testCase.wvd.t_diag);
            testCase.verifyEqual(group.variableWithName("ggg").valueAlongDimensionAtIndex("t",1),value);
        end

        function testRepeatedReservoirCallPreservesCompletedData(testCase)
            testCase.wvd.createReservoirGroup();
            group = testCase.wvd.diagfile.groupWithName("reservoir-damped-wave-geo");
            energy = group.variableWithName("E_1");
            sentinel = energy.valueAlongDimensionAtIndex("t",1) + 12345;
            energy.setValueAlongDimensionAtIndex(sentinel,"t",1);

            testCase.wvd.createReservoirGroup();

            testCase.verifyEqual(energy.valueAlongDimensionAtIndex("t",1),sentinel);
        end

        function testReservoirMissingAndOverwritePolicies(testCase)
            testCase.wvd.createReservoirGroup();
            group = testCase.wvd.diagfile.groupWithName("reservoir-damped-wave-geo");
            energy = group.variableWithName("E_1");
            expectedFirst = energy.valueAlongDimensionAtIndex("t",1);
            expectedSecond = energy.valueAlongDimensionAtIndex("t",2);
            sentinelFirst = expectedFirst + 12345;
            sentinelSecond = expectedSecond + 12345;
            energy.setValueAlongDimensionAtIndex(sentinelFirst,"t",1);
            energy.setValueAlongDimensionAtIndex(sentinelSecond,"t",2);
            TestDiagnosticGroupCompletion.markTimeIncomplete(group,2);

            testCase.wvd.createReservoirGroup(timeIndices=[1 2]);

            testCase.verifyEqual(energy.valueAlongDimensionAtIndex("t",1),sentinelFirst);
            testCase.verifyEqual(energy.valueAlongDimensionAtIndex("t",2),expectedSecond,"AbsTol",1e-12);

            testCase.wvd.createReservoirGroup(timeIndices=1,shouldOverwriteExisting=true);
            testCase.verifyEqual(energy.valueAlongDimensionAtIndex("t",1),expectedFirst,"AbsTol",1e-12);
        end

        function testLegacyReservoirGroupMigration(testCase)
            flowComponents = testCase.wvd.wvt.flowComponentWithName("wave");
            groupName = "legacy-reservoir";
            [triadVar,forcingVar,energyVar] = testCase.wvd.variablesForReservoirGroup(outputfile=testCase.wvd.diagfile,name=groupName,flowComponents=flowComponents);
            outputVariables = TestDiagnosticGroupCompletion.dictionaryVariables(triadVar,forcingVar,energyVar);
            for iVariable = 1:length(outputVariables)
                for outputIndex = 1:length(testCase.wvd.t_diag)
                    outputVariables{iVariable}.setValueAlongDimensionAtIndex(1,"t",outputIndex);
                end
            end

            testCase.wvd.createReservoirGroup(name=groupName,flowComponents=flowComponents);

            group = testCase.wvd.diagfile.groupWithName(groupName);
            testCase.verifyEqual(group.readVariables("t"),testCase.wvd.t_diag);
            testCase.verifyEqual(group.variableWithName("E_1").valueAlongDimensionAtIndex("t",1),1);
        end

        function testPartialSchemasAreRepaired(testCase)
            diagfile = testCase.wvd.diagfile;
            geostrophicGroup = diagfile.addGroup("geostrophic-flux");
            geostrophicGroup.addVariable("ggg",["j", "kRadial", "t"],type="double",isComplex=false);

            testCase.wvd.createGeostrophicFluxGroup(timeIndices=1);

            testCase.verifyTrue(all(geostrophicGroup.hasVariableWithName("ggw","ggw_tx","wwg_tx")));
            testCase.verifyEqual(geostrophicGroup.variableWithName("t").valueAlongDimensionAtIndex("t",1),testCase.wvd.t_diag(1));

            flowComponents = testCase.wvd.wvt.flowComponentWithName("wave");
            reservoirGroup = diagfile.addGroup("partial-reservoir");
            reservoirGroup.addAttribute("flow-components",reshape([flowComponents.name],[],1));
            reservoirGroup.addVariable("E_1","t",type="double",isComplex=false);

            testCase.wvd.createReservoirGroup(name="partial-reservoir",flowComponents=flowComponents,timeIndices=1);

            testCase.verifyTrue(reservoirGroup.hasVariableWithName("T_1_1_1"));
            testCase.verifyEqual(reservoirGroup.variableWithName("t").valueAlongDimensionAtIndex("t",1),testCase.wvd.t_diag(1));
        end

        function testReservoirDefinitionMismatchErrors(testCase)
            testCase.wvd.createReservoirGroup();
            differentFlowComponents = testCase.wvd.wvt.flowComponentWithName("wave");

            testCase.verifyError(@() testCase.wvd.createReservoirGroup(flowComponents=differentFlowComponents,timeIndices=[]),"WVDiagnostics:ReservoirGroupDefinitionMismatch");
        end
    end

    methods (Static, Access=private)
        function markTimeIncomplete(group,outputIndex)
            fillValue = netcdf.getConstant('NC_FILL_DOUBLE');
            group.variableWithName("t").setValueAlongDimensionAtIndex(fillValue,"t",outputIndex);
        end

        function outputVariables = dictionaryVariables(varargin)
            nVariables = 0;
            for iDictionary = 1:nargin
                nVariables = nVariables + length(varargin{iDictionary}.keys);
            end
            outputVariables = cell(nVariables,1);
            iOutputVariable = 0;
            for iDictionary = 1:nargin
                variableDictionary = varargin{iDictionary};
                variableNames = variableDictionary.keys;
                for iVariable = 1:length(variableNames)
                    iOutputVariable = iOutputVariable + 1;
                    outputVariables{iOutputVariable} = variableDictionary{variableNames(iVariable)};
                end
            end
        end
    end
end
