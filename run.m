clear all
close all
clc
%%
% Author: Yun Zhao (Monash University, Australia)
% Email: yun.zhao@monash.edu
% -------------------------------------------------------------------------
% To run this demo, one needs to install "signal processing toolbox" 
% in MATLAB.
%
% -------------------------------------------------------------------------
% Here we provide two showcases to enable reviewers and readers to get
% in touch with NeuroProcImager-Plus.
%
% Case 1: Whole-cortex model parameters estimation
%
% This showcase runs a specialized inference method to estimate parameters
% in the whole-cortex model from an example data file. The example data is
% a anaesthesia MEG recording of 78 source time series from one subject.
%
% To run Case 1, one needs to put "parameters estimation" in the pipeline 
% below.
%
% Case 2: Dynamic cortical stability
%
% This showcase performs stability analysis of the whole-cortex model in
% each time window and shows the result in the way we present it in the
% manuscript.
%
% To run Case 2, one needs to put "cortical stability" in the pipeline
% below.
%
% Note: Case 1 should be run prior to Case 2.
% -------------------------------------------------------------------------
% Here we show time statistics of running NeuroProcImager-Plus on a CPU
% workstation with 32 cores and 155 GB memory.
%
% Case 1 finished in around 5 minutes.
% Case 2 finished in around 12.5 hours.
%
% Run time highly depends on the hardware of the machine.
% -------------------------------------------------------------------------
%
%% Add source to MATLAB path
addpath(genpath('./Sources'))

%% Choose showcase(s) to do
tasks = ["parameters estimation"]; % which showcase to run. two showcases: "parameters estimation", "cortical stability".

%% Run
try
    main(tasks)
catch ME
    disp(ME)
end