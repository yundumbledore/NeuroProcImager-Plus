# NeuroProcImager-Plus
#### Yun Zhao, Levin Kuhlmann (Monash University, Australia), Email: yun.zhao@monash.edu, levin.kuhlmann@monash.edu

**NeuroProcImager-Plus** is associated with the manuscript *Cortical local dynamics, connectivity and stability correlates of global consciousness states*.

**NeuroProcImager-Plus**, an extension of **NeuroProcImager** ([Github link](https://github.com/yundumbledore/NeuroProcImager/tree/main), [Neuroimage Paper link](https://www.sciencedirect.com/science/article/pii/S1053811922007078)), explores neurophysiological mechanisms of brain functions regarding effective connectivity between brain regions, regional neurophysiological variables, and dynamic cortical stability analysis.


## Methods

### Whole-cortex model

## Demonstration
Here we provide a demonstration to enable reviewers and readers to get in touch with NeuroProcImager-Plus.

**Showcase: Dynamic cortical stability analysis under xenon-induced asymptotic loss of consciousness**

This showcase calculates and shows the time-evolving stability of the cerebral cortex during an xenon-induced anaesthesia experiment where a subject started with wakeful state and slowly lost responsiveness and finally recovered from unconsciousness.

To run this Case, one needs to download "regional_variable_estimates_S11.mat" from [Google Drive](https://drive.google.com/drive/folders/1i8ZqNcqIbl0AMgG1JY3nuSUMqaBCREqD?usp=sharing) and put it in the /Data folder. The time-evolving neurophysiological variables in each region were estimated using [NeuroProcImager](https://github.com/yundumbledore/NeuroProcImager/tree/main) by fitting a neural mass model to each MEG source time series. This demonstration further estimates the inter-regional connectivity parameters, calculates the Jacobi matrices of the whole-cortex model system, and shows the time-evolving eigenvalue spectrum and the number of positive eigenvalues.

This demonstration will run overnight because the recording is half an hour long and the sampling rate is 150Hz. Run time highly depends on the hardware of the CPU workstation.

## Adaptation to your data
Users can prepare EEG or MEG data of the whole brain or part of the brain area, and use our [NeuroProcImager](https://github.com/yundumbledore/NeuroProcImager/tree/main) to estimate the regional neurophysiological variables of each brain region (that is, the region corresponding to each MEG or EEG time series). Assuming there are N regions (i.e., N time series) of length t, and each region has n neurophysiological variables, the user needs to save the estimated values ​​in the format of n x t x N and name it "xi_hat_list".
