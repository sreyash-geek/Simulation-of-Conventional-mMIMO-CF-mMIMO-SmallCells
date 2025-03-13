**Simulation of Cellular mMIMO, Cell-free mMIMO & Small cells for SINR, SNR & Energy/Spectral Efficiency Analysis.** 

This project contains three MATLAB scripts designed to simulate, analyze, and visualize advanced wireless network scenarios. 
The codes demonstrate various aspects of wireless communication performance, from spectral/energy efficiency analysis, channel modeling and SINR computations to innovative 3D visualizations and animations. 

**Overview:**

This script simulates and compares 3 different network architectures:
1) Cellular Massive MIMO
2) Cellular Small Cells
3) Cell-Free Networks
   

**Key Features:**

_**1) Accurate SINR Computation:**_
Uses MMSE combining to reliably calculate SINRs for the 3 network architectures.

_**2) Spectral Efficiency Analysis:**_
Computes the spectral efficiency (bit/s/Hz) and provides a direct comparison across setups.

_**3) Channel Modeling:**_
Implements free-space propagation with realistic path loss and phase shift calculations.

_**4) Dynamic Visualization:**_
Includes both static and animated 3D plots that showcase how normalized SNR and phase differences evolve over space and time.

_**5) Energy Efficiency Metric:**_ Implements a power consumption model for each architecture and computes energy efficiency (bit/s/Hz per Watt).

_**6) Comparative Visualization:**_ Presents a bar chart to easily compare the energy efficiency across the different network architectures.


**Plots:**

1) **_CDF of SINR_** : Displays cumulative distribution function (CDF) plots of SINR to compare network performance across architectures.
2) **_CCDF Plot_** : Complementary CDF graphs give good insights on reliability to a threshold.
3) **_2D Spatial Distribution Plots_**: Visualizes the positions of antennas and UEs to relate network geometry with performance.
4) **_Phase difference surface plot_** : For a selected antenna (first), to visualize phase variations.
5) **_Animated 3D plot_** : Shows how the normalized SNR evolves as the target moves along a predefined trajectory.
6) **_Bar plot for Energy Efficiency_** : Compares the EE metric across the Cell-Free and Cellular massive MIMO architectures.
7) **_CDF Plot for Spectral Efficiency_** : Compares the Spectral Efficiency distribution across the Cell-Free and Cellular massive MIMO setups and all realizations.
