# Band-Limited White Noise Analysis and RLC Circuit Simulation

This MATLAB project as a part of Signal and System 2 course at HAW HH performs analysis on band-limited white noise, including signal generation, statistical analysis, power spectral density estimation, and system response through an RLC circuit filter. The project also compares theoretical and experimental results.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [License](#license)

## Overview
This project generates a Gaussian white noise signal and analyzes its properties in both time and frequency domains. It then applies an RLC filter to the signal and compares the theoretical response with the experimental response obtained via simulation.

## Features
- Generation of band-limited Gaussian white noise
- Time-domain visualization and histogram analysis
- Autocorrelation estimation
- Power spectral density (PSD) computation using different methods:
  - Direct Fourier Transform (DFT) of the autocorrelation function
  - FFT-based spectrum estimation
  - MATLAB's `pwelch()` function
- Analysis of the RLC circuit as a system:
  - Theoretical frequency response of the RLC circuit
  - Bode plot of the system
  - Simulation of discrete-time system response
- Comparison of theoretical and experimental system responses

## Installation
1. Ensure you have MATLAB installed on your system.
2. Clone this repository:
   ```sh
   git clone https://github.com/aunghtet-haw/whitenoise-analyis-ss2.git
   ```
3. Navigate to the project directory:
   ```sh
   cd whitenoise-analyis-ss2
   ```

## Usage
1. Open MATLAB and navigate to the project directory.
2. Run the main script:
   ```matlab
    whitenoise-analyis-ss2.m
   ```
3. The script will generate plots and numerical outputs related to noise analysis and system response.
4. Modify parameters such as noise variance, sampling rate, and RLC values to observe different effects.

## Results
This project generates various plots, including:
- Input signal time-domain waveform
- Histogram of input values
- Power spectral density using different methods
- Theoretical and experimental amplitude and phase responses of the RLC circuit
- Autocorrelation functions
- Comparison of system response with theoretical expectations

## License
This project is licensed under the MIT License. See the LICENSE file for details.

---
Feel free to contribute by submitting issues or pull requests!
