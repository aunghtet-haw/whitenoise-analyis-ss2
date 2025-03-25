# whitenoise-analyis-ss2
# **Band-Limited White Noise Analysis in MATLAB**

This MATLAB project generates, analyzes, and visualizes a **band-limited Gaussian white noise** signal. It demonstrates fundamental signal processing techniques, including power spectral estimation, autocorrelation, and power density spectrum analysis.

## **Features**
- **White Noise Generation**: Creates a Gaussian white noise signal with configurable mean and standard deviation.
- **Time-Domain Analysis**: Plots the noise signal and its statistical distribution.
- **Autocorrelation Estimation**: Computes and visualizes the autocorrelation function.
- **Power Spectrum Analysis**: Implements three methods for spectral estimation:
  1. **DFT of the autocorrelation function**
  2. **DFT of the signal**
  3. **MATLAB’s `pwelch()` function**
- **Power Calculation**: Compares different approaches to estimating the mean power of the signal.
- **Power Spectral Density (PSD)**: Uses `pwelch()` with different windowing techniques for spectral density estimation.

## **Installation & Usage**
### **Requirements**
- MATLAB (Tested on R2021a and later)

### **Run the Script**
1. Clone the repository:
   ```bash
   git clone https://github.com/aunghtet-haw/whitenoise-analyis-ss2.git
   cd whitenoise-analyis-ss2
   ```
2. Open MATLAB and run the script:
   ```matlab
    	whitenoiseanalysis.m
   ```

## **Generated Plots**
The script generates multiple plots, including:
- **White noise time-domain representation**
- **Histogram of sample values**
- **Autocorrelation function**
- **Power spectrum estimates using different methods**
- **Power spectral density with different windowing techniques**

## **License**
This project is open-source under the MIT License. Feel free to use and modify it for research or educational purposes.

