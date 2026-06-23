# Component-wise-ECG-Modelling-and-Parameter-Extraction

A MATLAB framework for **model‑based ECG morphology analysis and synthetic ECG generation** using **phase‑domain representation and Gaussian mixture modeling**.

Unlike many ECG modeling approaches that approximate the entire waveform with a single model, this project **models each physiological ECG component separately**:

• **P wave**  
• **QRS complex**  
• **T wave**

Each component is independently modeled using a **Gaussian mixture representation optimized via Particle Swarm Optimization (PSO)**.  
The framework then reconstructs **synthetic versions of each wave individually** as well as the **complete synthetic ECG waveform**.

This provides an interpretable parametric representation of the cardiac cycle.

---

# Overview

Electrocardiogram (ECG) signals consist of repeating cardiac cycles whose morphology varies due to heart rate variability, noise, and physiological changes. Direct modeling in the time domain is therefore difficult because each heartbeat has a different duration.

To address this, the ECG signal is transformed into the **cardiac phase domain**, where each beat is mapped into a normalized phase interval:

- start of beat → \( -\pi \)  
- **R‑peak → 0**  
- end of beat → \( \pi \)

This transformation aligns cardiac cycles regardless of their duration, allowing accurate estimation of the **mean ECG morphology**.

The mean ECG waveform is then decomposed into its physiological components (**P, QRS, and T waves**), and each component is modeled using **Gaussian kernels**.

---

# Key Features

• ECG modeling in the **cardiac phase domain**  
• **Automatic R‑peak detection** using the Pan–Tompkins algorithm  
• Optional **nonlinear phase alignment using Dynamic Time Warping (DTW)**  
• **Mean ECG morphology extraction** across heartbeats  
• **Separate modeling of P wave, QRS complex, and T wave**  
• **Gaussian mixture decomposition** of each ECG component  
• Parameter estimation via **Particle Swarm Optimization (PSO)**  
• Generation of **synthetic ECG components and full ECG waveform**

---

# Algorithm Pipeline

The processing pipeline of the framework is illustrated below.

1. **ECG Signal Loading**

The ECG signal is loaded from a MATLAB `.mat` dataset:

```
x  → ECG signal matrix  
fs → sampling frequency
```

Example:

```matlab
ecg = x(1,:);
```

---

2. **R‑Peak Detection**

R‑peaks are detected using the **Pan–Tompkins QRS detection algorithm**:

```matlab
[qrs_positions] = pantompkins_qrs(abs(ecg),fs);
```

These R‑peaks define the boundaries of individual cardiac cycles.

---

3. **Cardiac Phase Calculation**

Each ECG sample is mapped to a **cardiac phase** between consecutive R‑peaks.

```matlab
[Linearphase,~] = calculate_linear_phase_ver2(qrs_positions,length_sig,fs);
```

The phase is wrapped to the interval:

\[
[-\pi , \pi]
\]

where

\[
\theta = 0
\]

corresponds to the **R‑peak**.

---

4. **Nonlinear Phase Alignment (Optional)**

To improve alignment between heartbeats, **Dynamic Time Warping (DTW)** can be applied to obtain a nonlinear phase trajectory.

This step compensates for morphological variations between beats and improves mean ECG estimation.

Output:

```
NonlinearPhase
```

---

5. **Mean ECG Morphology Extraction**

The ECG waveform is averaged in the phase domain using **phase binning**.

Example configuration:

```
Number of phase bins = 200
```

For each phase bin:

• mean ECG amplitude is computed  
• standard deviation is estimated

Function used:

```
ecgsd_extractor_ver1
```

Outputs include:

• `ECGmean` – mean ECG using linear phase  
• `ECGmean_nonlinear_phase` – mean ECG using nonlinear phase  
• `ECGsd` – phase‑dependent variability

---

# Component‑Wise ECG Wave Modeling

The mean ECG waveform is divided into physiologically meaningful regions in the phase domain.

## P Wave Region

Phase interval:

\[
-\frac{\pi}{2} < \theta < -\frac{\pi}{6}
\]

The P wave is modeled using a **Gaussian mixture**:

```
Number of Gaussian kernels = 3
```

Synthetic output:

```
Synthetic_P
```

---

## QRS Complex Region

Phase interval:

\[
-\frac{\pi}{6} < \theta < \frac{\pi}{6}
\]

The QRS complex is modeled using:

```
Number of Gaussian kernels = 5
```

Synthetic output:

```
Synthetic_QRS
```

---

## T Wave Region

Phase interval:

\[
\frac{\pi}{6} < \theta < \frac{2\pi}{3}
\]

The T wave is modeled using:

```
Number of Gaussian kernels = 3
```

Synthetic output:

```
Synthetic_T
```

---

# Gaussian Kernel Model

Each ECG component is approximated using a set of Gaussian kernels:

\[
G(\theta) = a_i \, \exp \left(-\frac{(\theta-\theta_i)^2}{2b_i^2}\right)
\]

where:

• \(a_i\) → amplitude  
• \(b_i\) → width  
• \( \theta_i \) → phase location of the kernel

The sum of these kernels reconstructs the morphology of each ECG component.

---

# Parameter Estimation

Gaussian parameters are estimated by minimizing the reconstruction error between the **mean ECG component** and its **Gaussian approximation**.

Optimization method:

```
Particle Swarm Optimization (PSO)
```

Typical configuration:

```
SwarmSize = 200 – 2000
MaxIterations = 200
```

The optimization estimates the parameters:

```
ai_Pwave, bi_Pwave, tetai_Pwave
ai_QRSwave, bi_QRSwave, tetai_QRSwave
ai_Twave, bi_Twave, tetai_Twave
```

---

# Synthetic ECG Generation

After estimating the Gaussian parameters, the framework generates **synthetic ECG components**:

```
Synthetic_P
Synthetic_QRS
Synthetic_T
```

These synthetic waves are then combined to produce the complete ECG model:

```
Synthetic_ECG = Synthetic_P + Synthetic_QRS + Synthetic_T
```

This provides a **fully parametric representation of the ECG waveform**.

---

# Visualization Outputs

The script generates several diagnostic figures:

• ECG signal with detected **R‑peaks**  
• Mean ECG morphology (linear vs nonlinear phase)  
• Gaussian reconstruction of **P wave**  
• Gaussian reconstruction of **QRS complex**  
• Gaussian reconstruction of **T wave**  
• Synthetic ECG mean vs original ECG mean  
• Synthetic ECG vs original ECG signal

---

# Requirements

MATLAB toolboxes:

• Signal Processing Toolbox  
• Global Optimization Toolbox  
• Statistics Toolbox

Required helper functions:

```
pantompkins_qrs.m
calculate_linear_phase_ver2.m
ecgsd_extractor_ver1.m
```

Optional:

```
dtw()
```

for nonlinear phase alignment.

---

# Applications

This framework can be used for:

• ECG morphology analysis  
• synthetic ECG signal generation  
• ECG compression and parametric modeling  
• cardiac waveform simulation  
• biomedical signal processing research  
• model‑based ECG denoising

---


The result is a **compact and interpretable parametric model of the ECG waveform** suitable for analysis, simulation, and research applications.
:::
