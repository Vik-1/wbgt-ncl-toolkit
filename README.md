# WBGT NCL Toolkit                                                                           This NCL-based toolkit calculates **Wet-Bulb Globe Temperature (WBGT)** using:

- **Stull (2011)** approximation for wet-bulb temperature
- **Liljegren et al. (2008)** radiative-convective model for globe temperature
- **ISO 7243-compliant** WBGT formula for outdoor exposure

### 📦 Requirements

- NCL (version 6.6 or later recommended)
- NetCDF input files with:
  - Air temperature (°C)
  - Relative humidity (%)
  - Solar radiation (W/m²)
  - Wind speed (m/s)

### 🚀 Usage

```bash
ncl wbgt_improved.ncl

