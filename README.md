## Testing The Lamp-Post and Wind Reverberation Models with XMM-Newton Observations of NGC 5506

These papes contain the data and codes associated with the article entitled: **[Testing The Lamp-Post and Wind Reverberation Models with XMM-Newton Observations of NGC 5506](https://arxiv.org/abs/1908.09862)**, published in the Astrophysical Journal.

### Abstract
> The lamp-post geometry is often used to model X-ray data of accreting black holes. Despite its simple assumptions, it has proven to be powerful in inferring fundamental black hole properties such as the spin. Early results of X-ray reverberations showed support for such a simple picture, though wind-reverberation models have also been shown to explain the observed delays. Here, we analyze new and old XMM-Newton observations of the variable Seyfert-1 galaxy NGC 5506 to test these models. The source shows an emission line feature around 6.7 keV that is delayed relative to harder and softer energy bands. The spectral feature can be modeled with either a weakly relativistic disk line or by scattering in distant material. By modeling both the spectral and timing signatures, we find that the reflection fraction needed to explain the lags is *larger* than observed in the time-averaged spectrum, ruling out both a static lamp-post and simple wind reverberation models.


### Description
The analysis is organized into several python notebooks, which sometimes call outside functions either from my [toolset package `aztools`](https://zoghbi-a.github.io/aztools/) or the helper scripts: 
- [`spec_helpers.py`](https://github.com/zoghbi-a/Testing-Reverberation-In-NGC5506/blob/master/spec_helpers.py): This contrain a collection of functions used in the spectral modeling. These are documented individually.
- [`timing_helpers.py`](https://github.com/zoghbi-a/Testing-Reverberation-In-NGC5506/blob/master/timing_helpers.pyy): Contains the functions used in the timing analysis, and these are mostly used in the [Timing](timing) (see bellow) notebook, and they are also documented individually
- [`fit.tcl`](https://github.com/zoghbi-a/Testing-Reverberation-In-NGC5506/blob/master/fit.tcl): This contrains a collection of tcl functions called from xspec to do the modeling. They are mostly called from the spectral modeling notebooks (see below).

A quick description of each notebooks is as follows:

- [Data](data.md): This is the data extraction tools, and contains the code for downloading, reducing the data, and extracting the spectra from XMM, Suzaku and NuSTAR used in the data.
- [Spec](spec.md): does the consistent spectral modeling. This is where the spectral results in the paper come from.
- [Timing](timing.md): Contains the code for extracting and modeling the power spectra, covariance and lags.

### Data Products
All the data products are available through the Open Science Framework at the [following link](https://osf.io/3unf2/files/). There are three files:
- [`xmm_spec.tgz`](hhttps://osf.io/uvk4t/): contains the spectral products and modeling from individual observations.
- [`xmm_timing.tgz`](https://osf.io/kusrt/): contains all the timing products.
