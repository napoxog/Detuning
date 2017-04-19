# Detuning
Python script for De-tuning workflow simplification in Paradigm's QSI Software suite

The workflow basics can be found [here](http://www.epmag.com/seismic-amplitudes-benefit-seismic-trace-detuning-1017366#p=full)

Current implementation provides:
* Import of column-based file with extracted Seismic Wavelet
* Spline interpolation of the Wavelet
* Estimation of "tuned" composite amplitude using Wavelet 
* Estimation of "no-tune" amplitude using position of maximum of composite amplitude
* Estimation of transfer function to be used to create detune operator map
* Generation of the ASCII file with calculations results
* Vizualization of the estimations using [gnuplot](https://sourceforge.net/projects/gnuplot/files/gnuplot/) tool 
