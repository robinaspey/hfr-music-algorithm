# hfr-music-algorithm
This code is a work in progress an analyses time series data from direction finding radars as produced by a well known
west-coast US company. A sample time series file is included for running the analysis.

The purpose of this is to validate reading time series data by running the data through a band-pass Butterworth filter
and then moving a buffer along the data series (something like 256 bytes) on which an FFT is performed. The analysis
is for loop 1 and loop 2 direction finding antenna and the transmit/receive whip antenna which has a 16.14MHz carrier
frequency.

As part of the validation process the cumulative FFT (time integrated) is built up for the entire 8192 byte snapshot.
Since the signal does not decay significantly over the duration of the snapshot (even with the noise removed) further
work is required.
