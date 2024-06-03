Linear Predictive Coding (LPC) creates a model that can be used to predict the signal from it's linear combination of it's previous LPC coefficients.
The code is implemented in accordance to the LPC processor, which consists of 5 major blocks which are listed below: 

1. Pre-emphasis: It spectrally flattens the signal, to make it less susceptible to finite precision effects later in the signal processing.
2. Frame Blocking: The output of the Pre-emphasis block is blocked into frames of N samples each with adjacent frames separated by M samples. 
3. Windowing: Each frame is windowed to minimize the signal discontinuties at the beginning and end of the frame, typical window used is Hamming window.
4. Autocorrelation Analysis: Each frame of windowed signal is autocorrelated.
5. LPC Analysis: Each frame of p+1 autocorrelateds (p is order of LPC Analysis) is convereted into the LPC parameter set which may consist of LPC coefficients/reflection coefficients (PARCOR)/log area ratio coefficients/Cepstral coefficients.

NOTE: If in any case M>N or no overlapping between adjacent frames, some of the speech signal might be lost & LPC spectral estimates of adjacent frames will contain a noisy component whose magnitude is directly proportional to the 'M', this situation is intolerable in any LPC Analysis for speech recognition.

Values considered here: M = 100, N = 653 and p = 10.
