function [freq_error, correlation_peak] = blind_doppler_discriminator(block1, block2, fs)
    % Estimates frequency offset between two blocks using spectral cross-correlation.
    % This is a frequency-domain equivalent of a time-domain correlation.
    % The peak of the cross-correlation in the frequency domain corresponds
    % to the frequency offset.

    N = length(block1);
    
    % Compute FFTs
    FFT1 = fft(block1);
    FFT2 = fft(block2);
    
    % Spectral cross-correlation
    % Multiply the spectrum of one by the conjugate of the other
    X_spec = FFT1 .* conj(FFT2);
    
    % The IFFT of this product is the circular correlation
    x_corr = ifft(X_spec);
    
    % Find the peak of the correlation. The location of the peak gives the
    % time-domain shift, which relates to a phase ramp in the frequency domain.
    % A more direct method is to correlate the magnitude spectra.
    % Let's use the method proposed in many blind algorithms: find the frequency
    % shift that maximizes the correlation.
    
    % Simpler approach: find peak of frequency-domain correlation
    [~, idx] = max(abs(x_corr));
    
    % The peak location 'idx' corresponds to a time shift. A time shift of
    % 'tau' seconds corresponds to a linear phase ramp of e^(-j*2*pi*f*tau)
    % in the frequency domain. We are looking for the frequency offset.
    % A simpler interpretation is that the phase of the correlation at lag 0
    % is related to the frequency offset. Let's use an even more robust method
    % based on the FFT shift property. A frequency shift 'df' in time domain
    % becomes a circular shift in the frequency domain.
    
    % Simplified implementation for demonstration:
    % We will correlate the signals and find the phase of the peak.
    % The rate of change of phase gives frequency.
    % However, the paper suggests a "frequency-based Doppler discriminator".
    % Let's implement that directly. We find the shift required to align FFT1 and FFT2.
    
    % Re-implementing based on "spectral cross-correlation" for frequency offset:
    % Correlate the magnitude of the FFTs.
    corr_fft_mag = xcorr(abs(FFT1), abs(FFT2));
    [~, max_idx] = max(corr_fft_mag);
    
    % The center of the correlation result corresponds to zero shift
    center_idx = length(corr_fft_mag)/2 + 1;
    
    % The shift in FFT bins is the difference from the center
    bin_shift = max_idx - center_idx;
    
    % Convert bin shift to frequency error
    freq_resolution = fs / N;
    freq_error = bin_shift * freq_resolution;
    
    correlation_peak = max(abs(x_corr));
end
