function offset = code_phase_tracker(input_block, ref_signal)
    % Finds the code phase offset between the input block and a reference
    % signal using circular cross-correlation.
    
    % Correlation can be efficiently computed using FFTs.
    % corr(a, b) = ifft(fft(a) .* conj(fft(b)))
    
    N = length(input_block);
    
    corr_result = ifft(fft(input_block) .* conj(fft(ref_signal)));
    
    % Find the peak of the correlation magnitude
    [~, idx] = max(abs(corr_result));
    
    % The index 'idx' corresponds to the circular shift.
    % MATLAB's index is 1-based. A peak at idx=1 means 0 shift.
    offset = idx - 1;
    
    % Handle wrapping around (for negative shifts)
    if offset > N/2
        offset = offset - N;
    end
end
