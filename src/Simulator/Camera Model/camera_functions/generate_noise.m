function noise = generate_noise(waveform, snr)
    %the lower signal to noise ratio, the more noisy function
    waveLength = length(waveform);
    signalPower = (sum(abs(waveform).^2))/waveLength;
    noisePower =  signalPower/(10^(snr/10));
    %we want random amount of additon/subtraction: range (-1,1)
    randomArr = (2 .* rand(1,waveLength)) - 1;
    noise = sqrt(noisePower) * randomArr;
end
