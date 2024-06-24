function result = diffusion_step(dt, data, diffusion_rate, Flk)
    data_fft = fft(data);
    data_fft = data_fft .* exp(-dt * Flk' * diffusion_rate);
    result = real(ifft(data_fft));
end