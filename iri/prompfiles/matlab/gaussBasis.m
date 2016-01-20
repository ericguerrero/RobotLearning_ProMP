function val = gaussian1D(x, mu, sigma)
%  x     - Input vector
%  mu    - Mean
%  sigma - Standard deviation

    % Evaluate a 1D gaussian.
    val = exp(-(x - mu).^2 ./ (2 * sigma^2));
end