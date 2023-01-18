function G = antenna_gain(theta_Rx, theta_Tx, Go)
    % This function implements the antenna gain model from ITU-R (F.1336)
    % theta_Rx: elevation angle (degrees) of Rx relative to Tx w.r.t. ground plane
    % theta_Tx: antenna tilt angle (degrees) of Tx w.r.t. ground plane
    % Go      : the maximum gain in or near the horizontal plane (azimuth plane)
    %           (reads G naught)
    % return G: antenna gain in power (function should convert from dBi)

    if nargin < 3
        Go = 1;
    end
    theta = abs(theta_Rx - theta_Tx);
    assert((theta>=0) && (theta<=90),"Relative elevation angle must be within 0 and 90");
    theta_3 = 107.6*10^(-Go/10);
    
    if theta < theta_3
        G_dBi = Go - 12*(theta/theta_3)^2;
    else
        G_dBi = Go - 12 + 10*log10((theta/theta_3)^(-1.5));
    end
    
    G = 10^(G_dBi/10);
end