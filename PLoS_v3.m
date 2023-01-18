function p = PLoS_v3(r,h_tx,h_rx,a1,a2,a3)
    % This function implements the LoS probability model from the paper
    % "Blockage Modeling for Inter-layer UAVs Communications in Urban
    % Environments" 
    % param r    : horizontal distance between Tx and Rx (m)
    % param h_tx : height of Tx
    % param h_rx : height of Rx
    % param a1   : ratio of built up land area to total land area
    % param a2   : density of buildings (in /m^2) (note that the paper gives this value in /km^2)
    % param a3   : scale parameter that is used to characterize the building height
    %              distributions using a Rayleigh distribution
    %
    % returns p  : LoS probability between Tx and Rx 

    % The values of a1, a2 and a3 can be found from: Elevation Dependent Shadowing 
    % Model for Mobile Communications via High Altitude Platforms in Built-Up Areas
    % DOI:
    % For suburban       : a1 = 0.1, a2 = 7.5e-4, a3 = 8
    % For urban          : a1 = 0.3, a2 = 5e-4, a3 = 15
    % For dense urban    : a1 = 0.5, a2 = 3e-4, a3 = 20
    % For urban high-rise: a1 = 0.5, a2 = 3e-4, a3 = 50
    
    delta_h = h_tx - h_rx;
    %pow_factor = r*sqrt(a1*a2); % Note that this power factor has been
    %changed to follow the paper mentioned here, Blockage Modelling...
    pow_factor = 2*r*sqrt(a1*a2/pi) + a1;
    if (delta_h == 0)
        p = (1 - exp((-(h_tx)^2)/(2*a3^2)))^pow_factor;
    else
        if (delta_h<0)
            h1=h_rx;
            h2=h_tx;
        else 
            h1=h_tx;
            h2=h_rx;
        end
        delta_h = abs(delta_h);
        p = (1 - (sqrt(2*pi)*a3/delta_h)*abs((qfunc(h1/a3) - qfunc(h2/a3))))^pow_factor;
    end
end