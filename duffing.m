function xOut = duffing(tIn,xIn)

    global epsilon
    xOut = zeros(2,1);
    xOut(1) = xIn(2);
    xOut(2) = xIn(1) - xIn(1)^3 + epsilon*sin(tIn);
    
end