function sos=advanced_sos(theta,data,sim,Delta,Max_iter)  % sum of squares function for DRAM fitting   
    ymodel=advanced_modelfun(theta,data.xdata,sim,Delta,Max_iter);
    % get sum of squares
    sos = sum((ymodel - data.ydata).^2);
    disp(['a_act=',num2str(log(theta))])
    disp(['SOS=',num2str(sos)])
end