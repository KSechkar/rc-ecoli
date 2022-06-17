%% rc_ecoli_modelfun.m
% retrieves sum of squared errors between experimentally determined
% steady-state growth rate andribosomal mass fraction predicted by the model
% for parameter fitting using DRAM

%%

function sos=rc_ecoli_sos(theta,data,sim,Delta,Max_iter)
    % get the values
    ymodel=rc_ecoli_modelfun(theta,data.xdata,sim,Delta,Max_iter);

    % get sum of squares
    sos = sum((ymodel - data.ydata).^2);

    % print the results
    disp(['SOS=',num2str(sos)])
end