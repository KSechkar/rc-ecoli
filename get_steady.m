function last_x=get_steady(sim,Delta,Max_iter)
    % iterate integrations until steady state reached
    last_x_old=zeros(1,(size(sim.init_conditions,1)+size(sim.het.init_conditions,1))); % start with growth rate assumed to be zero
    for i=1:Max_iter
        % integrate
        sim = sim.simulate_model;
           
        last_x=cat(2,sim.x(end,:),sim.het.x(end,:));
        
        %[last_x_old; last_x]
        if(norm(last_x-last_x_old)<Delta) % if SS reached, quit
            disp(['SS reached in ',num2str(i),' iterations'])
            break
        else % if not, continue integrating
            if(i==Max_iter)
                disp('Warning! SS not reached yet')
            end
            last_x_old=last_x;
            % mRNAs
            sim.init_conditions('m_a')=last_x(1);
            sim.init_conditions('m_r')=last_x(2);
            % proteins
            sim.init_conditions('p_a')=last_x(3);
            sim.init_conditions('R')=last_x(4);
            % tRNAs
            sim.init_conditions('tc')=last_x(5);
            sim.init_conditions('tu')=last_x(6);
            % inactivated ribosomes
            sim.init_conditions('Bcm')=last_x(7);
            % heterologous
            x_het=last_x(8:end);
            for j=1:sim.het.num_genes
                % mRNA
                sim.het.init_conditions(['m_',sim.het.names{j}])=x_het(j);
                % protein
                sim.het.init_conditions(['p_',sim.het.names{j}])=x_het(sim.het.num_genes+j);
            end
        end
    end
end