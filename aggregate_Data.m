

function out = aggregate_Data(in)

% This function groups the different interesting epoch combinations in a
% Matrix. Five different stimuli were presented in random order, leading 
% to a number of 25 possible combinations of stimulus pairs, 
% out of which 18 were interesting to us.
% Those are sort then into a matrix for later analysis.
% 
% This example is used by the script "Data_Analysis". This analysis 
% is a common analysis pattern which is used in the lab. I adjusted 
% it to aggregate this special kind of stimulus.   
%  
% Some function and variable names have been modified to the original
% script.
%

% all interesting epoch-combinations (18 out of 25 possible ones)
N_EPOCHS = 18;

%figure epochlength out
epochlengths = zeros(1,length(in));

for ii=1:length(in)
    s = in(ii).istim;
    
    % Find beginning indices
    ds = diff(s);
    b_inds = find(ds)+1;
    duration = zeros(1,length(b_inds));
    
    % Calculate each epochduration for this ROI
    for kk =1:length(b_inds)-1
        duration(kk) = b_inds(kk+1)-b_inds(kk);
    end
    
    % Don't take the mean, because sometimes the same epoch is displayed
    % twice, so the epochlength seems to be double. Mode takes 'most
    % common' value
    epochlengths(ii) = mode(duration);
end

epochlength = round(mean(epochlengths));

% For each ROI (length(in)) the average of
% all interesting combinations (18 out of 25 possible ones) 
% During the duration of 2 epochs

rats = zeros(length(in),N_EPOCHS,2*epochlength);

                         

%% For each ROI
for ii=1:length(in)

    % df/f
    iratio = in(ii).iratio;
    iratio = iratio./mean(iratio) - 1;   
    d_iratio = iratio;
    
    % Now sort those epoch ratios in matrix "rats", for each epoch-combination 
    s = in(ii).istim;
    ds = diff(s);
    
    % Find indices for each trace
    % Problem: some combinations are left over (the ones which happen AFTER
    % a double-epoch was presented)
    b_inds = find(ds)+1;
    e_inds = b_inds+2*epochlength-1; % -1 because we want to go 'till the END
                                     % of the 2nd epoch.     

    % Assign an identifying number to each epoch-combination in this ROI 
    % this is just an easier way to identify a combination using switch()
    for kk=1:length(b_inds) 
       % the last epochs might not be displayed fully, then jump out
       if(e_inds(kk) > length(s))
           continue
       end

       pairs(kk) = cantor(s(b_inds(kk)),s(e_inds(kk)));
    end
    
    % Store all examined pairs of this ROI in temp 

    num_pairs = round(length(in(1).istim)/epochlength*2);
    temp = zeros(num_pairs,N_EPOCHS,2*epochlength);
    
    for jj=1:length(pairs)
        
        % If end of data is reached, jump out of loop.
        if(jj > numel(e_inds) || e_inds(jj) > length(d_iratio))
            continue;
        end
        
        b_ind = b_inds(jj);
        e_ind = e_inds(jj);
         
        % Check for interesting combinations. 
        switch pairs(jj)
            case cantor(0,1)
                temp(jj,1,:)=d_iratio(b_ind:e_ind)';
            case cantor(0,2)
                temp(jj,2,:)=d_iratio(b_ind:e_ind)';
            case cantor(0,3)
                temp(jj,3,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,0)
                temp(jj,4,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,2)
                temp(jj,5,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,3)
                temp(jj,6,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,4)
                temp(jj,7,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,0)
                temp(jj,8,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,1)
                temp(jj,9,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,3)
                temp(jj,10,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,4)
                temp(jj,11,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,0)
                temp(jj,12,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,1)
                temp(jj,13,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,2)
                temp(jj,14,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,4)
                temp(jj,15,:)=d_iratio(b_ind:e_ind)';
            case cantor(4,1)
                temp(jj,16,:)=d_iratio(b_ind:e_ind)';    
            case cantor(4,2)
                temp(jj,17,:)=d_iratio(b_ind:e_ind)';
            case cantor(4,3)
                temp(jj,18,:)=d_iratio(b_ind:e_ind)';
            otherwise
                fprintf('No matching found for epoch combination (%d %d). Not interested? \n', cantorX(pairs(jj)),cantorY(pairs(jj)));
        end
    end
    
    temp(temp == 0) = NaN;
    rats(ii,:,:) = mean(temp,1,'omitnan');
    name{ii}=in(ii).name;
end   

% No baseline substraction
out.rats = rats;
out.neuron = [in.cell];
out.flyID = [in.flyID];
out.name = name;
