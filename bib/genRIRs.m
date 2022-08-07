
function [RIRs, rot] = genRIRs(fs, c, PosMic, PosSource, dimRoom, beta, RT60, name_dataset)
    
    Fsbb = 16e3;
    mtype = 'omnidirectional';
    order = -1; 
    dim = 3;
    orientation = 0;
    hp_filter = 1;

    [M,~]=size(PosMic);
    [N,~]=size(PosSource);
    dist_mic = abs(norm(PosMic(1,:)- PosMic(2,:)))*100; %em cm
    
    rot = zeros(5,M);
    lag = zeros(M,1);
    
    
    if(RT60==0)
        length_rir =.1*fs;
    else
        length_rir =RT60*fs;
    end
    
    RIRs = cell(N,M);
    t = zeros(M,1);
    
    for i=1:N  
        parfor j=1:M      
            
            r = PosMic(j,:);
            s = PosSource(i,:);
            
            if ~isempty(beta)
                [RIRs{i,j}, ~,  drr] = rir_generator(c, fs, r, s, dimRoom, beta*ones(1,6), length_rir, mtype, order, dim, orientation, hp_filter);
            else
                [RIRs{i,j}, ~,  drr] = rir_generator(c, fs, r, s, dimRoom, RT60, RT60*fs, mtype, order, dim, orientation, hp_filter);
            end
            t(j) = abs(norm(s - r))/c;
            
            rot(:,j)=createlabel(RT60, t(j),drr, dist_mic, fs,j,N);
            
            
%             if ~ isempty(name_dataset)
%                 h = RIRs(:,j); 
%                 save_RIR(h,RT60,fs,rot(:,j),i,j,name_dataset);
%             end
        end  
    end
   

end

function [rot] = createlabel(RT60, t,drr, dist_mic, fs,j,N)

    rot(1)=RT60;
    if N < 2
        if j~=1
            lag = 0;%round((t-t1)*fs);
        else
            t1=t;
            lag=0;
        end 
        rot(2)=lag; 
        rot(3)=dist_mic;
    else
        rot(2)=0; 
        rot(3)=0;
    end
    rot(4)=fs; 
    rot(5) = drr;
end



function save_RIR(h,RT60,fs,rot,i,j,name_dataset)
    if ~isempty(name_dataset)    
        name_Folder = name_dataset;
        path_dataset= '../../../_dataset/';
        path = fullfile(path_dataset,name_Folder);
        path_dataset_RIR =fullfile(path,'/RIRs'); 
        mkdir (path)
        mkdir (path_dataset_RIR)
        filename = ['RIR_' num2str(RT60*1000) '_' num2str(fs) '_' num2str(i-1) '_' num2str(j-1) '.mat'];
        
        save([path_dataset_RIR '/' filename],'h','rot')
    end
end
